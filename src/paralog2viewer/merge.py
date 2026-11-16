# SPDX-License-Identifier: GPL-3.0-or-later
"""Merge normalized per-haplotype calls into one VCF.

The merged VCF has one record per allele and one haploid sample column per
haplotype, named by the exact ``HP`` tag value. ``GT`` is 1 when the
haplotype carries the allele, 0 when its reads cover the REF span and do not
carry it, and ``.`` when its reads do not cover the REF span.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from pathlib import Path

import pysam

from paralog2viewer.call import (
    SAMPLE_FIELDS,
    is_explicit,
    make_header,
    read_bed_interval,
)
from paralog2viewer.datalayout import allele_key

# Minimum read depth at every REF base for a haplotype to count as covering
# a site (GT 0 instead of .). Default 1: any read spanning the site.
MIN_COVERAGE_DEPTH = 1

_BAD_NAME = re.compile(r"[\s,:;=<>\[\]]")
INFO_COPY = ("SVTYPE", "SVLEN")


def _as_tuple(value: object) -> tuple[str, ...]:
    """Return an INFO value as a tuple of strings (None gives empty)."""
    if value is None:
        return ()
    if isinstance(value, tuple):
        return tuple(str(v) for v in value)
    return (str(value),)


def read_depth(
    bam: Path,
    spans: dict[str, set[int]] | None = None,
) -> dict[str, dict[int, int]]:
    """Return per-base read depth of a BAM.

    Deletions count as covered; reference skips (spliced gaps) do not.

    Args:
        bam: Haplotype BAM (indexed or not).
        spans: Optional contig -> set of 0-based positions to keep; other
            positions are dropped to bound memory. None keeps every position.

    Returns:
        Contig -> {0-based position: depth}; only positions with reads.
    """
    depth: dict[str, dict[int, int]] = {}
    with pysam.AlignmentFile(str(bam), "rb", check_sq=False) as handle:
        for column in handle.pileup(
            stepper="nofilter",
            min_base_quality=0,
            max_depth=10_000_000,
            truncate=False,
            until_eof=True,
        ):
            if spans is not None and column.reference_pos not in spans.get(
                column.reference_name, ()
            ):
                continue
            count = sum(1 for p in column.pileups if not p.is_refskip)
            if count:
                depth.setdefault(column.reference_name, {})[
                    column.reference_pos
                ] = count
    return depth


def is_covered(
    depth: dict[str, dict[int, int]],
    chrom: str,
    start0: int,
    length: int,
    min_depth: int = MIN_COVERAGE_DEPTH,
) -> bool:
    """Return True if every base of the REF span has depth >= min_depth."""
    track = depth.get(chrom, {})
    return all(
        track.get(pos, 0) >= min_depth for pos in range(start0, start0 + length)
    )


def merge_haplotype_calls(
    hap_dirs: Sequence[Path],
    reference: Path,
    work_dir: Path,
    *,
    min_depth: int = MIN_COVERAGE_DEPTH,
) -> Path:
    """Merge ``calls.vcf.gz`` of all haplotypes into ``merged.vcf.gz``.

    Args:
        hap_dirs: Haplotype directories; each name is the exact ``HP`` value
            and must be unique. Each holds ``calls.vcf.gz``, ``reads.bam``
            and ``region.bed``. Sample columns follow this order.
        reference: Reference FASTA with ``.fai``; all its contigs go in the
            header.
        work_dir: Output directory.
        min_depth: Depth needed at every REF base for GT 0 instead of ``.``.

    Returns:
        Path of ``<work_dir>/merged.vcf.gz`` (tabix index beside it).

    Raises:
        ValueError: On duplicate or unusable HP names, a sample name that
            differs from the directory name, or a REF that differs from the
            FASTA.
    """
    names = [Path(d).name for d in hap_dirs]
    if len(set(names)) != len(names):
        raise ValueError("HP values are not unique across haplotypes")
    for name in names:
        if not name or _BAD_NAME.search(name):
            raise ValueError(f"HP value unusable as VCF sample name: {name!r}")

    # Alleles: key -> shared fields and per-haplotype FORMAT values.
    # The key carries the 0-based exclusive stop: for symbolic ALT it is the
    # record END, so calls with equal POS/ALT but different END stay apart.
    alleles: dict[tuple[str, int, str, str, int], dict] = {}
    regions: dict[str, str] = {}
    for hap_dir in map(Path, hap_dirs):
        hap = hap_dir.name
        chrom, start0, end0 = read_bed_interval(hap_dir / "region.bed")
        regions[hap] = f"{chrom}:{start0 + 1}-{end0}"
        with pysam.VariantFile(str(hap_dir / "calls.vcf.gz")) as vcf:
            if list(vcf.header.samples) != [hap]:
                raise ValueError(f"calls.vcf.gz sample differs from {hap}")
            for rec in vcf:
                alt = rec.alts[0]
                stop = (
                    rec.pos - 1 + len(rec.ref)
                    if is_explicit(rec.ref, alt)
                    else rec.stop
                )
                key = (rec.chrom, rec.pos, rec.ref, alt, stop)
                entry = alleles.setdefault(
                    key, {"src": set(), "info": {}, "haps": {}}
                )
                entry["src"].update(_as_tuple(rec.info.get("SRC")))
                for tag in INFO_COPY:
                    if tag in rec.info and tag not in entry["info"]:
                        entry["info"][tag] = rec.info[tag]
                sample = rec.samples[hap]
                entry["haps"][hap] = {
                    f: sample[f]
                    for f in SAMPLE_FIELDS
                    if f in rec.format and sample.get(f) is not None
                }

    # Depth is only needed at the REF span of each allele.
    needed: dict[str, set[int]] = {}
    for chrom, pos, ref, _alt, _stop in alleles:
        needed.setdefault(chrom, set()).update(
            range(pos - 1, pos - 1 + len(ref))
        )
    depths = {
        Path(d).name: read_depth(Path(d) / "reads.bam", needed)
        for d in hap_dirs
    }
    merged = Path(work_dir) / "merged.vcf.gz"
    Path(work_dir).mkdir(parents=True, exist_ok=True)
    with pysam.FastaFile(str(reference)) as fasta:
        contigs = dict(zip(fasta.references, fasta.lengths, strict=True))
        order = {name: i for i, name in enumerate(contigs)}
        header = make_header(contigs, names)
        ordered = sorted(alleles, key=lambda k: (order[k[0]], *k[1:]))
        with pysam.VariantFile(str(merged), "wz", header=header) as out:
            for key in ordered:
                chrom, pos, ref, alt, stop = key
                if (
                    fasta.fetch(chrom, pos - 1, pos - 1 + len(ref)).upper()
                    != ref.upper()
                ):
                    raise ValueError(
                        f"REF differs from reference for "
                        f"{allele_key(chrom, pos, ref, alt)}"
                    )
                entry = alleles[key]
                rec = out.new_record(
                    contig=chrom,
                    start=pos - 1,
                    stop=stop,
                    alleles=(ref, alt),
                )
                rec.info["SRC"] = tuple(sorted(entry["src"]))
                hap_regions = sorted({regions[h] for h in entry["haps"]})
                rec.info["REGION"] = tuple(hap_regions)
                for tag, value in entry["info"].items():
                    rec.info[tag] = value
                if not is_explicit(ref, alt) and "SVTYPE" not in rec.info:
                    rec.info["SVTYPE"] = "SV"
                for hap in names:
                    sample = rec.samples[hap]
                    if hap in entry["haps"]:
                        sample["GT"] = (1,)
                        for field, value in entry["haps"][hap].items():
                            sample[field] = value
                    elif is_covered(depths[hap], chrom, pos - 1, len(ref),
                                    min_depth):  # fmt: skip
                        sample["GT"] = (0,)
                    else:
                        sample["GT"] = (None,)
                out.write(rec)
    pysam.tabix_index(str(merged), preset="vcf", force=True)
    return merged
