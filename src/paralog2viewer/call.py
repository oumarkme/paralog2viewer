# SPDX-License-Identifier: GPL-3.0-or-later
"""Run Clair3 (SNV/indel) and Sawfish (SV) per haplotype.

Per-haplotype layout (written by the realign step)::

    <hap_dir>/reads.bam(.bai)  realigned reads of one haplotype
    <hap_dir>/region.bed       one 0-based half-open calling interval
    <hap_dir>/expected_cn.bed  expected copy number for Sawfish

This module adds ``clair3.vcf.gz``, ``sawfish.vcf.gz`` and ``calls.vcf.gz``
(bgzip plus tabix). The name of ``hap_dir`` is the exact ``HP`` tag value and
is used as the VCF sample name. Calls are haploid. Ported from the
CALL_CLAIR3, CALL_SAWFISH and MERGE_HAPLOTYPE_CALLS processes of the old
Nextflow pipeline.
"""

from __future__ import annotations

import collections
import functools
import hashlib
import json
import logging
import shutil
import tempfile
from pathlib import Path

import pysam

from paralog2viewer.proc import ToolError, run

# Clair3 model directory. The model is not in the image: users download it
# and bind-mount it here (README, "Clair3 model").
CLAIR3_MODEL = "/opt/clair3_model"
CLAIR3_MODEL_FILES = ("pileup.pt", "full_alignment.pt")
# SHA-256 of r1041_e82_400bps_sup_v520 (Clair3 2.0.3), the validated model.
VALIDATED_MODEL = {
    "pileup.pt": (
        "ec73e7fd3d955d8d64efbab29ba81f45c6129774a6c18beecc67c5829611fb4b"
    ),
    "full_alignment.pt": (
        "3d45d34401613df5c4dcd31f3c444e13338a5cb0c8beb395d4f6ab71682f25d3"
    ),
}

log = logging.getLogger(__name__)

BASES = frozenset("ACGTN")


@functools.cache
def model_checksums(model_dir: str = CLAIR3_MODEL) -> dict[str, str]:
    """Returns the SHA-256 of each model file and warns on an unknown model.

    Raises:
        FileNotFoundError: A model file is missing.
    """
    sums = {}
    for name in CLAIR3_MODEL_FILES:
        digest = hashlib.sha256()
        with open(Path(model_dir) / name, "rb") as handle:
            for block in iter(lambda: handle.read(1 << 20), b""):
                digest.update(block)
        sums[name] = digest.hexdigest()
    if sums != VALIDATED_MODEL:
        log.warning(
            "Clair3 model in %s is not the validated "
            "r1041_e82_400bps_sup_v520; calls may differ",
            model_dir,
        )
    return sums


# Header lines shared by calls.vcf.gz and the merged VCF.
FORMAT_LINES = (
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Haploid genotype: '
    '1 carries the allele, 0 reference, . not covered">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth for '
    'REF and ALT">',
    '##FORMAT=<ID=AF,Number=1,Type=Float,Description="Alternate allele '
    'fraction">',
    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
    '##FORMAT=<ID=SR,Number=1,Type=Integer,Description="Supporting-read '
    "count for Sawfish alleles (ALT element of the Sawfish AD field; field "
    'name not verified)">',
)
INFO_LINES = (
    '##INFO=<ID=SRC,Number=.,Type=String,Description="Caller: clair3 or '
    'sawfish">',
    '##INFO=<ID=REGION,Number=.,Type=String,Description="Calling region '
    'chrom:start-end (1-based, closed)">',
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Structural variant '
    'type">',
    '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Structural variant '
    'length">',
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of '
    'symbolic alleles">',
)
SAMPLE_FIELDS = ("DP", "AD", "AF", "GQ", "SR")


def read_bed_interval(bed: Path) -> tuple[str, int, int]:
    """Read the first interval of a BED file.

    Args:
        bed: BED file; only the first line is used.

    Returns:
        Contig, 0-based start and exclusive end.

    Raises:
        ValueError: If the interval is empty or invalid.
    """
    with open(bed, encoding="utf-8") as handle:
        fields = handle.readline().strip().split("\t")
    if len(fields) < 3:
        raise ValueError(f"Invalid callable BED line in {bed}")
    chrom, start0, end0 = fields[0], int(fields[1]), int(fields[2])
    if not 0 <= start0 < end0:
        raise ValueError("Invalid callable BED interval")
    return chrom, start0, end0


def make_header(
    contigs: dict[str, int], samples: list[str]
) -> pysam.VariantHeader:
    """Build the VCF header used for calls and merged output.

    Args:
        contigs: Contig names and lengths, in output order.
        samples: Sample names (haplotype ``HP`` values).

    Returns:
        A new header with contig, INFO, FORMAT and sample lines.
    """
    header = pysam.VariantHeader()
    for name, length in contigs.items():
        header.contigs.add(name, length=length)
    for line in INFO_LINES + FORMAT_LINES:
        header.add_line(line)
    for sample in samples:
        header.add_sample(sample)
    return header


def is_explicit(ref: str, alt: str) -> bool:
    """Return True if REF and ALT are plain sequence (not symbolic/BND)."""
    return set(ref.upper() + alt.upper()) <= BASES


def clair3_command(
    bam: Path,
    reference: Path,
    bed: Path,
    hap: str,
    out_dir: Path,
    *,
    threads: int,
    model_path: str = CLAIR3_MODEL,
) -> list[str]:
    """Return the run_clair3.sh command (options as in the Nextflow version).

    Args:
        bam: Haplotype BAM.
        reference: Reference FASTA.
        bed: Calling BED.
        hap: Sample name.
        out_dir: Clair3 output directory.
        threads: Thread count.
        model_path: Clair3 model directory.
    """
    return [
        "run_clair3.sh",
        f"--bam_fn={bam}",
        f"--ref_fn={reference}",
        f"--threads={threads}",
        "--platform=ont",
        f"--model_path={model_path}",
        f"--output={out_dir}",
        f"--bed_fn={bed}",
        f"--sample_name={hap}",
        "--snp_min_af=0.10",
        "--indel_min_af=0.13",
        "--var_pct_full=0.2",
        "--ref_pct_full=0.2",
        "--var_pct_phasing=0.8",
        "--haploid_precise",
        "--no_phasing_for_fa",
        "--chunk_size=250000000",
    ]


def sawfish_commands(
    bam: Path,
    reference: Path,
    expected_cn: Path,
    discover_dir: Path,
    joint_dir: Path,
    *,
    threads: int,
) -> list[list[str]]:
    """Return the sawfish discover and joint-call commands (as before)."""
    return [
        [
            "sawfish", "discover", "--threads", str(threads),
            "--ref", str(reference), "--bam", str(bam),
            "--expected-cn", str(expected_cn),
            "--output-dir", str(discover_dir),
            "--disable-cnv", "--cov-regex", ".",
        ],
        [
            "sawfish", "joint-call", "--threads", str(threads),
            "--sample", str(discover_dir), "--output-dir", str(joint_dir),
            "--report-supporting-reads", "--treat-single-copy-as-haploid",
        ],
    ]  # fmt: skip


def filter_clair3_vcf(
    source: Path, target: Path, bed: Path, hap: str, outside_tsv: Path
) -> int:
    """Keep Clair3 records whose start lies in the calling BED.

    Args:
        source: Clair3 ``merge_output.vcf.gz``.
        target: Output bgzip VCF (tabix-indexed).
        bed: Calling BED (0-based half-open).
        hap: Expected sample name.
        outside_tsv: Table that receives the excluded records.

    Returns:
        Number of records excluded for lying outside the BED.

    Raises:
        ValueError: On sample or contig mismatch.
    """
    chrom, start0, end0 = read_bed_interval(bed)
    excluded = 0
    with pysam.VariantFile(str(source)) as vcf:
        if list(vcf.header.samples) != [hap]:
            raise ValueError("Clair3 VCF sample identity mismatch")
        with (
            pysam.VariantFile(str(target), "wz", header=vcf.header) as out,
            open(outside_tsv, "w", encoding="utf-8") as bad,
        ):
            bad.write("CHROM\tPOS\tREF\tALT\n")
            for record in vcf:
                if record.contig != chrom:
                    raise ValueError("Clair3 VCF contig differs from BED")
                if start0 <= record.start < end0:
                    out.write(record)
                else:
                    excluded += 1
                    alts = ",".join(record.alts or [])
                    bad.write(
                        f"{record.contig}\t{record.pos}\t{record.ref}\t{alts}\n"
                    )
    pysam.tabix_index(str(target), preset="vcf", force=True)
    return excluded


def call_clair3(
    hap_dir: Path,
    reference: Path,
    *,
    threads: int = 4,
    model_path: str = CLAIR3_MODEL,
) -> Path:
    """Run Clair3 on one haplotype and write ``clair3.vcf.gz``.

    Scratch goes to a temporary directory inside ``hap_dir`` (never
    the shared parent itself) and is removed afterwards.

    Returns:
        Path of ``<hap_dir>/clair3.vcf.gz``.
    """
    hap = hap_dir.name
    scratch = Path(tempfile.mkdtemp(prefix="clair3.", dir=hap_dir))
    try:
        run(
            clair3_command(
                hap_dir / "reads.bam",
                reference,
                hap_dir / "region.bed",
                hap,
                scratch / "clair3",
                threads=threads,
                model_path=model_path,
            )
        )
        target = hap_dir / "clair3.vcf.gz"
        filter_clair3_vcf(
            scratch / "clair3" / "merge_output.vcf.gz",
            target,
            hap_dir / "region.bed",
            hap,
            hap_dir / "clair3.outside_region.tsv",
        )
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
    return target


def call_sawfish(hap_dir: Path, reference: Path, *, threads: int = 4) -> Path:
    """Run Sawfish on one haplotype and write ``sawfish.vcf.gz``.

    Returns:
        Path of ``<hap_dir>/sawfish.vcf.gz``.

    Raises:
        ToolError: If Sawfish wrote no VCF.
        ValueError: If the VCF sample is not the haplotype name.
    """
    hap = hap_dir.name
    scratch = Path(tempfile.mkdtemp(prefix="sawfish.", dir=hap_dir))
    try:
        discover, joint = sawfish_commands(
            hap_dir / "reads.bam",
            reference,
            hap_dir / "expected_cn.bed",
            scratch / "discover",
            scratch / "joint",
            threads=threads,
        )
        run(discover)
        run(joint)
        produced = scratch / "joint" / "genotyped.sv.vcf.gz"
        if not produced.exists():
            raise ToolError(f"Sawfish wrote no {produced.name}")
        target = hap_dir / "sawfish.vcf.gz"
        shutil.copyfile(produced, target)
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
    with pysam.VariantFile(str(target)) as vcf:
        if list(vcf.header.samples) != [hap]:
            raise ValueError("Sawfish VCF sample identity mismatch")
    # Re-index so a .tbi always exists whatever index Sawfish wrote.
    pysam.tabix_index(str(target), preset="vcf", force=True)
    return target


def _sample_fields(
    record: pysam.VariantRecord, hap: str, alt_index: int
) -> dict[str, object]:
    """Return DP, AD, AF, GQ (and SR for Sawfish) of one ALT allele."""
    sample = record.samples[hap]
    fields: dict[str, object] = {}
    for key in ("DP", "GQ"):
        if key in record.format and sample.get(key) is not None:
            fields[key] = sample[key]
    if "AD" in record.format and sample.get("AD") is not None:
        values = tuple(sample["AD"])
        if len(values) == len(record.alts or ()) + 1:
            fields["AD"] = (values[0], values[alt_index])
    if "AF" in record.format and sample.get("AF") is not None:
        value = sample["AF"]
        if isinstance(value, tuple):
            value = value[alt_index - 1] if len(value) >= alt_index else None
        if value is not None:
            fields["AF"] = value
    return fields


def combine_calls(
    hap_dir: Path, reference: Path, *, clair3_vcf: Path, sawfish_vcf: Path
) -> dict[str, int]:
    """Combine Clair3 and Sawfish calls of one haplotype into calls.vcf.gz.

    Keeps PASS records whose haploid GT selects one ALT, splits ALTs, checks
    REF against the FASTA, normalizes explicit alleles with
    ``bcftools norm -f ref -c e -m -any``, drops explicit alleles outside the
    calling BED, removes duplicate alleles (Clair3 wins) and keeps passing
    Sawfish symbolic and breakend alleles with their true ALT.

    Args:
        hap_dir: Haplotype directory (name is the sample name).
        reference: Reference FASTA.
        clair3_vcf: Filtered Clair3 VCF.
        sawfish_vcf: Sawfish VCF.

    Returns:
        Counts per category; also written to ``calls.audit.json``.

    Raises:
        ValueError: On sample or contig mismatch, or REF differing from
            the FASTA.
    """
    hap = hap_dir.name
    chrom, start0, end0 = read_bed_interval(hap_dir / "region.bed")
    counts: collections.Counter[str] = collections.Counter()
    raw = hap_dir / "calls.raw.vcf"
    normalized = hap_dir / "calls.normalized.vcf"
    with pysam.FastaFile(str(reference)) as fasta:
        header = make_header({chrom: fasta.get_reference_length(chrom)}, [hap])
        with pysam.VariantFile(str(raw), "w", header=header) as out:
            for caller, source in (
                ("clair3", clair3_vcf),
                ("sawfish", sawfish_vcf),
            ):
                with pysam.VariantFile(str(source)) as vcf:
                    if list(vcf.header.samples) != [hap]:
                        raise ValueError(f"{caller} sample differs from {hap}")
                    for record in vcf:
                        _collect(
                            record, caller, hap, chrom, (start0, end0),
                            fasta, out, counts,
                        )  # fmt: skip
    run(
        [
            "bcftools", "norm", "-f", reference, "-c", "e", "-m", "-any",
            "-Ov", "-o", normalized, raw,
        ]
    )  # fmt: skip
    alleles: dict[tuple[str, int, str, str], pysam.VariantRecord] = {}
    with pysam.VariantFile(str(normalized)) as source:
        header = source.header.copy()
        for record in source:
            alt = record.alts[0]
            if is_explicit(record.ref, alt) and (
                record.chrom != chrom
                or record.start < start0
                or record.stop > end0
            ):
                counts["normalized_outside_bounds"] += 1
                continue
            key = (record.chrom, record.pos, record.ref, alt)
            if key in alleles:
                counts["duplicate_allele"] += 1
                continue
            alleles[key] = record.copy()
    target = hap_dir / "calls.vcf.gz"
    with pysam.VariantFile(str(target), "wz", header=header) as out:
        for key in sorted(alleles, key=lambda k: (k[1], k[2], k[3])):
            out.write(alleles[key])
    pysam.tabix_index(str(target), preset="vcf", force=True)
    counts["alleles_written"] = len(alleles)
    raw.unlink()
    normalized.unlink()
    with open(hap_dir / "calls.audit.json", "w", encoding="utf-8") as handle:
        json.dump(
            {"hap": hap, "chrom": chrom, "start0": start0, "end0": end0,
             "counts": dict(counts)},
            handle, indent=2, sort_keys=True,
        )  # fmt: skip
        handle.write("\n")
    return dict(counts)


def _collect(
    record: pysam.VariantRecord,
    caller: str,
    hap: str,
    chrom: str,
    bounds: tuple[int, int],
    fasta: pysam.FastaFile,
    out: pysam.VariantFile,
    counts: collections.Counter[str],
) -> None:
    """Write passing haploid ALT alleles of one caller record, one each."""
    start0, end0 = bounds
    genotype = tuple(record.samples[hap].get("GT") or ())
    for alt_index, alt in enumerate(record.alts or (), 1):
        counts[f"{caller}:seen"] += 1
        if genotype != (alt_index,):
            counts[f"{caller}:non_alt_or_nonhaploid"] += 1
            continue
        if list(record.filter) != ["PASS"]:
            counts[f"{caller}:filtered"] += 1
            continue
        if record.chrom != chrom:
            raise ValueError(f"{caller} contig differs from reporting BED")
        explicit = is_explicit(record.ref, alt)
        if explicit and (record.start < start0 or record.stop > end0):
            counts[f"{caller}:outside_bounds"] += 1
            continue
        span = len(record.ref)
        if (
            explicit
            and fasta.fetch(chrom, record.start, record.start + span).upper()
            != record.ref.upper()
        ):
            raise ValueError(
                f"{caller} REF differs from reference at {chrom}:{record.pos}"
            )
        item = out.new_record(
            contig=chrom,
            start=record.start,
            stop=record.stop,
            alleles=(record.ref, alt),
        )
        item.filter.add("PASS")
        item.info["SRC"] = caller
        item.samples[hap]["GT"] = (1,)
        fields = _sample_fields(record, hap, alt_index)
        if caller == "sawfish" and "AD" in fields:
            fields["SR"] = fields["AD"][1]
        for key, value in fields.items():
            item.samples[hap][key] = value
        if not explicit:
            item.info["SVTYPE"] = str(record.info.get("SVTYPE", "SV"))
            svlen = record.info.get("SVLEN")
            if svlen is not None:
                item.info["SVLEN"] = (
                    svlen[0] if isinstance(svlen, tuple) else svlen
                )
        out.write(item)
        counts[f"{caller}:{'sequence' if explicit else 'symbolic'}"] += 1


def call_haplotype(
    hap_dir: Path,
    reference: Path,
    *,
    threads: int = 4,
    model_path: str = CLAIR3_MODEL,
) -> Path:
    """Call one haplotype with Clair3 and Sawfish and combine the calls.

    Returns:
        Path of ``<hap_dir>/calls.vcf.gz``.
    """
    clair3 = call_clair3(
        hap_dir, reference, threads=threads, model_path=model_path
    )
    sawfish = call_sawfish(hap_dir, reference, threads=threads)
    combine_calls(hap_dir, reference, clair3_vcf=clair3, sawfish_vcf=sawfish)
    return hap_dir / "calls.vcf.gz"
