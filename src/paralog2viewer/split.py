# SPDX-License-Identifier: GPL-3.0-or-later
"""Select confident primary Paraphase reads per gene region and haplotype.

Port of the former ``SPLIT_HAPLOTYPES`` process. Reads are selected by read
ID directly from the Paraphase BAM (confident primary ``RN``, ``HP`` and
``YC`` tags) and written to one ``source_reads.bam`` per haplotype.

Work directory layout written here (``W`` is the run work directory)::

    W/split_audit.json
    W/hap/<HP>/source_reads.bam
"""

from __future__ import annotations

import collections
import csv
import json
import re
from collections.abc import Mapping
from dataclasses import dataclass
from importlib import resources
from pathlib import Path
from typing import Any

import pysam

# Paraphase YC colors that mark a confident haplotype assignment.
CONFIDENT_COLORS = frozenset({"166,206,227", "178,223,138", "177,156,217"})
_SAFE_LABEL = re.compile(r"[A-Za-z0-9_]+")
_PHASE_REGION = re.compile(r"38:([^:,]+):(\d+)-(\d+)")

NativeRegions = Mapping[str, tuple[str, int, int]]


@dataclass(frozen=True)
class HaplotypeJob:
    """One resolved haplotype ready for realignment.

    Attributes:
        region: Gene region key in the Paraphase JSON (``RN`` tag value).
        hp: Exact ``HP`` tag value; unique across all regions of a sample.
        chrom: GRCh38 contig of the target locus.
        core_start: 0-based start of the target locus (half-open).
        core_end: 0-based exclusive end of the target locus.
        max_read_length: Longest selected primary read, in bases.
        source_bam: BAM with the selected Paraphase records.
        sample: Sample label written to the read-group ``DS`` field.
    """

    region: str
    hp: str
    chrom: str
    core_start: int
    core_end: int
    max_read_length: int
    source_bam: Path
    sample: str = ""


def load_native_regions(path: Path | None = None) -> dict[str, Any]:
    """Read the native (gene2) locus table.

    Args:
        path: Tab-separated file with columns ``region``, ``chrom``,
            ``start`` and ``end``; start and end are 1-based and inclusive,
            as in the Paraphase-era table. Defaults to the packaged
            ``resources/native_regions.tsv``.

    Returns:
        Mapping from region key to ``(chrom, start, end)``.
    """
    if path is None:
        text = (
            resources.files("paralog2viewer")
            .joinpath("resources", "native_regions.tsv")
            .read_text()
        )
    else:
        text = Path(path).read_text()
    rows = csv.DictReader(
        (line for line in text.splitlines() if not line.startswith("#")),
        delimiter="\t",
    )
    return {
        row["region"]: (row["chrom"], int(row["start"]), int(row["end"]))
        for row in rows
    }


def parse_phase_region(text: str) -> list[tuple[str, int, int]]:
    """Parse a Paraphase ``phase_region`` string (GRCh38 only).

    Args:
        text: ``38:<chrom>:<start>-<end>``, 1-based inclusive; clustered
            genes list several loci separated by commas.

    Returns:
        ``(chrom, start, end)`` per locus, 1-based inclusive, in input order.

    Raises:
        ValueError: If any locus is not GRCh38 ``38:chrom:start-end``.
    """
    loci = []
    for item in text.split(","):
        match = _PHASE_REGION.fullmatch(item.strip())
        if not match:
            raise ValueError(f"unexpected phase_region {text!r}")
        loci.append((match.group(1), int(match.group(2)), int(match.group(3))))
    return loci


def dictionary_mismatches(
    header: dict[str, Any], ref: pysam.FastaFile
) -> list[str] | None:
    """Compare BAM @SQ lines with a FASTA.

    Args:
        header: BAM header as a dict (``AlignmentFile.header.to_dict()``).
        ref: Open reference FASTA.

    Returns:
        None if every @SQ contig exists in the FASTA with equal length;
        otherwise the names of the offending contigs (empty list when the
        header has no @SQ lines).
    """
    sq = header.get("SQ", [])
    bad = [
        row["SN"]
        for row in sq
        if row["SN"] not in ref.references
        or row["LN"] != ref.get_reference_length(row["SN"])
    ]
    return bad if bad or not sq else None


def _check_dictionary(header: dict[str, Any], ref: pysam.FastaFile) -> None:
    """Raise ValueError unless the BAM @SQ lines match the FASTA."""
    if dictionary_mismatches(header, ref) is not None:
        raise ValueError(
            "BAM sequence dictionary does not match the configured GRCh38 FASTA"
        )


def _is_confident(read: pysam.AlignedSegment) -> bool:
    """True for a mapped primary record with RN, HP and a confident YC."""
    return (
        not read.is_unmapped
        and not read.is_secondary
        and not read.is_supplementary
        and read.has_tag("RN")
        and read.has_tag("HP")
        and read.has_tag("YC")
        and read.get_tag("YC") in CONFIDENT_COLORS
    )


def _is_gene2(
    region: str, hap: str, smn2: set[str], native: NativeRegions
) -> bool:
    """True when the haplotype is the paralog/pseudogene copy."""
    return (
        (region == "pms2" and "pms2cl" in hap)
        or (region == "strc" and "strcp1" in hap)
        or (region in {"ikbkg", "ncf1"} and "pseudo" in hap)
        or (region == "smn1" and hap in smn2)
        or (region in native and "_gene2hap" in hap)
    )


def _native_target(
    region: str,
    hap: str,
    details: dict[str, Any],
    bam: pysam.AlignmentFile,
    native: NativeRegions,
) -> tuple[tuple[str, int, int, Any] | None, str]:
    """Resolve the native (gene2) target of one haplotype.

    Returns:
        ``((chrom, start0, end0, boundary), "used")`` when the JSON has a
        complete native boundary and a confident primary read of this
        haplotype lies in the native region (0-based half-open); otherwise
        ``(None, status)`` and the caller keeps the primary locus.
    """
    boundary = details.get("boundary_gene2")
    if not boundary or None in boundary:
        return None, "missing_native_boundary"
    chrom, lo, hi = native[region]
    for read in bam.fetch(chrom, lo - 1, hi):
        if (
            _is_confident(read)
            and read.get_tag("RN") == region
            and read.get_tag("HP") == hap
        ):
            return (chrom, lo - 1, hi, boundary), "used"
    return None, "no_native_primary_reads"


def _plan_jobs(
    data: dict[str, Any],
    bam_path: Path,
    ref: pysam.FastaFile,
    assignments: dict[tuple[str, str], set[str]],
    native: NativeRegions,
    audit: dict[str, Any],
    work_dir: Path,
    sample: str,
) -> list[dict[str, Any]]:
    """Decide the target locus of every resolved haplotype (job dicts)."""
    jobs: list[dict[str, Any]] = []
    # Separate handle: indexed fetch must not share the streaming handle.
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for region, result in sorted(data.items()):
            final = result.get("final_haplotypes") or {}
            audit["genes"].append(
                {
                    "gene": region,
                    "resolved_haplotypes": len(set(final.values())),
                    "status": "resolved" if final else "no_resolved_haplotypes",
                }
            )
            if not final:
                continue
            try:
                loci = parse_phase_region(result["phase_region"])
            except ValueError as err:
                raise ValueError(
                    f"Unexpected assembly or phase region: {region}"
                ) from err
            # Several loci must share one contig; the job window is their
            # covering span (reads are chosen by read ID, not by window).
            if len({c for c, _lo, _hi in loci}) != 1:
                raise ValueError(
                    f"phase_region spans several contigs: {region}"
                )
            chrom = loci[0][0]
            lo = min(locus[1] for locus in loci)
            hi = max(locus[2] for locus in loci)
            if chrom not in ref.references or not (
                0 < lo < hi <= ref.get_reference_length(chrom)
            ):
                raise ValueError(f"Invalid phase coordinates: {region}")
            smn2 = set((result.get("smn2_haplotypes") or {}).values())
            smn2.update((result.get("smn_del78_haplotypes") or {}).values())
            for hap in sorted(set(final.values())):
                if not _SAFE_LABEL.fullmatch(
                    region
                ) or not _SAFE_LABEL.fullmatch(hap):
                    raise ValueError(
                        f"Unsafe region or haplotype label: {region}:{hap}"
                    )
                selected, fetch_lo, fetch_hi = chrom, lo - 1, hi
                details = result["haplotype_details"][hap]
                bounds = details.get("boundary")
                status = "not_applicable"
                if _is_gene2(region, hap, smn2, native):
                    target, status = _native_target(
                        region, hap, details, bam, native
                    )
                    if target:
                        selected, fetch_lo, fetch_hi, bounds = target
                if (
                    not 0
                    <= fetch_lo
                    < fetch_hi
                    <= ref.get_reference_length(selected)
                ):
                    raise ValueError(
                        f"Invalid GRCh38 target locus: {region}:{hap}"
                    )
                directory = work_dir / "hap" / hap
                directory.mkdir(parents=True, exist_ok=True)
                jobs.append(
                    {
                        "sample": sample,
                        "region": region,
                        "hap": hap,
                        "projection": "gene2"
                        if status == "used"
                        else "primary",
                        "native_projection_status": status,
                        "chrom": selected,
                        "core_start0": fetch_lo,
                        "core_end0": fetch_hi,
                        "paraphase_boundary": bounds,
                        "assigned_read_ids": len(assignments[(region, hap)]),
                        "source_bam": str(
                            (directory / "source_reads.bam").resolve()
                        ),
                    }
                )
    return jobs


def _read_membership(
    jobs: list[dict[str, Any]],
    assignments: dict[tuple[str, str], set[str]],
) -> dict[str, list[int]]:
    """Map read ID to the indexes of the jobs that select it.

    Raises:
        ValueError: An HP value is used by two regions, or one read has
            conflicting primary HP assignments within a region.
    """
    owners: dict[str, str] = {}
    for job in jobs:
        other = owners.setdefault(job["hap"], job["region"])
        if other != job["region"]:
            raise ValueError(
                f"HP value {job['hap']} used by regions {other} and "
                f"{job['region']}; Genotype must be unique"
            )
    membership: dict[str, list[int]] = collections.defaultdict(list)
    unique: dict[tuple[str, str], str] = {}
    for index, job in enumerate(jobs):
        for read_id in assignments[(job["region"], job["hap"])]:
            previous = unique.setdefault((job["region"], read_id), job["hap"])
            if previous != job["hap"]:
                raise ValueError(
                    f"Conflicting primary HP assignments: "
                    f"{job['region']}:{read_id}"
                )
            membership[read_id].append(index)
    return membership


def _extract_reads(
    bam: Path,
    header: dict[str, Any],
    jobs: list[dict[str, Any]],
    membership: dict[str, list[int]],
    sample: str,
    threads: int,
) -> int:
    """Write every job's ``source_reads.bam`` and fill its count fields.

    Returns:
        Number of records scanned in the Paraphase BAM.
    """
    writers: dict[int, pysam.AlignmentFile] = {}
    primary_ids: dict[int, set[str]] = collections.defaultdict(set)
    counts: dict[int, collections.Counter[str]] = collections.defaultdict(
        collections.Counter
    )
    max_len: dict[int, int] = collections.defaultdict(int)
    for index, job in enumerate(jobs):
        if not job["assigned_read_ids"]:
            continue
        out_header = dict(header)
        out_header["RG"] = [
            {"ID": job["hap"], "SM": job["hap"], "PL": "ONT", "DS": sample}
        ]
        out_header["PG"] = list(header.get("PG", [])) + [
            {"ID": "paralog2viewer_extract", "PN": "paralog2viewer"}
        ]
        writers[index] = pysam.AlignmentFile(
            job["source_bam"], "wb", header=out_header
        )
    scanned = 0
    try:
        # Fresh handle for this second full pass: reusing the handle of the
        # assignment scan truncates iteration on Paraphase BAMs.
        with pysam.AlignmentFile(str(bam), "rb", threads=threads) as extract:
            for read in extract.fetch(until_eof=True):
                scanned += 1
                if read.is_secondary or read.query_name not in membership:
                    continue
                for index in membership[read.query_name]:
                    # A read chimeric across paralogous loci can carry several
                    # non-supplementary records, each with the full sequence.
                    # Keep only the first so FASTQ extraction sees it once.
                    if (
                        not read.is_supplementary
                        and read.query_name in primary_ids[index]
                    ):
                        continue
                    job = jobs[index]
                    read.set_tag("HP", job["hap"], value_type="Z")
                    read.set_tag("RN", job["region"], value_type="Z")
                    read.set_tag("RG", job["hap"], value_type="Z")
                    writers[index].write(read)
                    if read.is_supplementary:
                        counts[index]["supplementary"] += 1
                    else:
                        counts[index]["primary"] += 1
                        primary_ids[index].add(read.query_name)
                        max_len[index] = max(
                            max_len[index], read.query_length or 0
                        )
    finally:
        for writer in writers.values():
            writer.close()
    for index, job in enumerate(jobs):
        job["primary_records"] = counts[index]["primary"]
        job["supplementary_records"] = counts[index]["supplementary"]
        job["matched_primary_read_ids"] = len(primary_ids[index])
        job["max_read_length"] = max_len[index]
        ready = bool(primary_ids[index] and max_len[index])
        job["status"] = "ready" if ready else "no_matching_primary"
    return scanned


def split_haplotypes(
    bam: Path,
    paraphase_json: Path,
    reference: Path,
    work_dir: Path,
    threads: int,
    *,
    sample: str = "",
    native_regions: NativeRegions | None = None,
) -> list[HaplotypeJob]:
    """Write one source BAM per resolved haplotype.

    Every gene in the Paraphase JSON is processed. Reads are chosen by read
    ID from confident primary ``RN``/``HP``/``YC`` records; all non-secondary
    records of a chosen read (including supplementary ones) are copied with
    ``HP``, ``RN`` and ``RG`` set to the haplotype. Gene2 haplotypes use the
    native locus only when the JSON has a complete native boundary and
    confident primary reads lie in the native region.

    Args:
        bam: Paraphase BAM with its ``.bai`` next to it.
        paraphase_json: Paraphase JSON for the same sample.
        reference: GRCh38 FASTA (indexed); must match the BAM @SQ lines.
        work_dir: Run work directory ``W``; writes ``W/split_audit.json`` and
            ``W/hap/<HP>/source_reads.bam``.
        threads: Decompression threads for reading the BAM (>= 1).
        sample: Sample label for the read-group ``DS`` field.
        native_regions: Gene2 native loci, region -> (chrom, 1-based start,
            1-based inclusive end); defaults to the packaged table.

    Returns:
        Jobs with at least one matched primary read. ``core_start`` is
        0-based and ``core_end`` exclusive (half-open).

    Raises:
        ValueError: Dictionary mismatch, unsafe labels, invalid coordinates,
            conflicting HP assignments, or an HP value used by two regions.
    """
    bam, work_dir = Path(bam), Path(work_dir)
    native = load_native_regions() if native_regions is None else native_regions
    data = json.loads(Path(paraphase_json).read_text())
    work_dir.mkdir(parents=True, exist_ok=True)
    audit: dict[str, Any] = {
        "sample": sample,
        "assembly": "GRCh38",
        "source_bam": str(bam),
        "source_json": str(paraphase_json),
        "reference": str(reference),
        "assignment_rule": "confident primary Paraphase RN+HP+YC; reads and "
        "full sequence taken directly from the Paraphase BAM",
        "genes": [],
        "haplotypes": [],
    }
    assignments: dict[tuple[str, str], set[str]] = collections.defaultdict(set)
    with (
        pysam.FastaFile(str(reference)) as ref,
        pysam.AlignmentFile(str(bam), "rb", threads=threads) as handle,
    ):
        header = handle.header.to_dict()
        _check_dictionary(header, ref)
        for read in handle.fetch(until_eof=True):
            if _is_confident(read):
                assignments[(read.get_tag("RN"), read.get_tag("HP"))].add(
                    read.query_name
                )
        jobs = _plan_jobs(
            data, bam, ref, assignments, native, audit, work_dir, sample
        )
    audit["haplotypes"] = jobs
    membership = _read_membership(jobs, assignments)
    audit["paraphase_bam_records_scanned"] = _extract_reads(
        bam, header, jobs, membership, sample, threads
    )
    ready: list[HaplotypeJob] = []
    for job in jobs:
        if job["assigned_read_ids"]:
            pysam.quickcheck(job["source_bam"])
        if job["status"] == "ready":
            ready.append(
                HaplotypeJob(
                    region=job["region"],
                    hp=job["hap"],
                    chrom=job["chrom"],
                    core_start=job["core_start0"],
                    core_end=job["core_end0"],
                    max_read_length=job["max_read_length"],
                    source_bam=Path(job["source_bam"]),
                    sample=sample,
                )
            )
    (work_dir / "split_audit.json").write_text(
        json.dumps(audit, indent=2) + "\n"
    )
    return ready
