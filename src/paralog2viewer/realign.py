# SPDX-License-Identifier: GPL-3.0-or-later
"""Cut windows, realign reads with minimap2 and derive calling BEDs.

Port of the former ``PREPARE_HAPLOTYPE_BAM`` process. The source reads of one
haplotype are realigned from raw sequence (Paraphase CIGARs are not trusted)
to a window around the target locus. The calling BED is the union of the
complete mapped reference spans of the realigned reads.

Work directory layout (``W`` is the run work directory, ``<HP>`` the exact
``HP`` tag value)::

    W/hap/<HP>/source_reads.bam    consumed (written by split.py)
    W/hap/<HP>/reads.bam, .bai     realigned, sorted, indexed
    W/hap/<HP>/region.bed          0-based half-open calling region
    W/hap/<HP>/expected_cn.bed     contig-wide copy number 1 for this HP
    W/hap/<HP>/alignment_audit.json  realignment losses (internal)
"""

from __future__ import annotations

import functools
import json
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pysam

from paralog2viewer.proc import pipeline, run
from paralog2viewer.split import HaplotypeJob

# Alignments within this many bases of the artificial window edge are not
# trusted to show the read's full reference span.
EDGE_MARGIN = 100


@dataclass(frozen=True)
class PreparedHaplotype:
    """Realigned reads of one haplotype.

    Attributes:
        region: Gene region key.
        hp: Exact ``HP`` tag value.
        chrom: Contig of the calling region.
        bam: Coordinate-sorted BAM (``.bai`` next to it) in full-contig
            GRCh38 coordinates.
        bed: One-line 0-based half-open calling BED.
        expected_cn: BED with copy number 1 over the whole contig.
        audit: Realignment audit JSON.
        start: 0-based start of the calling region.
        end: 0-based exclusive end of the calling region.
    """

    region: str
    hp: str
    chrom: str
    bam: Path
    bed: Path
    expected_cn: Path
    audit: Path
    start: int
    end: int


@functools.cache
def _minimap2_version() -> str:
    """Return the installed minimap2 version string (looked up once)."""
    out = subprocess.run(
        ["minimap2", "--version"], capture_output=True, text=True, check=True
    )
    return out.stdout.strip()


def prepare_haplotype(
    job: HaplotypeJob,
    reference: Path,
    work_dir: Path,
    threads: int,
) -> PreparedHaplotype | None:
    """Realign one haplotype's reads and derive its calling BED.

    The window is the target locus extended by two maximum read lengths on
    each side, clipped to the contig. Reads with any alignment within
    ``EDGE_MARGIN`` bases of the window edge are dropped whole (including
    split partners). Output coordinates are shifted back to the contig.

    Args:
        job: Haplotype from ``split_haplotypes``.
        reference: GRCh38 FASTA (indexed).
        work_dir: Run work directory ``W``; files go to ``W/hap/<HP>/``.
        threads: Threads for minimap2 and samtools sort (>= 1).

    Returns:
        The prepared haplotype, or None when no primary read overlaps the
        target locus (only the audit is written then).

    Raises:
        ValueError: The source BAM has no primary reads, or an alignment
            falls outside the window.
        ToolError: samtools or minimap2 failed.
    """
    hap_dir = Path(work_dir) / "hap" / job.hp
    hap_dir.mkdir(parents=True, exist_ok=True)
    reads_bam = hap_dir / "reads.bam"
    region_bed = hap_dir / "region.bed"
    expected_cn = hap_dir / "expected_cn.bed"
    for path in (reads_bam, hap_dir / "reads.bam.bai", region_bed, expected_cn):
        path.unlink(missing_ok=True)
    with pysam.AlignmentFile(str(job.source_bam), "rb") as source:
        if not any(not r.flag & 0x900 for r in source.fetch(until_eof=True)):
            raise ValueError(f"No primary reads in {job.source_bam}")
    with pysam.FastaFile(str(reference)) as fasta:
        contig_length = fasta.get_reference_length(job.chrom)
        # One query length on each side covers any read spanning the locus;
        # a second protects alignments with long deletions.
        t_start = max(0, job.core_start - 2 * job.max_read_length)
        t_end = min(contig_length, job.core_end + 2 * job.max_read_length)
        sequence = fasta.fetch(job.chrom, t_start, t_end)
    target_fa = hap_dir / "target.fa"
    aligned = hap_dir / "target_local.bam"
    target_fa.write_text(f">{job.chrom}\n{sequence}\n")
    try:
        pipeline(
            [
                ["samtools", "fastq", "-F", "0x900", job.source_bam],
                [
                    "minimap2",
                    "-ax",
                    "map-ont",
                    "-Y",
                    "--secondary=no",
                    "-t",
                    str(threads),
                    target_fa,
                    "-",
                ],
                [
                    "samtools",
                    "sort",
                    "-@",
                    str(threads),
                    "-o",
                    str(aligned),
                    "-",
                ],
            ],
        )
        window = {
            "chrom": job.chrom,
            "contig_length": contig_length,
            "core_start0": job.core_start,
            "core_end0": job.core_end,
            "target_start0": t_start,
            "target_end0": t_end,
            "max_read_length": job.max_read_length,
        }
        audit_info = _shift_and_filter(job, aligned, hap_dir, window)
    finally:
        target_fa.unlink(missing_ok=True)
        aligned.unlink(missing_ok=True)
    audit = {
        "sample": job.sample,
        "region": job.region,
        "hap": job.hp,
        "assembly": "GRCh38",
        "source": "full-length reads from the Paraphase BAM by primary HP",
        "aligner": f"minimap2 {_minimap2_version()} map-ont",
        **window,
        **audit_info,
    }
    audit_path = hap_dir / "alignment_audit.json"
    audit_path.write_text(json.dumps(audit, indent=2) + "\n")
    if audit["status"] != "ready":
        reads_bam.unlink(missing_ok=True)
        return None
    start, end = audit["report_start0"], audit["report_end0"]
    region_bed.write_text(f"{job.chrom}\t{start}\t{end}\n")
    expected_cn.write_text(f"{job.chrom}\t0\t{contig_length}\t{job.hp}\t1\n")
    run(["samtools", "index", "-@", str(threads), reads_bam])
    run(["samtools", "quickcheck", reads_bam])
    return PreparedHaplotype(
        region=job.region,
        hp=job.hp,
        chrom=job.chrom,
        bam=reads_bam,
        bed=region_bed,
        expected_cn=expected_cn,
        audit=audit_path,
        start=start,
        end=end,
    )


def _shift_and_filter(
    job: HaplotypeJob,
    aligned: Path,
    hap_dir: Path,
    window: dict[str, Any],
) -> dict[str, Any]:
    """Drop edge reads, shift to contig coordinates, write ``reads.bam``.

    The input is sorted in window coordinates; adding one constant offset
    keeps it sorted, so no second sort is needed.

    Returns:
        Audit fields: counts, ``report_start0``/``report_end0`` and status.
    """
    offset = window["target_start0"]
    local_end = window["target_end0"] - offset
    out_header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": job.chrom, "LN": window["contig_length"]}],
        "RG": [{"ID": job.hp, "SM": job.hp, "PL": "ONT", "DS": job.sample}],
        "PG": [
            {
                "ID": "paralog2viewer_realign",
                "PN": "minimap2",
                "VN": _minimap2_version(),
                "CL": "map-ont -Y --secondary=no",
            }
        ],
    }
    counts = {
        "mapped_primary": 0,
        "mapped_supplementary": 0,
        "unmapped": 0,
        "edge_alignments": 0,
    }
    edge_ids: set[str] = set()
    with pysam.AlignmentFile(str(aligned), "rb") as src:
        for read in src.fetch(until_eof=True):
            if not read.is_unmapped and (
                read.reference_start <= EDGE_MARGIN
                or read.reference_end >= local_end - EDGE_MARGIN
            ):
                edge_ids.add(read.query_name)
                counts["edge_alignments"] += 1
    mapped_ids: set[str] = set()
    core_ids: set[str] = set()
    min_start: int | None = None
    max_end: int | None = None
    with (
        pysam.AlignmentFile(str(aligned), "rb") as src,
        pysam.AlignmentFile(
            str(hap_dir / "reads.bam"), "wb", header=out_header
        ) as out,
    ):
        for read in src.fetch(until_eof=True):
            # A read touching the window edge lacks a trustworthy full span;
            # exclude the whole read, including its SA partners.
            if read.query_name in edge_ids:
                continue
            if read.is_unmapped:
                counts["unmapped"] += 1
                continue
            if read.reference_name != job.chrom:
                raise ValueError(
                    f"Unexpected target contig: {read.reference_name}"
                )
            if read.reference_start < 0 or read.reference_end > local_end:
                raise ValueError("Alignment exceeds target FASTA")
            read.reference_start += offset
            if read.has_tag("SA"):
                entries = []
                for entry in read.get_tag("SA").rstrip(";").split(";"):
                    fields = entry.split(",")
                    if fields[0] != job.chrom:
                        raise ValueError(
                            "Supplementary alignment has wrong contig"
                        )
                    fields[1] = str(int(fields[1]) + offset)
                    entries.append(",".join(fields))
                read.set_tag("SA", ";".join(entries) + ";")
            read.set_tag("RG", job.hp, value_type="Z")
            read.set_tag("HP", job.hp, value_type="Z")
            read.set_tag("RN", job.region, value_type="Z")
            out.write(read)
            r_start, r_end = read.reference_start, read.reference_end
            min_start = (
                r_start if min_start is None else min(min_start, r_start)
            )
            max_end = r_end if max_end is None else max(max_end, r_end)
            if read.is_supplementary:
                counts["mapped_supplementary"] += 1
            else:
                counts["mapped_primary"] += 1
                mapped_ids.add(read.query_name)
                if r_start < job.core_end and r_end > job.core_start:
                    core_ids.add(read.query_name)
    ready = bool(core_ids) and min_start is not None
    return {
        **counts,
        "mapped_primary_read_ids": len(mapped_ids),
        "edge_read_ids_excluded": len(edge_ids),
        "core_overlap_primary_read_ids": len(core_ids),
        "report_start0": min_start,
        "report_end0": max_end,
        "status": "ready" if ready else "no_mapped_core_reads",
    }
