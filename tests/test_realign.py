# SPDX-License-Identifier: GPL-3.0-or-later
import json
import random
import shutil
from pathlib import Path

import pysam
import pytest

from paralog2viewer import realign
from paralog2viewer.realign import prepare_haplotype
from paralog2viewer.split import HaplotypeJob

pytestmark = pytest.mark.skipif(
    shutil.which("samtools") is None or shutil.which("minimap2") is None,
    reason="samtools and minimap2 required",
)

CONTIG = 9000
COMP = str.maketrans("ACGT", "TGCA")


@pytest.fixture
def ref(tmp_path):
    rng = random.Random(7)
    seq = "".join(rng.choice("ACGT") for _ in range(CONTIG))
    path = tmp_path / "ref.fa"
    path.write_text(f">chrT\n{seq}\n")
    pysam.faidx(str(path))
    return path


def ref_seq(ref: Path) -> str:
    with pysam.FastaFile(str(ref)) as fasta:
        return fasta.fetch("chrT")


def make_source(path: Path, ref: Path, reads: list[tuple]) -> Path:
    """Reads are (name, ref_start, ref_end, reverse, flag).

    The stored SEQ is the reference substring; Paraphase-style placement is
    deliberately wrong (position 0) to show that coordinates are not reused.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    seq = ref_seq(ref)
    header = {
        "SQ": [{"SN": "chrT", "LN": CONTIG}],
        "RG": [{"ID": "hp1", "SM": "hp1"}],
    }
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for name, start, end, reverse, flag in reads:
            read = pysam.AlignedSegment(out.header)
            read.query_name = name
            read.flag = flag | (16 if reverse else 0)
            read.reference_id = 0
            read.reference_start = 0
            read.mapping_quality = 60
            read.query_sequence = seq[start:end]
            read.query_qualities = pysam.qualitystring_to_array(
                "I" * (end - start)
            )
            read.cigartuples = [(0, end - start)]
            read.set_tag("RN", "g")
            read.set_tag("HP", "hp1")
            out.write(read)
    return path


def make_job(source: Path, core=(5000, 5800), sample="S") -> HaplotypeJob:
    return HaplotypeJob(
        region="g",
        hp="hp1",
        chrom="chrT",
        core_start=core[0],
        core_end=core[1],
        max_read_length=1500,
        source_bam=source,
        sample=sample,
    )


def test_realigns_and_derives_bed_from_read_spans(tmp_path, ref):
    reads = [
        ("a", 5000, 5800, False, 0),
        ("b", 5200, 5900, True, 0),  # reverse-strand record
        ("c", 5000, 5100, False, 2048),  # supplementary record: not realigned
    ]
    source = make_source(
        tmp_path / "hap" / "hp1" / "source_reads.bam", ref, reads
    )
    work = tmp_path
    got = prepare_haplotype(make_job(source), ref, work, 2)
    hap = work / "hap" / "hp1"
    assert got is not None
    assert (got.start, got.end) == (5000, 5900)
    assert got.bed.read_text() == "chrT\t5000\t5900\n"
    assert got.expected_cn.read_text() == f"chrT\t0\t{CONTIG}\thp1\t1\n"
    assert got.bam == hap / "reads.bam"
    assert (hap / "reads.bam.bai").exists()
    assert sorted(p.name for p in hap.iterdir()) == [
        "alignment_audit.json",
        "expected_cn.bed",
        "reads.bam",
        "reads.bam.bai",
        "region.bed",
        "source_reads.bam",
    ]
    with pysam.AlignmentFile(str(got.bam), "rb") as bam:
        assert bam.header.to_dict()["HD"]["SO"] == "coordinate"
        records = {r.query_name: r for r in bam}
    assert set(records) == {"a", "b"}
    assert records["a"].reference_start == 5000
    assert records["b"].reference_start == 5200
    assert records["b"].is_reverse
    assert records["b"].reference_end == 5900
    for rec in records.values():
        assert (rec.get_tag("HP"), rec.get_tag("RG")) == ("hp1", "hp1")
        assert rec.get_tag("RN") == "g"
    audit = json.loads(got.audit.read_text())
    assert audit["status"] == "ready"
    assert audit["mapped_primary"] == 2
    assert audit["target_start0"] == 2000
    assert audit["sample"] == "S"


def test_read_touching_window_edge_is_dropped_whole(tmp_path, ref):
    # Window starts at 2000; a read at 2050 lies within 100 bp of the edge.
    reads = [("a", 5000, 5800, False, 0), ("edge", 2050, 3200, False, 0)]
    source = make_source(tmp_path / "s.bam", ref, reads)
    got = prepare_haplotype(make_job(source), ref, tmp_path, 1)
    assert got is not None
    assert (got.start, got.end) == (5000, 5800)
    audit = json.loads(got.audit.read_text())
    assert audit["edge_read_ids_excluded"] == 1
    with pysam.AlignmentFile(str(got.bam), "rb") as bam:
        assert [r.query_name for r in bam] == ["a"]


def test_no_core_overlap_returns_none_and_no_bed(tmp_path, ref):
    source = make_source(tmp_path / "s.bam", ref, [("a", 3000, 3500, False, 0)])
    got = prepare_haplotype(make_job(source), ref, tmp_path, 1)
    hap = tmp_path / "hap" / "hp1"
    assert got is None
    assert not (hap / "region.bed").exists()
    assert not (hap / "reads.bam").exists()
    audit = json.loads((hap / "alignment_audit.json").read_text())
    assert audit["status"] == "no_mapped_core_reads"
    assert audit["report_start0"] == 3000


def test_source_without_primary_reads_raises(tmp_path, ref):
    source = make_source(
        tmp_path / "s.bam", ref, [("a", 5000, 5800, False, 2048)]
    )
    with pytest.raises(ValueError, match="No primary reads"):
        prepare_haplotype(make_job(source), ref, tmp_path, 1)


def test_minimap2_version_is_looked_up_once(monkeypatch) -> None:
    calls = []

    def fake_run(cmd, **kwargs):
        calls.append(cmd)
        return realign.subprocess.CompletedProcess(cmd, 0, "2.30-r1\n", "")

    realign._minimap2_version.cache_clear()
    monkeypatch.setattr(realign.subprocess, "run", fake_run)
    assert realign._minimap2_version() == realign._minimap2_version()
    assert len(calls) == 1
    realign._minimap2_version.cache_clear()
