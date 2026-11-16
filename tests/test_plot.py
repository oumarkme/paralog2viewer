# SPDX-License-Identifier: GPL-3.0-or-later
"""Tests for paralog2viewer.plot with synthetic BAM/VCF data."""

from __future__ import annotations

import random
from pathlib import Path

import pysam
import pytest
from matplotlib.axes import Axes
from PIL import Image

from paralog2viewer import plot
from paralog2viewer.datalayout import allele_key

random.seed(7)
REF_SEQ = "".join(random.choice("ACGT") for _ in range(3000))
KNOWN_HEADER = (
    "functional_gene\tparalog\ttranscript\thgvs\tlabel\tchrom\tpos\tref"
    "\talt\tevidence_type\treference\n"
)


def _alt_base(pos1: int) -> str:
    return "A" if REF_SEQ[pos1 - 1] != "A" else "C"


def _write_bam(
    folder: Path,
    start: int,
    n_reads: int,
    alt_pos: tuple[int, ...],
    length: int = 600,
) -> None:
    """Write n reads of `length` from `start` (0-based), alt at alt_pos."""
    folder.mkdir(parents=True)
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": len(REF_SEQ)}],
    }
    path = folder / "reads.bam"
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for i in range(n_reads):
            s = start + i * 10
            seq = list(REF_SEQ[s : s + length])
            for p in alt_pos:
                if s <= p - 1 < s + length:
                    seq[p - 1 - s] = _alt_base(p)
            read = pysam.AlignedSegment()
            read.query_name = f"{folder.name}_r{i}"
            read.query_sequence = "".join(seq)
            read.flag = 0
            read.reference_id = 0
            read.reference_start = s
            read.mapping_quality = 60
            read.cigarstring = f"{length}M"
            read.query_qualities = pysam.qualitystring_to_array("I" * length)
            out.write(read)
    pysam.index(str(path))


def _write_vcf(
    path: Path, samples: list[str], records: list[tuple[int, list[str]]]
) -> None:
    """records: (pos, GT per sample); ALT is a SNV unless pos is negative."""
    lines = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=chr1,length=3000>",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(samples),
    ]
    for pos, gts in records:
        lines.append(
            f"chr1\t{pos}\t.\t{REF_SEQ[pos - 1]}\t{_alt_base(pos)}\t.\tPASS"
            "\t.\tGT\t" + "\t".join(gts)
        )
    raw = path.with_suffix("")
    raw.write_text("\n".join(lines) + "\n")
    pysam.tabix_index(str(raw), preset="vcf", force=True, keep_original=False)


@pytest.fixture()
def data(tmp_path: Path) -> dict:
    fa = tmp_path / "ref.fa"
    fa.write_text(f">chr1\n{REF_SEQ}\n")
    pysam.faidx(str(fa))
    known = tmp_path / "known.tsv"
    pos = 900
    known.write_text(
        "# comment\n"
        + KNOWN_HEADER
        + f"G\tP\tNM_1\tc.1A>C\tG c.1A>C\tchr1\t{pos}\t{REF_SEQ[pos - 1]}"
        f"\t{_alt_base(pos)}\tparalog-specific site\tPMID:1\n"
    )
    # R_hap1: ref. R_hap2: alt at 500 (P/LP), 600 (shared with hap1 as
    # well, so not drawn) and 900 (pseudogene-type). Reads span 100-1090.
    _write_bam(tmp_path / "h1", 100, 40, (600,))
    _write_bam(tmp_path / "h2", 100, 40, (500, 600, 900))
    # R2: two haplotypes at different loci.
    _write_bam(tmp_path / "g1", 100, 10, ())
    _write_bam(tmp_path / "g2", 1800, 10, ())
    vcf = tmp_path / "m.vcf.gz"
    samples = ["R_hap1", "R_hap2", "R2_gene1hap1", "R2_gene2hap1"]
    _write_vcf(
        vcf,
        samples,
        [
            (500, ["0", "1", ".", "."]),
            (600, ["1", "1", ".", "."]),
            (900, [".", "1", ".", "."]),
        ],
    )
    return {
        "vcf": vcf,
        "fa": fa,
        "known": known,
        "dirs": {
            "R_hap1": tmp_path / "h1",
            "R_hap2": tmp_path / "h2",
            "R2_gene1hap1": tmp_path / "g1",
            "R2_gene2hap1": tmp_path / "g2",
        },
        "keys": {allele_key("chr1", 500, REF_SEQ[499], _alt_base(500))},
        "out": tmp_path / "figures",
    }


def test_figures_named_per_region_and_locus(data: dict) -> None:
    paths = plot.plot_regions(
        data["vcf"],
        data["dirs"],
        data["keys"],
        data["known"],
        data["fa"],
        data["out"],
        workers=2,
    )
    assert [p.name for p in paths] == [
        "R.jpeg",
        "R2_locus1.jpeg",
        "R2_locus2.jpeg",
    ]
    with Image.open(paths[0]) as img:
        assert img.format == "JPEG"
        assert img.info["dpi"][0] == pytest.approx(300, abs=1)
        assert abs(img.size[0] - 180 / 25.4 * 300) <= 1
        assert img.size[1] / 300 * 25.4 <= 247.0


def test_marks_follow_the_rules(
    data: dict, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls = []
    original = Axes.scatter

    def spy(self, x, y, *args, **kwargs):  # noqa: ANN001
        calls.append((kwargs["marker"], sorted(set(x)), len(x)))
        return original(self, x, y, *args, **kwargs)

    monkeypatch.setattr(Axes, "scatter", spy)
    plot.plot_regions(
        data["vcf"],
        {k: v for k, v in data["dirs"].items() if k.startswith("R_")},
        data["keys"],
        data["known"],
        data["fa"],
        data["out"],
        workers=1,
    )
    by_marker = {m: (xs, n) for m, xs, n in calls}
    # Dots for 500 and 900 (hap2 only); 600 is shared by all: absent.
    assert by_marker["o"][0] == [500, 900]
    assert 600 not in by_marker["o"][0]
    assert by_marker["D"] == ([500], 40)  # P/LP: all 40 hap2 reads
    assert by_marker["s"][0] == [900]  # pseudogene square


def test_alt_state_snv_deletion_insertion_and_no_cover() -> None:
    def make(seq: str, cigar: str, start: int = 0):
        r = pysam.AlignedSegment()
        r.query_sequence = seq
        r.reference_start = start
        r.cigarstring = cigar
        return r

    r = make("ACGTACGT", "8M")
    m = {i: i for i in range(8)}
    assert plot.alt_state(r, m, 2, "G", "G") is True
    assert plot.alt_state(r, m, 2, "G", "T") is False
    assert plot.alt_state(r, m, 20, "G", "T") is None
    # Deletion of "CG" after anchor A at 0: query ACGT -> A + T.
    d = make("ATAC", "1M2D3M")
    md = {0: 0, 3: 1, 4: 2, 5: 3}
    assert plot.alt_state(d, md, 0, "ACG", "A") is True
    assert plot.alt_state(r, m, 0, "ACG", "A") is False
    # Insertion "TT" after anchor at 0.
    ins = make("ATTCG", "1M2I2M")
    mi = {0: 0, 1: 3, 2: 4}
    assert plot.alt_state(ins, mi, 0, "A", "ATT") is True
    assert plot.alt_state(r, m, 0, "A", "ATT") is False


def test_alt_state_delins_compares_replacement_bases() -> None:
    """GAGT>GC (PMS2-like) matches the read with C, not a plain deletion."""

    def make(seq: str, cigar: str):
        r = pysam.AlignedSegment()
        r.query_sequence = seq
        r.reference_start = 0
        r.cigarstring = cigar
        return r

    # Reference: G A G T T T T. Read with delins: G C then T T T.
    delins = make("GCTTT", "1M1I3D3M")
    m_delins = {0: 0, 4: 2, 5: 3, 6: 4}
    # Read lacking only "AG": G T T T T.
    plain = make("GTTTT", "1M2D4M")
    m_plain = {0: 0, 3: 1, 4: 2, 5: 3, 6: 4}
    # Reference-like read.
    wild = make("GAGTTTT", "7M")
    m_wild = {i: i for i in range(7)}
    assert plot.alt_state(delins, m_delins, 0, "GAGT", "GC") is True
    assert plot.alt_state(plain, m_plain, 0, "GAGT", "GC") is False
    assert plot.alt_state(wild, m_wild, 0, "GAGT", "GC") is False
    # A plain deletion of the same span is still recognised.
    assert plot.alt_state(plain, m_plain, 0, "GAG", "G") is True


def test_sv_state_long_deletion_and_split_alignment() -> None:
    def make(cigar: str, start: int):
        r = pysam.AlignedSegment()
        r.reference_start = start
        r.cigarstring = cigar
        r.query_sequence = "A" * 100
        return r

    deleted = make("50M60D50M", 1000)
    assert plot.sv_state([deleted], 1060, 1100) is True
    plain = make("100M", 1000)
    assert plot.sv_state([plain], 1050, 1051) is False
    assert plot.sv_state([plain], 5000, 5001) is None
    assert plot.sv_state([plain, make("100M", 3000)], 3050, 3051) is True


def test_cluster_loci_splits_non_overlapping_and_joins_overlaps() -> None:
    loci = plot.cluster_loci(
        {
            "a": ("chr1", 0, 100),
            "b": ("chr1", 50, 200),
            "c": ("chr1", 300, 400),
            "d": ("chr2", 0, 10),
        }
    )
    assert [lc["haps"] for lc in loci] == [["a", "b"], ["c"], ["d"]]
    assert loci[0]["end0"] == 200


def test_missing_contig_and_missing_bam_raise(data: dict) -> None:
    other = data["fa"].parent / "other.fa"
    other.write_text(">chrX\nACGT\n")
    pysam.faidx(str(other))
    with pytest.raises(ValueError, match="contigs"):
        plot.plot_regions(
            data["vcf"], data["dirs"], set(), data["known"], other, data["out"]
        )
    bad = {"R_hap1": data["fa"].parent / "nope"}
    with pytest.raises(FileNotFoundError):
        plot.plot_regions(
            data["vcf"], bad, set(), data["known"], data["fa"], data["out"]
        )


def test_single_haplotype_draws_no_colored_dots(
    data: dict, monkeypatch: pytest.MonkeyPatch
) -> None:
    markers = []
    original = Axes.scatter

    def spy(self, x, y, *args, **kwargs):  # noqa: ANN001
        markers.append(kwargs["marker"])
        return original(self, x, y, *args, **kwargs)

    monkeypatch.setattr(Axes, "scatter", spy)
    plot.plot_regions(
        data["vcf"],
        {"R_hap2": data["dirs"]["R_hap2"]},
        set(),
        data["known"],
        data["fa"],
        data["out"],
    )
    assert "o" not in markers
    assert markers == ["s"]
