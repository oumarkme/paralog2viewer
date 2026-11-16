# SPDX-License-Identifier: GPL-3.0-or-later
"""Tests for paralog2viewer.plot with synthetic BAM/VCF data."""

from __future__ import annotations

import json
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


def test_read_hits_follow_the_rules(data: dict) -> None:
    alleles = plot.read_alleles(data["vcf"], data["keys"])
    hap2 = plot.load_reads(data["dirs"]["R_hap2"] / "reads.bam")
    hap1 = plot.load_reads(data["dirs"]["R_hap1"] / "reads.bam")

    def count(hap: str, reads: plot.HapReads) -> dict[int, int]:
        carriers = [v for v in alleles if hap in v.haps]
        positions = [v.pos for v in carriers]
        out: dict[int, int] = {}
        for group in reads.rows:
            for v in plot.read_hits(group, carriers, positions):
                out[v.pos] = out.get(v.pos, 0) + 1
        return out

    # Reads span 100-1090: all 40 carry 500; only reads reaching base 900
    # (read index >= 20) carry 900.
    assert count("R_hap2", hap2) == {500: 40, 600: 40, 900: 20}
    assert count("R_hap1", hap1) == {600: 40}


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


def test_single_haplotype_marks_only_the_known_site(
    data: dict, monkeypatch: pytest.MonkeyPatch
) -> None:
    markers = set()
    original = Axes.scatter

    def spy(self, x, y, *args, **kwargs):  # noqa: ANN001
        markers.add(kwargs["marker"])
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
    assert markers == {"s"}


def test_hap_names_roles_and_order() -> None:
    table = plot.load_roles()
    assert plot.parse_hap_name("pms2_pms2clhap2") == ("pms2", "pms2clhap", 2)
    assert plot.parse_hap_name("3p7delhap1") == ("", "3p7delhap", 1)
    assert plot.parse_hap_name("weird") == ("", "", 0)
    pseudo = plot.resolve_role("pms2_pms2clhap1", table)
    assert (pseudo.symbol, pseudo.role) == ("PMS2CL", "pseudogene")
    fusion = plot.resolve_role("CYP2D6_fusionhap1", table)
    assert (fusion.symbol, fusion.role) == ("CYP2D6", "fusion")
    copy = plot.resolve_role("AMY1A_hap3", table)
    assert (copy.symbol, copy.role) == ("AMY1A", "unresolved")
    assert plot.resolve_role("new_gene_xhap1", table).role == "unresolved"
    roles = {
        h: plot.resolve_role(h, table)
        for h in ("pms2_pms2clhap1", "pms2_pms2hap2", "pms2_pms2hap1")
    }
    assert plot.order_haps(list(roles), roles) == [
        "pms2_pms2hap1",
        "pms2_pms2hap2",
        "pms2_pms2clhap1",
    ]
    assert all(r.role in plot.ROLE_ORDER for r in table.values())


def test_read_cap_is_deterministic_and_a_subset(data: dict) -> None:
    reads = plot.load_reads(data["dirs"]["R_hap2"] / "reads.bam")
    first = plot.sample_reads(reads.rows, 10, "R_hap2")
    again = plot.sample_reads(reads.rows, 10, "R_hap2")
    other = plot.sample_reads(reads.rows, 10, "R_hap1")
    names = [g[0].query_name for g in first]
    assert names == [g[0].query_name for g in again]
    assert names != [g[0].query_name for g in other]
    assert len(set(names)) == 10
    assert plot.sample_reads(reads.rows, 100, "k") == list(reads.rows)


def test_pack_rows_never_overlaps_within_a_row() -> None:
    spans = [(0, 100), (50, 150), (400, 500), (110, 200), (0, 90)]
    rows = plot.pack_rows(spans, gap=10)
    for row in set(rows):
        mine = sorted(s for s, r in zip(spans, rows, strict=True) if r == row)
        assert all(
            b[0] >= a[1] + 10 for a, b in zip(mine, mine[1:], strict=False)
        )
    assert rows[4] == 0
    assert max(rows) + 1 <= 3


def test_gene_models_pick_default_transcript_and_overlap(
    tmp_path: Path,
) -> None:
    def row(acc: str, strand: str, start: int, exons: list, sym: str) -> str:
        starts = ",".join(str(a) for a, _ in exons) + ","
        ends = ",".join(str(b) for _, b in exons) + ","
        return "\t".join(
            ["0", acc, "chr1", strand, str(start), str(exons[-1][1]),
             str(start), str(exons[-1][1]), str(len(exons)), starts, ends,
             "0", sym, "cmpl", "cmpl", "0"]
        )  # fmt: skip

    table = tmp_path / "refGene.txt"
    table.write_text(
        "\n".join(
            [
                row("NM_1.2", "-", 100, [(100, 200), (300, 400)], "AAA"),
                row(
                    "NM_2",
                    "-",
                    100,
                    [(100, 200), (300, 400), (500, 900)],
                    "AAA",
                ),
                row("NM_3", "+", 1000, [(1000, 1100)], "BBB"),
                row("NM_4", "+", 5000, [(5000, 5100)], "FAR"),
            ]
        )
        + "\n"
    )
    default = tmp_path / "default.tsv"
    default.write_text("Name\tChr\tDefaultTranscript\nAAA\t1\tNM_1.2\n")
    chosen = plot.read_gene_models(
        table, plot.read_default_transcripts(default), "chr1", 0, 1200
    )
    assert [(t.symbol, t.accession) for t in chosen] == [
        ("AAA", "NM_1.2"),
        ("BBB", "NM_3"),
    ]
    assert chosen[0].strand == "-"
    assert chosen[0].exons == ((100, 200), (300, 400))
    # Without a default, the longest transcript is used.
    longest = plot.read_gene_models(table, None, "1", 0, 1200)
    assert longest[0].accession == "NM_2"
    assert (
        plot.read_gene_models(table, None, "chr1", 0, 1200, max_genes=1)[
            0
        ].symbol
        == "AAA"
    )


def test_span_ignores_stub_reads(tmp_path: Path) -> None:
    _write_bam(tmp_path / "s", 100, 20, (), length=600)
    reads = plot.load_reads(tmp_path / "s" / "reads.bam")
    stub = pysam.AlignedSegment()
    stub.query_name = "stub"
    stub.reference_start = 2500
    stub.cigarstring = "30M"
    reads.rows.append([stub])
    start, end = plot.span(reads)
    assert (start, end) == (100, 100 + 19 * 10 + 600)


def test_load_hap_info_reads_truncation_and_copies(tmp_path: Path) -> None:
    path = tmp_path / "p.json"
    path.write_text(
        json.dumps(
            {
                "pms2": {
                    "two_copy_haplotypes": ["12"],
                    "final_haplotypes": {"12": "pms2_pms2hap1"},
                    "haplotype_details": {
                        "pms2_pms2hap1": {"boundary": [1, 2]},
                        "pms2_pms2clhap1": {"is_truncated": ["5p"]},
                    },
                },
                "other": {"haplotype_details": None},
            }
        )
    )
    info = plot.load_hap_info(path)
    assert info["pms2_pms2hap1"]["two_copy"] is True
    assert info["pms2_pms2clhap1"]["is_truncated"] == ["5p"]
    assert plot._truncation_text(["5p", "3p"]) == "5' and 3' end truncated"
    assert set(info) == {"pms2_pms2hap1", "pms2_pms2clhap1"}


def test_figure_height_is_bounded_with_many_haplotypes(tmp_path: Path) -> None:
    fa = tmp_path / "ref.fa"
    fa.write_text(f">chr1\n{REF_SEQ}\n")
    pysam.faidx(str(fa))
    known = tmp_path / "known.tsv"
    known.write_text(KNOWN_HEADER)
    haps = [f"AMY_hap{i}" for i in range(1, 13)]
    dirs = {}
    for i, h in enumerate(haps):
        _write_bam(tmp_path / h, 100 + i * 5, 80, (500 + i,))
        dirs[h] = tmp_path / h
    vcf = tmp_path / "m.vcf.gz"
    _write_vcf(
        vcf,
        haps,
        [
            (500 + i, ["1" if j == i else "0" for j in range(12)])
            for i in range(12)
        ],
    )
    refgene = tmp_path / "refGene.txt"
    refgene.write_text(
        "0\tNM_9\tchr1\t+\t200\t2500\t200\t2500\t2\t200,2000,\t500,2500,\t0\tGENEX"
        "\tcmpl\tcmpl\t0,0,\n"
    )
    info = tmp_path / "p.json"
    info.write_text(
        json.dumps(
            {
                "AMY": {
                    "haplotype_details": {
                        h: {"is_truncated": ["3p"]} for h in haps
                    }
                }
            }
        )
    )
    out = plot.plot_regions(
        vcf,
        dirs,
        set(),
        known,
        fa,
        tmp_path / "figs",
        gene_models=refgene,
        hap_info=info,
        max_reads=25,
    )
    with Image.open(out[0]) as img:
        assert img.size[1] / 300 * 25.4 <= 247.0 + 0.1
        assert abs(img.size[0] - 180 / 25.4 * 300) <= 1


def _loci(haps: list[list[str]]) -> list[dict]:
    return [{"haps": h} for h in haps]


def test_locus_file_names_use_symbols_and_fall_back() -> None:
    table = plot.load_roles()
    names = ["pms2_pms2hap1", "pms2_pms2clhap1", "AMY1A_hap1", "AMY1A_hap2"]
    roles = {h: plot.resolve_role(h, table) for h in names}
    assert plot.locus_file_names(
        "pms2", _loci([["pms2_pms2hap1"], ["pms2_pms2clhap1"]]), roles
    ) == ["pms2_PMS2.jpeg", "pms2_PMS2CL.jpeg"]
    assert plot.locus_file_names("pms2", _loci([["pms2_pms2hap1"]]), roles) == [
        "pms2.jpeg"
    ]
    # Unresolved haplotypes carry no symbol: locus numbers.
    assert plot.locus_file_names(
        "AMY1A", _loci([["AMY1A_hap1"], ["AMY1A_hap2"]]), roles
    ) == ["AMY1A_locus1.jpeg", "AMY1A_locus2.jpeg"]
    # Two loci with the same symbol collide: locus numbers.
    twins = _loci([["pms2_pms2hap1"], ["pms2_pms2hap1"]])
    assert plot.locus_file_names("pms2", twins, roles) == [
        "pms2_locus1.jpeg",
        "pms2_locus2.jpeg",
    ]


def test_split_region_figures_named_by_symbol(tmp_path: Path) -> None:
    fa = tmp_path / "ref.fa"
    fa.write_text(f">chr1\n{REF_SEQ}\n")
    pysam.faidx(str(fa))
    known = tmp_path / "known.tsv"
    known.write_text(KNOWN_HEADER)
    dirs = {"pms2_pms2hap1": tmp_path / "a", "pms2_pms2clhap1": tmp_path / "b"}
    dirs_start = {"pms2_pms2hap1": 100, "pms2_pms2clhap1": 1800}
    for hap, folder in dirs.items():
        _write_bam(folder, dirs_start[hap], 10, ())
    vcf = tmp_path / "m.vcf.gz"
    _write_vcf(vcf, list(dirs), [(500, ["1", "0"])])
    out = plot.plot_regions(
        vcf, dirs, set(), known, fa, tmp_path / "figs", workers=1
    )
    assert [p.name for p in out] == ["pms2_PMS2.jpeg", "pms2_PMS2CL.jpeg"]


def test_same_role_haplotypes_are_colored_by_carrier(
    data: dict, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Two gene haplotypes: ticks take per-haplotype colors, not role blue."""
    from matplotlib.collections import LineCollection

    colors = []
    original = LineCollection.__init__

    def spy(self, segments, *args, **kwargs):  # noqa: ANN001
        colors.append(kwargs.get("colors"))
        original(self, segments, *args, **kwargs)

    monkeypatch.setattr(LineCollection, "__init__", spy)
    dirs = {
        "pms2_pms2hap1": data["dirs"]["R_hap1"],
        "pms2_pms2hap2": data["dirs"]["R_hap2"],
    }
    vcf = data["vcf"].parent / "same.vcf.gz"
    _write_vcf(vcf, list(dirs), [(500, ["0", "1"]), (700, ["1", "0"])])
    plot.plot_regions(
        vcf, dirs, set(), data["known"], data["fa"], data["out"], workers=1
    )
    assert plot.HAP_COLORS[0] in colors
    assert plot.HAP_COLORS[1] in colors
    assert plot.ROLE_COLOR["gene"] not in colors
