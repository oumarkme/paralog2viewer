# SPDX-License-Identifier: GPL-3.0-or-later
"""Tests for paralog2viewer.rccx with synthetic Paraphase, VCF and FASTA."""

from __future__ import annotations

import csv
import json
import random
import shutil
from pathlib import Path

import pysam
import pytest
from PIL import Image

from paralog2viewer import rccx

BLASTN = shutil.which("blastn")
N_SITES = 20
BOUNDARY = (32013300, 32046200)
FAMILIES = ("RP", "C4", "CYP21", "TNX")


# ------------------------------------------------------------ gene table


def test_gene_table_matches_spec():
    table = rccx.load_gene_table()
    boxes = {g.name: (g.start, g.end) for g in table.genes}
    assert boxes["WHR1"] == (31971091, 31982821)
    assert boxes["C4A"] == (31981991, 32004877)
    assert boxes["CYP21A1P"] == (32005689, 32008399)
    assert boxes["TNXA"] == (32008614, 32012472)
    assert boxes["WHR1B"] == (32013270, 32013787)
    assert boxes["C4B"] == (32014730, 32035973)
    assert boxes["CYP21A2"] == (32038298, 32042321)
    assert boxes["TNXB"] == (32041153, 32115334)
    assert table.families() == list(FAMILIES)
    prox, dist = table.pair("TNX")
    assert (prox.name, dist.name) == ("TNXA", "TNXB")
    assert table.gene("WHR1").display == "RP1"
    assert table.marker_window == ("chr6", 32037527, 32040421)
    assert table.promoter_sites == (32038297, 32038310, 32038313, 32038320)


# ------------------------------------------------------------ precompute


def test_blast_mismatches_skip_gaps_and_ambiguity():
    text = "1\t8\t1\t8\tACGTAC-T\tACCTANGT\n1\t1\t1\t1\tA\tT\n"
    diffs = rccx._blast_mismatches(text, 1000, 5000)
    # First HSP: column 3 mismatches (1002, subject 5002), column 6 has N
    # in the subject (1005), column 7 is a query gap. The second HSP adds
    # 1000 but must not overwrite 1002.
    assert diffs[1002] == ("G", "C", 5002)
    assert diffs[1000] == ("A", "T", 5000)
    assert 1005 not in diffs


def test_blast_mismatches_skip_low_identity_alignments():
    # 20 columns, 1 mismatch (95 %) is kept; 10 columns, 3 mismatches
    # (70 %) is a repeat hit and is dropped at min_identity 0.95.
    good = "1\t20\t1\t20\t" + "A" * 19 + "C\t" + "A" * 19 + "G\n"
    bad = "1\t10\t1\t10\tAAAAAAAAAA\tAAAAAAACCC\n"
    assert set(rccx._blast_mismatches(good + bad, 100, 0, 0.95)) == {119}
    assert len(rccx._blast_mismatches(bad, 100, 0, 0.0)) == 3


def test_blast_mismatch_subject_position_after_gaps_and_strand():
    # A query gap consumes a subject base; a subject gap does not.
    gapped = "1\t5\t1\t5\tAC-GT\tACTCT\n"
    assert rccx._blast_mismatches(gapped, 100, 500) == {102: ("G", "C", 503)}
    sgap = "1\t4\t1\t3\tACGT\tAC-A\n"
    assert rccx._blast_mismatches(sgap, 100, 500) == {103: ("T", "A", 502)}
    # Reverse-strand hit (send < sstart) or unknown start: position 0.
    reverse = "1\t3\t9\t7\tACG\tACT\n"
    assert rccx._blast_mismatches(reverse, 100, 500)[102][2] == 0
    assert rccx._blast_mismatches(gapped, 100)[102][2] == 0


@pytest.mark.skipif(BLASTN is None, reason="blastn not available")
def test_build_rccx_sites_finds_planted_substitutions(tmp_path: Path):
    rng = random.Random(1)
    proximal = "".join(rng.choice("ACGT") for _ in range(2000))
    distal = list(proximal)
    planted = [150, 600, 1000, 1700]
    for off in planted:
        distal[off] = "A" if proximal[off] != "A" else "C"
    pad = lambda n: "".join(rng.choice("ACGT") for _ in range(n))  # noqa: E731
    contig = pad(1000) + proximal + pad(3000) + "".join(distal) + pad(1000)
    fasta = tmp_path / "ref.fa"
    fasta.write_text(f">chr6\n{contig}\n")
    pysam.faidx(str(fasta))
    genes = tmp_path / "genes.tsv"
    genes.write_text(
        "kind\tname\tfamily\trole\tchrom\tstart\tend\tdisplay\n"
        "gene\tPROX\tX\tproximal\tchr6\t1001\t3000\tPROX\n"
        "gene\tDIST\tX\tdistal\tchr6\t6001\t8000\tDIST\n"
    )
    out = tmp_path / "sites.tsv"
    n = rccx.build_rccx_sites(fasta, genes, out, blastn=BLASTN)
    sites = rccx.load_sites(out)["X"]
    assert n == len(sites) == 4
    assert [s.pos for s in sites] == [6001 + o for o in planted]
    for s, off in zip(sites, planted, strict=True):
        assert (s.distal_base, s.proximal_base) == (distal[off], proximal[off])
        assert (s.distal_gene, s.proximal_gene) == ("DIST", "PROX")
        assert s.proximal_pos == 1001 + off
    header = out.read_text().splitlines()[0].split("\t")
    assert tuple(header) == rccx.SITE_COLUMNS


def test_build_rccx_sites_rejects_gene_outside_contig(tmp_path: Path):
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr6\nACGTACGTAC\n")
    pysam.faidx(str(fasta))
    genes = tmp_path / "genes.tsv"
    genes.write_text(
        "kind\tname\tfamily\trole\tchrom\tstart\tend\tdisplay\n"
        "gene\tP\tX\tproximal\tchr6\t1\t5\tP\n"
        "gene\tD\tX\tdistal\tchr6\t6\t50\tD\n"
    )
    with pytest.raises(ValueError, match="outside the reference"):
        rccx.build_rccx_sites(fasta, genes, tmp_path / "s.tsv")


# ---------------------------------------------------------- vote logic


def _vote(sites, votes, span):
    return rccx.vote_haplotype(sites, dict(enumerate(votes)), span)


def _flat_sites(n: int):
    return [rccx.Site("X", "chr6", i, "D", "P", "A", "G") for i in range(n)]


def test_vote_thresholds():
    sites = _flat_sites(20)
    span = (0, 19)
    assert (
        _vote(sites, ["proximal"] * 18 + ["distal"] * 2, span)["call"]
        == "proximal_type"
    )
    assert (
        _vote(sites, ["distal"] * 18 + ["proximal"] * 2, span)["call"]
        == "distal_type"
    )
    # Opposite majorities in both halves, 9/10 and 9/10 -> chimeric.
    chim = ["distal"] * 9 + ["proximal"] + ["proximal"] * 9 + ["distal"]
    res = _vote(sites, chim, span)
    assert res["call"] == "chimeric"
    assert (res["left_majority"], res["right_majority"]) == (
        "distal",
        "proximal",
    )
    # Alternating votes: no clean breakpoint.
    assert _vote(sites, ["proximal", "distal"] * 10, span)["call"] == "mixed"
    # Fewer than 5 informative sites in a half blocks chimeric.
    few = ["other"] * 5 + ["distal"] * 5 + ["proximal"] * 10
    assert _vote(sites, few, span)["call"] == "mixed"
    # Fewer than 3 sites, and no informative vote.
    assert _vote(sites, ["proximal"] * 20, (0, 1))["call"] == "not_determined"
    res = _vote(sites, ["missing"] * 20, span)
    assert res["call"] == "not_determined" and res["n_informative"] == 0


def _like(votes, span=(0, 19), n=20):
    return rccx.vote_haplotype(
        _flat_sites(n), dict(enumerate(votes)), span, like_calls=True
    )


def _mix(n_prox: int, n_dist: int):
    """Votes with the distal ones spread evenly, so no half is a chimera."""
    n = n_prox + n_dist
    return [
        "distal" if (i + 1) * n_dist // n > i * n_dist // n else "proximal"
        for i in range(n)
    ]


def test_like_bands():
    # 19/20 = 0.95 and 1/20 = 0.05 stay confident calls.
    res = _like(_mix(19, 1))
    assert res["call"] == "proximal_type" and res["like_pct"] == ""
    res = _like(_mix(1, 19))
    assert res["call"] == "distal_type" and res["n_matching"] == 19
    # 14/18 = 0.78 proximal-like, 78 %, 14 of 18 covered sites.
    res = _like(_mix(14, 4) + ["other"] * 2)
    assert (res["call"], res["like_pct"]) == ("proximal_like", 78)
    assert (res["n_matching"], res["n_informative"]) == (14, 18)
    # 3/18 = 0.17 proximal -> distal-like with the distal share, 83 %.
    res = _like(_mix(3, 15) + ["missing"] * 2)
    assert (res["call"], res["like_pct"]) == ("distal_like", 83)
    assert (res["n_matching"], res["n_informative"]) == (15, 18)
    # 0.18 proximal -> distal-like 82 %; exactly 0.5 -> mixed 50 %.
    res = _like(_mix(9, 41), n=50, span=(0, 49))
    assert (res["call"], res["like_pct"]) == ("distal_like", 82)
    res = _like(_mix(10, 10))
    assert (res["call"], res["like_pct"]) == ("mixed", 50)
    # Bounds: 0.9 and 0.1 are confident, 0.55 and 0.45 lean.
    assert _like(_mix(18, 2))["call"] == "proximal_type"
    assert _like(_mix(2, 18))["call"] == "distal_type"
    assert _like(_mix(11, 9))["call"] == "proximal_like"
    assert _like(_mix(9, 11))["like_pct"] == 55


def test_like_keeps_not_determined_and_chimeric_first():
    assert _like(["proximal"] * 20, (0, 1))["call"] == "not_determined"
    assert _like(["missing"] * 20)["call"] == "not_determined"
    chim = ["distal"] * 9 + ["proximal"] + ["proximal"] * 9 + ["distal"]
    res = _like(chim)
    assert res["call"] == "chimeric" and res["like_pct"] == ""
    # Without like_calls (CYP21) the same votes stay mixed.
    assert (
        _vote(_flat_sites(20), _mix(14, 4) + ["other"] * 2, (0, 19))["call"]
        == "mixed"
    )


def test_vote_margin_blocks_weak_chimera():
    sites = _flat_sites(20)
    # Left half (positions < 9) is 5 distal / 4 proximal: margin 0.56 < 0.65.
    votes = ["distal"] * 5 + ["proximal"] * 4 + ["proximal"] * 11
    assert _vote(sites, votes, (0, 19))["call"] == "mixed"


def test_site_vote_rules():
    site = rccx.Site("X", "chr6", 5, "D", "P", "A", "G")
    vote = rccx._site_vote
    assert vote(site, [("A", "G", "1")], False, None) == "proximal"
    assert vote(site, [("A", "T", "1")], False, None) == "other"
    assert vote(site, [("A", "G", "0")], False, None) == "distal"
    assert vote(site, [("A", "G", ".")], True, None) == "missing"
    # A second allele at the site that is ALT does not hide the match.
    assert (
        vote(site, [("A", "T", "1"), ("A", "G", "1")], False, None)
        == "proximal"
    )
    # No record: reference only when provably callable.
    assert vote(site, None, True, None) == "distal"
    assert vote(site, None, False, None) == "missing"
    assert vote(site, None, False, True) == "distal"
    assert vote(site, None, True, False) == "missing"


def test_depths_count_hap_bam_reads(tmp_path: Path):
    bam_path = tmp_path / "reads.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr6", "LN": 1000}],
    }
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out:
        for i in range(3):
            seg = pysam.AlignedSegment(out.header)
            seg.query_name = f"r{i}"
            seg.query_sequence = "A" * 100
            seg.reference_id = 0
            seg.reference_start = 99  # 1-based 100..199
            seg.cigarstring = "100M"
            seg.query_qualities = pysam.qualitystring_to_array("I" * 100)
            out.write(seg)
    pysam.index(str(bam_path))
    got = rccx._depths(bam_path, "chr6", [99, 100, 199, 200])
    assert got == {99: 0, 100: 3, 199: 3, 200: 0}


def test_proximal_frame_vote_rules():
    site = rccx.Site("X", "chr6", 5, "D", "P", "A", "G", 50)
    vote = rccx._proximal_frame_vote
    # In the proximal copy the proximal base is the reference.
    assert vote(site, [("G", "A", "1")], False, None) == "distal"
    assert vote(site, [("G", "T", "1")], False, None) == "other"
    assert vote(site, [("G", "A", "0")], False, None) == "proximal"
    assert vote(site, [("G", "A", ".")], True, None) == "missing"
    assert vote(site, None, False, True) == "proximal"
    assert vote(site, None, True, False) == "missing"
    assert vote(site, None, True, None) == "proximal"
    assert vote(site, None, False, None) == "missing"


def test_frame_vote_picks_the_frame_of_the_reads():
    site = rccx.Site("X", "chr6", 5, "D", "P", "A", "G", 50)
    ref_d = [("A", "G", "0")]
    alt_p = [("G", "A", "1")]
    vote = rccx._frame_vote
    # BAM: the frame with the larger depth (>= 3) decides.
    assert vote(site, (ref_d, True, 9), (None, True, 0), 3, True) == "distal"
    assert vote(site, (None, True, 0), (alt_p, True, 9), 3, True) == "distal"
    assert vote(site, (None, True, 0), (None, True, 9), 3, True) == "proximal"
    assert vote(site, (None, True, 4), (None, True, 9), 3, True) == "proximal"
    assert vote(site, (None, True, 9), (None, True, 4), 3, True) == "distal"
    assert vote(site, (None, True, 1), (None, True, 2), 3, True) == "missing"
    # Without BAM: the frame with a genotype decides ...
    nocall = [("A", "G", ".")]
    assert vote(site, (nocall, True, None), (alt_p, True, None), 3, False) == (
        "distal"
    )
    assert vote(site, (ref_d, True, None), (None, True, None), 3, False) == (
        "distal"
    )
    # ... and a site without a genotype is reference only in one frame.
    only_d = vote(site, (None, True, None), (None, False, None), 3, False)
    only_p = vote(site, (None, False, None), (None, True, None), 3, False)
    both = vote(site, (None, True, None), (None, True, None), 3, False)
    assert (only_d, only_p, both) == ("distal", "proximal", "missing")
    # No proximal position: the distal frame only, as before.
    legacy = rccx.Site("X", "chr6", 5, "D", "P", "A", "G")
    assert vote(legacy, (None, True, 9), (None, True, 9), 3, True) == "distal"


# ----------------------------------------------------- end-to-end calls


def _site_positions(table: rccx.GeneTable, family: str) -> list[int]:
    _prox, dist = table.pair(family)
    lo, hi = max(dist.start, BOUNDARY[0]), min(dist.end, BOUNDARY[1])
    step = (hi - lo + 1) / N_SITES
    return [int(lo + (i + 0.5) * step) for i in range(N_SITES)]


def _write_inputs(tmp_path: Path):
    table = rccx.load_gene_table()
    sites_tsv = tmp_path / "sites.tsv"
    positions: dict[str, list[int]] = {}
    with open(sites_tsv, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        # Legacy 7-column table without proximal_pos: distal frame only.
        writer.writerow(rccx.SITE_COLUMNS[:-1])
        for family in FAMILIES:
            prox, dist = table.pair(family)
            positions[family] = _site_positions(table, family)
            for pos in positions[family]:
                writer.writerow(
                    [family, "chr6", pos, dist.name, prox.name, "A", "G"]
                )
    # Patterns by haplotype: share of proximal votes along the site index.
    patterns = {
        "g_hap1": lambda i: 1,
        "g_hap2": lambda i: 0,
        "g_hap3": lambda i: int(i >= N_SITES // 2),
        "g_hap4": lambda i: i % 2,
        "g_hap5": lambda i: 1,
        "g_hap6": lambda i: 0,
    }
    haps = list(patterns)
    lines = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=chr6,length=170000000>",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(haps),
    ]
    recs = []
    for family in FAMILIES:
        for i, pos in enumerate(positions[family]):
            gts = [str(patterns[h](i)) for h in haps]
            recs.append((pos, "A", "G", gts))
    # A private ALT at the first C4 site: g_hap2 carries it, so its vote
    # there is "other" and must not count.
    first = positions["C4"][0]
    recs.append((first, "A", "T", ["0", "1", "0", "0", "0", "0"]))
    for pos, ref, alt, gts in sorted(recs):
        lines.append(
            f"chr6\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t" + "\t".join(gts)
        )
    vcf = tmp_path / "merged.vcf"
    vcf.write_text("\n".join(lines) + "\n")
    merged = Path(pysam.tabix_index(str(vcf), preset="vcf", force=True))

    c4 = table.pair("C4")[1]
    details = {h: {"boundary": list(BOUNDARY)} for h in haps}
    # g_hap5 only spans the first two C4 sites (and none of RP/TNX/CYP21).
    details["g_hap5"]["boundary"] = [c4.start, positions["C4"][1]]
    region = {
        "haplotype_details": details,
        "alleles_final": [
            ["g_hap2", "g_hap1"],
            ["g_hap3"],
            ["g_hap5", "g_hap4", "g_hap6"],
        ],
        "annotated_alleles": ["WT", "del", "tri"],
        "hap_variants": {
            "g_hap1": ["P31L"],
            "g_hap2": [],
            "g_hap3": [],
            "g_hap4": ["Q319X"],
            "g_hap5": [],
            "g_hap6": [],
        },
    }
    pjson = tmp_path / "x.paraphase.json"
    pjson.write_text(json.dumps({"rccx": region}))
    return pjson, merged, sites_tsv


@pytest.fixture
def calls(tmp_path: Path):
    pjson, merged, sites = _write_inputs(tmp_path)
    return rccx.classify_rccx(pjson, merged, sites)


def _by_hap(rows):
    return {r["hap"]: r for r in rows}


def test_module_structure(calls):
    by = _by_hap(calls)
    assert [r["hap"] for r in calls] == [
        "g_hap1",
        "g_hap2",
        "g_hap3",
        "g_hap4",
        "g_hap5",
        "g_hap6",
    ]
    assert (by["g_hap1"]["structure"], by["g_hap1"]["deleted_modules"]) == (
        "bimodular",
        0,
    )
    assert by["g_hap1"]["module_position"] == "proximal"
    assert by["g_hap2"]["module_position"] == "distal"
    assert (
        by["g_hap3"]["structure"],
        by["g_hap3"]["deleted_modules"],
        by["g_hap3"]["module_position"],
    ) == ("monomodular", 1, "sole")
    assert by["g_hap4"]["structure"] == "trimodular"
    assert by["g_hap4"]["allele"] == 3 and by["g_hap4"]["allele_tag"] == "tri"
    assert by["g_hap1"]["allele_tag"] == "WT"


def test_cyp21_paraphase_markers(calls):
    by = _by_hap(calls)
    assert by["g_hap1"]["cyp21_paraphase"] == "CYP21A1P"
    assert by["g_hap1"]["cyp21_marker_hits"] == 1
    assert by["g_hap2"]["cyp21_paraphase"] == "CYP21A2"
    assert by["g_hap2"]["cyp21_marker_hits"] == 0


def test_blast_vote_calls(calls):
    by = _by_hap(calls)
    assert by["g_hap1"]["c4_call"] == "proximal_type"
    assert by["g_hap1"]["c4_gene"] == "C4A"
    assert by["g_hap2"]["c4_call"] == "distal_type"
    assert by["g_hap2"]["c4_gene"] == "C4B"
    assert by["g_hap3"]["c4_call"] == "chimeric"
    assert by["g_hap3"]["c4_gene"] == ""
    assert by["g_hap4"]["c4_call"] == "mixed"
    assert by["g_hap5"]["c4_call"] == "not_determined"
    assert by["g_hap5"]["c4_n_sites"] == 2
    assert by["g_hap1"]["rp_gene"] == "WHR1"
    assert by["g_hap2"]["tnx_gene"] == "TNXB"
    assert by["g_hap1"]["cyp21_blast_gene"] == "CYP21A1P"


def test_like_columns_are_appended_for_rp_c4_and_tnx_only(calls):
    keys = list(calls[0])
    tail = keys[-12:]
    assert tail == [
        f"{pre}_{name}"
        for pre in ("rp", "c4", "tnx")
        for name in (
            "like_gene",
            "like_pct",
            "matching_sites",
            "covered_sites",
        )
    ]
    assert not any(k.startswith("cyp21_blast_like") for k in keys)
    by = _by_hap(calls)
    # Confident call: counts filled, no leaning gene or percent.
    row = by["g_hap1"]
    assert (row["c4_like_gene"], row["c4_like_pct"]) == ("", "")
    assert row["c4_matching_sites"] == row["c4_covered_sites"]
    # Alternating votes are mixed; confident calls keep their old names.
    assert by["g_hap4"]["c4_call"] == "mixed"
    assert by["g_hap4"]["c4_like_pct"] in (50, "")
    assert by["g_hap1"]["rp_call"] == "proximal_type"
    assert by["g_hap2"]["cyp21_blast_call"] == "distal_type"


def _set_call(row, pre, call, like_gene, pct, matching, covered):
    row[f"{pre}_call"] = call
    row[f"{pre}_gene"] = ""
    row[f"{pre}_like_gene"] = like_gene
    row[f"{pre}_like_pct"] = pct
    row[f"{pre}_matching_sites"] = matching
    row[f"{pre}_covered_sites"] = covered


def synthetic_rows(calls):
    """Copy of the fixture rows with one of each C4/TNX call kind."""
    rows = [dict(r) for r in calls]
    plan = [
        ("c4", "proximal_like", "C4A", 78, 14, 18),
        ("c4", "distal_like", "C4B", 83, 15, 18),
        ("c4", "mixed", "", 50, 9, 18),
        ("tnx", "proximal_like", "TNXA", 61, 11, 18),
        ("tnx", "distal_like", "TNXB", 72, 13, 18),
        ("tnx", "mixed", "", 50, 9, 18),
    ]
    for row, (pre, *rest) in zip(rows[1:], plan, strict=False):
        _set_call(row, pre, *rest)
    return rows


def test_like_labels_table_cells_and_boxes(calls, tmp_path: Path):
    rows = synthetic_rows(calls)
    table = rccx.load_gene_table()
    assert rccx._label_for(rows[1], "C4") == ("like", "C4A", "78")
    assert rccx._label_for(rows[3], "C4") == ("mixed", "", "50")
    assert rccx._label_for(rows[0], "C4")[0] == "gene"
    assert rccx._table_cells(rows[1], table, True)[6] == "C4A-like (78%, 14/18)"
    assert rccx._table_cells(rows[3], table, True)[6] == "mixed (50%, 9/18)"
    assert rccx._table_cells(rows[0], table, True)[6] == "C4A"
    with rccx.plt.rc_context(rccx._set_font()):
        fig = rccx._build_figure(rows, table)
        fig.canvas.draw()
        schema = fig.axes[0]
        texts = {t.get_text() for t in schema.texts}
        assert {"C4A-like\n78%", "C4B-like\n83%", "mixed\n50%"} <= texts
        # TNXA-like is wider than a 3-module box: short names (A-like).
        assert {"A-like\n61%", "B-like\n72%"} <= texts
        thick = [
            p
            for p in schema.patches
            if p.get_linewidth() == rccx.LIKE_LINEWIDTH
        ]
        assert len(thick) == 4
        assert all(p.get_facecolor() == (1.0, 1.0, 1.0, 1.0) for p in thick)
        edges = {
            rccx.matplotlib.colors.to_hex(p.get_edgecolor()) for p in thick
        }
        assert edges == {"#0072b2", "#56b4e9", "#cc79a7", "#000000"}
        # Every box label stays inside its box.
        for text in schema.texts:
            if "-like" in text.get_text() or text.get_text() == "mixed\n50%":
                bb = text.get_window_extent()
                inside = [
                    p.get_window_extent()
                    for p in schema.patches
                    if p.get_window_extent().contains(
                        *(
                            (bb.x0 + bb.x1) / 2,
                            (bb.y0 + bb.y1) / 2,
                        )
                    )
                ]
                assert any(
                    b.x0 <= bb.x0
                    and bb.x1 <= b.x1
                    and b.y0 <= bb.y0
                    and bb.y1 <= b.y1
                    for b in inside
                ), (text.get_text(), bb, inside)
        legend_text = " ".join(
            t.get_text() for t in fig.axes[-1].get_legend().get_texts()
        )
        assert "-like: majority type below 90 %" in legend_text
        rccx.plt.close(fig)
    out = rccx.plot_rccx_module(rows, tmp_path / "like.jpeg")
    assert out.exists()


def test_private_allele_is_not_counted(calls):
    row = _by_hap(calls)["g_hap2"]
    assert row["c4_n_sites"] == N_SITES
    assert row["c4_n_informative"] == N_SITES - 1
    assert row["c4_call"] == "distal_type"


def test_rows_are_json_free_scalars_and_table_roundtrip(calls, tmp_path):
    out = tmp_path / "t" / "calls.tsv"
    rccx.write_rccx_table(calls, out)
    back = list(csv.DictReader(open(out, encoding="utf-8"), delimiter="\t"))
    assert len(back) == len(calls)
    assert back[0]["hap"] == "g_hap1" and back[0]["c4_gene"] == "C4A"
    with pytest.raises(ValueError):
        rccx.write_rccx_table([], out)


def test_multiallelic_record_is_rejected(tmp_path: Path):
    pjson, _merged, sites = _write_inputs(tmp_path)
    vcf = tmp_path / "multi.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr6,length=170000000>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tg_hap1\n"
        "chr6\t32020000\t.\tA\tG,T\t.\tPASS\t.\tGT\t1\n"
    )
    with pytest.raises(ValueError, match="one record per allele"):
        rccx.classify_rccx(pjson, vcf, sites)


# ------------------------------------------------------------ the figure


def test_figure_geometry_and_legibility(calls, tmp_path: Path):
    out = rccx.plot_rccx_module(calls, tmp_path / "figures" / "m.jpeg")
    with Image.open(out) as img:
        assert img.format == "JPEG"
        assert round(img.info["dpi"][0]) == 150
        assert img.width == round(rccx.WIDTH_MM / rccx.MM * 150)
        assert img.height / 150 * rccx.MM <= rccx.MAX_HEIGHT_MM

    table = rccx.load_gene_table()
    with rccx.plt.rc_context(rccx._set_font()):
        fig = rccx._build_figure(calls, table)
        fig.canvas.draw()
        texts = [t for t in fig.findobj(rccx.plt.Text) if t.get_text()]
        assert texts and all(t.get_fontsize() >= 7 for t in texts)
        assert all(
            fig.bbox.x0 <= t.get_window_extent().x0
            and t.get_window_extent().x1 <= fig.bbox.x1
            for t in texts
        )
        schema_ax = fig.axes[0]
        boxes = [t.get_window_extent() for t in schema_ax.texts]
        assert all(fig.bbox.x0 <= b.x0 and b.x1 <= fig.bbox.x1 for b in boxes)
        for i, a in enumerate(boxes):
            for b in boxes[i + 1 :]:
                assert not a.overlaps(b), "schematic texts overlap"
        names = {t.get_text() for t in schema_ax.texts}
        assert "deleted module" in names and "chimeric" in names
        assert "g_hap3" in names
        rccx.plt.close(fig)


def test_legend_markers_are_unique_and_not_gene_colors(calls):
    table = rccx.load_gene_table()
    with rccx.plt.rc_context(rccx._set_font()):
        fig = rccx._build_figure(calls, table)
        legend = fig.axes[-1].get_legend()
        keys = [
            (
                p.get_facecolor(),
                p.get_edgecolor(),
                p.get_hatch(),
                p.get_linestyle(),
                p.get_linewidth(),
            )
            for p in legend.get_patches()
        ]
        rccx.plt.close(fig)
    assert len(keys) == len(set(keys)), "two legend entries look alike"
    gene_rgba = {
        rccx.matplotlib.colors.to_rgba(c) for c in rccx.GENE_COLORS.values()
    }
    assert not {k[0] for k in keys} & gene_rgba, "legend color is a gene color"


def test_conflict_flag_and_empty_rows(tmp_path: Path):
    assert rccx._conflict(
        {"cyp21_blast_gene": "CYP21A2", "cyp21_paraphase": "CYP21A1P"}
    )
    assert not rccx._conflict(
        {"cyp21_blast_gene": "", "cyp21_paraphase": "CYP21A1P"}
    )
    with pytest.raises(ValueError, match="no haplotype rows"):
        rccx.plot_rccx_module([], tmp_path / "x.jpeg")


# ------------------------------------------------ reads in the proximal copy


def _write_flat_bam(path: Path, start: int, end: int, n_reads: int = 3):
    """Write a BAM with n_reads reads covering 1-based [start, end]."""
    path.parent.mkdir(parents=True, exist_ok=True)
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr6", "LN": 170000000}],
    }
    length = end - start + 1
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for i in range(n_reads):
            seg = pysam.AlignedSegment(out.header)
            seg.query_name = f"r{i}"
            seg.query_sequence = "A" * length
            seg.reference_id = 0
            seg.reference_start = start - 1
            seg.cigarstring = f"{length}M"
            seg.query_qualities = pysam.qualitystring_to_array("I" * length)
            out.write(seg)
    pysam.index(str(path))


def _write_frame_inputs(tmp_path: Path):
    """Three haplotypes: reads in the distal copy, or in the proximal copy.

    h_dist: reads and genotype 0 at the distal sites (distal type).
    h_prox: reads in the proximal copy, reference there (proximal type).
    h_conv: reads in the proximal copy carrying the distal base (distal
    type), as after realignment of a gene-converted haplotype.
    The Paraphase JSON has no phased alleles (phasing failed).
    """
    table = rccx.load_gene_table()
    positions = _site_positions(table, "C4")
    offset = 30000
    sites_tsv = tmp_path / "sites.tsv"
    with open(sites_tsv, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(rccx.SITE_COLUMNS)
        for pos in positions:
            writer.writerow(
                ["C4", "chr6", pos, "C4B", "C4A", "A", "G", pos - offset]
            )
    haps = ["h_dist", "h_prox", "h_conv"]
    lines = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=chr6,length=170000000>",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(haps),
    ]
    recs = [(pos - offset, "G", "A", [".", "0", "1"]) for pos in positions]
    recs += [(pos, "A", "G", ["0", ".", "."]) for pos in positions]
    for pos, ref, alt, gts in sorted(recs):
        lines.append(
            f"chr6\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t" + "\t".join(gts)
        )
    vcf = tmp_path / "merged.vcf"
    vcf.write_text("\n".join(lines) + "\n")
    merged = Path(pysam.tabix_index(str(vcf), preset="vcf", force=True))
    bam_dir = tmp_path / "hap"
    lo, hi = positions[0] - 50, positions[-1] + 50
    _write_flat_bam(bam_dir / "h_dist" / "reads.bam", lo, hi)
    for hap in ("h_prox", "h_conv"):
        _write_flat_bam(bam_dir / hap / "reads.bam", lo - offset, hi - offset)
    region = {
        "haplotype_details": {h: {"boundary": list(BOUNDARY)} for h in haps},
        "alleles_final": [],
        "annotated_alleles": [],
        "hap_variants": {"h_dist": [], "h_prox": ["P31L"], "h_conv": []},
    }
    pjson = tmp_path / "frame.paraphase.json"
    pjson.write_text(json.dumps({"rccx": region}))
    return pjson, merged, sites_tsv, bam_dir


@pytest.mark.parametrize("use_bam", [True, False])
def test_reads_realigned_to_the_proximal_copy_are_voted_there(
    tmp_path: Path, use_bam: bool
):
    pjson, merged, sites, bam_dir = _write_frame_inputs(tmp_path)
    rows = rccx.classify_rccx(
        pjson, merged, sites, hap_bam_dir=bam_dir if use_bam else None
    )
    by = _by_hap(rows)
    assert {h: r["c4_call"] for h, r in by.items()} == {
        "h_dist": "distal_type",
        "h_prox": "proximal_type",
        "h_conv": "distal_type",
    }
    assert all(r["c4_n_informative"] == N_SITES for r in rows)
    assert by["h_prox"]["coverage_source"] == ("bam" if use_bam else "vcf")


def test_unassigned_haplotypes_are_placed_by_cyp21_marker(tmp_path: Path):
    pjson, merged, sites, _bam = _write_frame_inputs(tmp_path)
    by = _by_hap(rccx.classify_rccx(pjson, merged, sites))
    assert {h: r["structure"] for h, r in by.items()} == dict.fromkeys(
        by, "unassigned"
    )
    assert {h: r["module_position"] for h, r in by.items()} == {
        "h_dist": "distal",
        "h_prox": "proximal",
        "h_conv": "distal",
    }


def test_rp_like_label_uses_display_name(calls):
    rows = [dict(r) for r in calls]
    _set_call(rows[1], "rp", "distal_like", "WHR1B", 67, 2, 3)
    table = rccx.load_gene_table()
    assert rccx._label_for(rows[1], "RP") == ("like", "WHR1B", "67")
    assert rccx._table_cells(rows[1], table, True)[5] == "RP2-like (67%, 2/3)"
    with rccx.plt.rc_context(rccx._set_font()):
        fig = rccx._build_figure(rows, table)
        fig.canvas.draw()
        texts = {t.get_text() for t in fig.axes[0].texts}
        rccx.plt.close(fig)
    assert "RP2-like\n67%" in texts


@pytest.mark.parametrize("use_bam", [True, False])
def test_failed_haplotype_is_not_determined_not_distal(
    tmp_path: Path, use_bam: bool
):
    """In the JSON but neither a VCF sample nor a BAM: no vote at all."""
    pjson, merged, sites, bam_dir = _write_frame_inputs(tmp_path)
    data = json.loads(pjson.read_text())
    data["rccx"]["haplotype_details"]["h_gone"] = {"boundary": list(BOUNDARY)}
    data["rccx"]["hap_variants"]["h_gone"] = []
    pjson.write_text(json.dumps(data))
    rows = rccx.classify_rccx(
        pjson, merged, sites, hap_bam_dir=bam_dir if use_bam else None
    )
    gone = _by_hap(rows)["h_gone"]
    assert gone["c4_call"] == "not_determined"
    assert gone["c4_n_informative"] == 0
    assert _by_hap(rows)["h_dist"]["c4_call"] == "distal_type"


def test_cyp21_symbols_follow_the_marker_hits(tmp_path: Path):
    path = tmp_path / "p.json"
    path.write_text('{"rccx": {"hap_variants": {"h1": ["m1"], "h2": []}}}')
    prox, dist = rccx.load_gene_table().pair("CYP21")
    assert rccx.cyp21_symbols(path) == {"h1": prox.name, "h2": dist.name}
    assert prox.name != dist.name
    path.write_text("{}")
    assert rccx.cyp21_symbols(path) == {}


def _is_distal(hap: str) -> bool:
    return hap in {"d1", "d2", "d3"}


def test_unlinked_distal_module_joins_the_allele_without_a_distal_end():
    # 2 + 3 module pattern: bimodular (p4 + d1) and two proximal modules
    # with the lone distal module unlinked -> 2 + 3 modules.
    groups, placed, left = rccx._infer_unlinked(
        [["p4", "d1"], ["p2", "p5"]], ["d3"], _is_distal
    )
    assert groups == [["p4", "d1"], ["p2", "p5", "d3"]]
    assert placed == {"d3"} and left == []


def test_unlinked_modules_stay_unlinked_when_ambiguous():
    # Both alleles already have a distal end: nothing to decide on.
    groups, placed, left = rccx._infer_unlinked(
        [["p1", "d1"], ["p2", "d2"]], ["d3"], _is_distal
    )
    assert placed == set() and left == ["d3"] and len(groups) == 2
    # An unlinked proximal module could sit on either allele.
    _g, placed, left = rccx._infer_unlinked(
        [["p1", "d1"], ["p2"]], ["p3"], _is_distal
    )
    assert placed == set() and left == ["p3"]


def test_one_linked_allele_and_a_lone_distal_module_make_two_alleles():
    groups, placed, left = rccx._infer_unlinked(
        [["p1", "d1"]], ["d2"], _is_distal
    )
    assert groups == [["p1", "d1"], ["d2"]] and placed == {"d2"}
    assert left == []


def test_module_text_shows_order_and_tied_range():
    row = {"module_order": "1-2", "n_modules": 3, "module_position": "proximal"}
    assert rccx._module_text(row) == "1-2/3 proximal"
    row = {"module_order": "3", "n_modules": 3, "module_position": "distal"}
    assert rccx._module_text(row) == "3/3 distal"
    row = {"module_order": "1", "n_modules": 1, "module_position": "sole"}
    assert rccx._module_text(row) == "sole"


def test_rows_carry_module_order_and_assignment(calls):
    by = _by_hap(calls)
    assert (by["g_hap1"]["module_order"], by["g_hap2"]["module_order"]) == (
        "1",
        "2",
    )
    assert by["g_hap1"]["assignment"] == "paraphase"
