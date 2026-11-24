# SPDX-License-Identifier: GPL-3.0-or-later
"""Tests for paralog2viewer.summary with synthetic run files."""

from __future__ import annotations

import base64
import csv
import datetime
import gzip
import json
import re
from pathlib import Path

from paralog2viewer import summary

JPEG = b"\xff\xd8\xff\xe0synthetic-jpeg-bytes\xff\xd9"
WHEN = datetime.datetime(2026, 1, 2, 3, 4, tzinfo=datetime.UTC)
COLUMNS = (
    "Chr Start Ref Alt Gene AAChange Genotype ClinVar InterVar Consequence "
    "Max_AF REVEL_score CADD_phred"
).split()


def _row(**kwargs: str) -> dict[str, str]:
    row = dict.fromkeys(COLUMNS, "")
    row["Max Allele Freq"] = row.pop("Max_AF")
    row.update(kwargs)
    return row


def _inter(label: str) -> str:
    return f" InterVar: {label} PVS1=0 PS=[0] PM=[0] "


# ---------------------------------------------------------- flag rule


def test_flag_rule_each_branch_and_the_negatives() -> None:
    flag = summary.flag_reasons
    assert flag(_row(ClinVar="Pathogenic/Likely pathogenic"))
    assert flag(_row(ClinVar="Likely_pathogenic"))
    assert flag(_row(InterVar=_inter("Likely pathogenic")))
    # Not flagged: benign, conflicting, VUS without a score.
    assert not flag(_row(ClinVar="Benign/Likely benign"))
    assert not flag(
        _row(ClinVar="Conflicting_interpretations_of_pathogenicity")
    )
    assert not flag(_row(InterVar=_inter("Uncertain significance")))
    # High impact needs a rare allele; a missing frequency counts as rare.
    frameshift = {"Consequence": "frameshift_variant;splice_region_variant"}
    assert flag(_row(**frameshift, **{"Max Allele Freq": "0.0002"}))
    assert flag(_row(**frameshift))
    assert not flag(_row(**frameshift, **{"Max Allele Freq": "0.05"}))
    assert not flag(_row(Consequence="missense_variant"))
    # VUS: REVEL or CADD, rare only; boundaries are exclusive.
    vus = {"InterVar": _inter("Uncertain significance")}
    assert flag(_row(**vus, REVEL_score="0.71"))
    assert flag(_row(**vus, CADD_phred="25.1"))
    assert not flag(_row(**vus, REVEL_score="0.7", CADD_phred="25"))
    assert not flag(
        _row(**vus, REVEL_score="0.9", **{"Max Allele Freq": "0.02"})
    )
    assert not flag(_row(**vus, REVEL_score="."))


# ----------------------------------------------------------- the page


def _run_files(tmp_path: Path, *, with_figures: bool = True) -> dict:
    tmp_path.mkdir(parents=True, exist_ok=True)
    work = tmp_path / "work"
    work.mkdir()
    audit = {
        "genes": [
            {"gene": "A", "resolved_haplotypes": 2, "status": "resolved"},
            {
                "gene": "<b>B",
                "resolved_haplotypes": 0,
                "status": "no_resolved_haplotypes",
            },
            {"gene": "C", "resolved_haplotypes": 2, "status": "resolved"},
            {"gene": "D_x", "resolved_haplotypes": 1, "status": "resolved"},
        ],
        "haplotypes": [
            {
                "region": r,
                "hap": h,
                "matched_primary_read_ids": n,
                "status": "ready",
            }
            for r, h, n in (
                ("A", "A_hap1", 10),
                ("A", "A_hap2", 5),
                ("C", "C_hap1", 3),
                ("C", "C_hap2", 4),
                ("D_x", "D_x_hap1", 1),
            )
        ],
    }
    (work / "split_audit.json").write_text(json.dumps(audit))
    paraphase = tmp_path / "s.json"
    paraphase.write_text(
        json.dumps(
            {
                "A": {
                    "final_haplotypes": {"a": "A_hap1", "b": "A_hap2"},
                    "total_cn": 2,
                },
                "C": {
                    "final_haplotypes": {"a": "C_hap1", "b": "C_hap2"},
                    "total_cn": 3,
                },
                "D_x": {"final_haplotypes": {"a": "D_x_hap1"}},
            }
        )
    )
    rows = [
        _row(
            Chr="1",
            Start="100",
            Ref="C",
            Alt="T",
            Gene="G1",
            Genotype="A_hap1",
            ClinVar="Pathogenic",
            AAChange="p.R1X",
        ),
        _row(
            Chr="1", Start="200", Ref="A", Alt="G", Gene="G1", Genotype="A_hap2"
        ),
        _row(
            Chr="1", Start="300", Ref="G", Alt="-", Gene="G2", Genotype="A_hap1"
        ),
        _row(
            Chr="1",
            Start="400",
            Ref="T",
            Alt="C",
            Gene="<i>",
            Genotype="C_hap1",
            InterVar=_inter("Pathogenic"),
        ),
    ]
    table = tmp_path / "S.paralog.tsv.gz"
    with gzip.open(table, "wt", newline="") as handle:
        handle.write("#sex=Male\n")
        writer = csv.DictWriter(handle, rows[0], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    figures = tmp_path / "figures"
    if with_figures:
        figures.mkdir()
        (figures / "A.jpeg").write_bytes(JPEG)
        (figures / "D_x_SYM.jpeg").write_bytes(JPEG)
    return {
        "sample": "S",
        "version": "1.0.0",
        "work_dir": work,
        "failures": {},
        "table_path": table,
        "paraphase_json": paraphase,
        "figures_dir": figures,
        "generated": WHEN,
    }


def test_model_status_order_and_cn(tmp_path: Path) -> None:
    model = summary.build_model(**_run_files(tmp_path))
    status = {g["name"]: g["status"] for g in model["genes"]}
    assert status == {"C": "review", "A": "ok", "D_x": "ok"}
    # Review first, then OK (flagged genes before the rest).
    assert [g["name"] for g in model["genes"]] == ["C", "A", "D_x"]
    # A gene without haplotypes is listed, but does not raise the status:
    # most genes of a targeted run have none.
    assert model["no_haps"] == ["<b>B"]
    assert model["overall"] == "review" and not model["problems"]
    assert model["reviews"] == ["C: CN 2 resolved vs 3 expected"]
    by_name = {g["name"]: g for g in model["genes"]}
    assert by_name["A"]["reads"] == 15 and by_name["A"]["variants"] == 3
    # A figure named <gene>_<suffix> belongs to the gene D_x, not to "D".
    assert [p.name for p in by_name["D_x"]["figures"]] == ["D_x_SYM.jpeg"]


def test_every_flagged_row_has_a_reason_and_known_hit_matches(
    tmp_path: Path,
) -> None:
    files = _run_files(tmp_path)
    # Deletion in the table (ANNOVAR 300 G>-) equals VCF-style 299 AG>A.
    hit = ("G1", "chr1", 299, "AG", "A", "G1 del", "c.2del", "PMID:1")
    files["known_hits"] = [hit]
    files["hap_symbol"] = {"A_hap1": "G1", "A_hap2": "A"}
    model = summary.build_model(**files)
    flagged = {g["name"]: g["flagged"] for g in model["genes"]}
    assert all(r["_reasons"] for rows in flagged.values() for r in rows)
    reasons = sorted(r for row in flagged["A"] for r in row["_reasons"])
    assert reasons == [
        "ClinVar Pathogenic",
        "known variant G1 del (c.2del; PMID:1)",
    ]
    assert len(flagged["C"]) == 1 and len(flagged["D_x"]) == 0
    # The same hit on a haplotype of another gene symbol is ignored.
    files["hap_symbol"] = {"A_hap1": "OTHER"}
    model = summary.build_model(**files)
    assert sum(len(g["flagged"]) for g in model["genes"]) == 2


def test_page_escapes_embeds_figures_and_has_no_external_requests(
    tmp_path: Path,
) -> None:
    files = _run_files(tmp_path)
    files["failures"] = {"figures": "boom <script>x</script>\nline2"}
    out = summary.write_summary(
        tmp_path / "out" / summary.SUMMARY_NAME, **files
    )
    page = out.read_text(encoding="utf-8")
    assert "<b>B" not in page and "&lt;b&gt;B" in page
    assert "<i>" not in page and "<script>x" not in page
    assert "failed: figures: boom &lt;script&gt;" in page
    b64 = base64.b64encode(JPEG).decode()
    # One copy per figure file (two figures), the Download link reuses it.
    assert page.count(b64) == 2
    assert 'download="A.jpeg"' in page and 'download="D_x_SYM.jpeg"' in page
    assert " src=" not in page.replace("data-src=", "")
    assert "prefers-color-scheme" not in page and "color-scheme:light;" in page
    for rule in ("viewport", "REVEL > 0.7"):
        assert rule.replace(">", "&gt;") in page or rule in page
    for forbidden in ("http://", "https://", "<link", "@import"):
        assert forbidden not in page
    assert "RCCX" not in page
    assert not list(out.parent.glob("*.part"))


def test_ok_run_and_review_run(tmp_path: Path) -> None:
    files = _run_files(tmp_path, with_figures=False)
    audit = json.loads((files["work_dir"] / "split_audit.json").read_text())
    audit["genes"] = [g for g in audit["genes"] if g["gene"] in ("A", "D_x")]
    (files["work_dir"] / "split_audit.json").write_text(json.dumps(audit))
    model = summary.build_model(**files)
    assert model["overall"] == "ok" and not model["problems"]
    assert all(not g["figures"] for g in model["genes"])
    # A CN mismatch alone gives REVIEW, not PROBLEM.
    data = json.loads(files["paraphase_json"].read_text())
    data["A"]["total_cn"] = 4
    files["paraphase_json"].write_text(json.dumps(data))
    assert summary.build_model(**files)["overall"] == "review"


def test_failed_haplotype_marks_its_gene_a_problem(tmp_path: Path) -> None:
    files = _run_files(tmp_path)
    files["failures"] = {"haplotype A_hap2": "clair3 failed"}
    model = summary.build_model(**files)
    gene = next(g for g in model["genes"] if g["name"] == "A")
    assert gene["status"] == "problem"
    assert model["problems"][0].startswith("failed: haplotype A_hap2")


def test_rccx_section_and_missing_inputs(tmp_path: Path) -> None:
    files = _run_files(tmp_path)
    (files["work_dir"] / "rccx_calls.tsv").write_text(
        "hap\tallele_tag\trp_call\trp_reason\n"
        "H1\ta1\tproximal_type\tlow depth\n"
    )
    (files["figures_dir"] / "RCCX_module.jpeg").write_bytes(JPEG)
    page = summary.render_html(summary.build_model(**files))
    assert "RESEARCH USE ONLY" not in page and "low depth" not in page
    assert "RCCX_module.jpeg" in page and 'id="rccx-fig"' not in page
    assert 'id="frccx" alt="RCCX_module" loading="lazy" src="data:' in page
    # No table, no audit, no figures directory: still a page.
    (files["work_dir"] / "split_audit.json").unlink()
    files.update(table_path=None, figures_dir=None, paraphase_json=None)
    model = summary.build_model(**files)
    assert model["genes"] == [] and model["overall"] == "problem"
    assert "split_audit.json missing" in summary.render_html(model)


def test_flags_use_the_figure_p_lp_sets_and_both_clinvar_columns() -> None:
    from paralog2viewer import annotate

    assert annotate.PLP_CLINVAR == summary.PLP_CLINVAR
    sig = _row(ClinVar_SIG="Likely_pathogenic")
    assert summary.flag_reasons(sig) == ["ClinVar Likely_pathogenic"]
    assert not summary.flag_reasons(_row(ClinVar="Pathogenic, risk factor"))
    assert not summary.flag_reasons(_row(InterVar=_inter("Pathogenic?")))


def test_deletion_matches_annovar_row_with_leading_anchor_trimmed() -> None:
    # VCF REF AGA ALT A at 99 is the ANNOVAR row Start 100 Ref GA Alt -.
    assert summary._site("chr1", 99, "AGA", "A") == summary._site(
        "1", 100, "GA", "-"
    )
    assert summary._site("chr1", 99, "A", "AT") == summary._site(
        "1", 99, "-", "T"
    )


def _body(page: str) -> str:
    return page.split("</style>", 1)[1]


def test_red_marks_problems_and_is_absent_from_an_ok_page(
    tmp_path: Path,
) -> None:
    red = re.compile(r'class="[^"]*\bred\b')
    # Problem fixture: failure, CN mismatch, flagged rows, RCCX notice.
    bad = _run_files(tmp_path / "bad")
    bad["failures"] = {"figures": "boom"}
    (bad["work_dir"] / "rccx_calls.tsv").write_text("hap\trp_reason\nH\tx\n")
    page = summary.render_html(summary.build_model(**bad))
    assert len(red.findall(_body(page))) >= 6
    assert '<span class="t note" title="CN ' in page
    for text in ("failed: figures: boom", "2 / 3"):
        assert text in page
    assert "var(--red)" in page and "#D55E00" in page
    # All-OK fixture: no table, no flags, matching CN, no failures.
    ok = _run_files(tmp_path / "ok", with_figures=False)
    audit = json.loads((ok["work_dir"] / "split_audit.json").read_text())
    audit["genes"] = [g for g in audit["genes"] if g["gene"] in ("A", "D_x")]
    (ok["work_dir"] / "split_audit.json").write_text(json.dumps(audit))
    ok["table_path"] = None
    model = summary.build_model(**ok)
    assert model["overall"] == "ok"
    assert not red.findall(_body(summary.render_html(model)))


def test_plp_tier_takes_the_higher_of_clinvar_and_intervar() -> None:
    assert summary.plp_tier(_row(ClinVar="Pathogenic")) == "P"
    assert summary.plp_tier(_row(ClinVar="Likely_pathogenic")) == "LP"
    assert (
        summary.plp_tier(
            _row(ClinVar="Likely_pathogenic", InterVar=_inter("Pathogenic"))
        )
        == "P"
    )
    assert summary.plp_tier(_row(InterVar=_inter("Likely pathogenic"))) == "LP"
    assert summary.plp_tier(_row(ClinVar="Uncertain_significance")) == ""


def test_tiers_are_counted_per_locus_role(tmp_path: Path) -> None:
    files = _run_files(tmp_path)
    files["hap_role"] = {"A_hap1": "gene", "A_hap2": "pseudogene"}
    model = summary.build_model(**files)
    tiers = {g["name"]: dict(g["tiers"]) for g in model["genes"]}
    assert all(t in summary.TIERS for g in tiers.values() for t, _ in g)
    assert all(r in summary.ROLES for g in tiers.values() for _, r in g)
    # A haplotype without a role is unresolved, never silently a gene.
    assert ("P", "unresolved") in tiers["C"]
    page = summary.render_html(model)
    assert "OK</span>" not in page and "REVIEW" not in page
    assert 'data-tier="P"' in page and "Pathogenic: " in page


def test_rccx_section_says_it_is_the_in_house_pipeline(tmp_path: Path) -> None:
    files = _run_files(tmp_path)
    (files["work_dir"] / "rccx_calls.tsv").write_text("hap\trp_reason\nH\tx\n")
    page = summary.render_html(summary.build_model(**files))
    assert "IN-HOUSE PIPELINE" in page and "not the Paraphase result" in page


def test_same_variant_on_two_haplotypes_counts_twice(tmp_path: Path) -> None:
    files = _run_files(tmp_path)
    model = summary.build_model(**files)
    for gene in model["genes"]:
        flagged = gene["flagged"]
        counted = sum(
            n for (_t, r), n in gene["tiers"].items() if r != "copy number"
        )
        assert counted == len(flagged)


def test_cn_mismatch_counts_as_a_note(tmp_path: Path) -> None:
    model = summary.build_model(**_run_files(tmp_path))
    for gene in model["genes"]:
        cn = [n for n in gene["notes"] if n.startswith("CN ")]
        assert gene["tiers"].get(("note", "copy number"), 0) == len(cn)
    assert any(g["tiers"].get(("note", "copy number")) for g in model["genes"])
