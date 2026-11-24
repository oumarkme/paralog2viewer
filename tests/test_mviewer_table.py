# SPDX-License-Identifier: GPL-3.0-or-later
"""Tests for the 365-column MViewer table writer (synthetic data only)."""

import gzip
import json
from pathlib import Path

import pysam
import pytest

from paralog2viewer.mviewer_table import (
    DEFAULT_IGNORED_COLUMNS,
    _aachange_refgene,
    compare_tables,
    load_header,
    max_allele_freq,
    net_float,
    variant_type,
    write_mviewer_table,
)

HEADER = load_header()
COL = {name: i for i, name in enumerate(HEADER)}
VCF_TEXT = """##fileformat=VCFv4.2
##contig=<ID=chr10,length=1000>
##contig=<ID=chr2,length=1000>
##contig=<ID=chrX,length=1000>
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tgeneA_hap10\tgeneA_hap2
chr2\t100\t.\tG\tA\t.\tPASS\t.\tGT:DP\t1:12\t1:8
chr2\t200\trsX\tCT\tC\t50\tPASS\t.\tGT:DP\t0:9\t1:7
chr2\t300\t.\tA\t<DEL>\t.\tPASS\tEND=350;SVTYPE=DEL\tGT\t.\t1
chr10\t50\t.\tT\tTG\t.\tPASS\t.\tGT\t1\t0
chrX\t10\t.\tC\tT\t.\t.\t.\tGT\t0\t.
"""


def annovar(chrom: str, start: int, ref: str, alt: str, **extra: str) -> dict:
    """Return a native-named ANNOVAR dict with all header columns '.'."""
    cols = {name: "." for name in HEADER[5:258]}
    cols.update({"Func.refGene": "exonic", "Gene.refGene": "GENEA"})
    cols.update(
        {
            "Chr": chrom,
            "Start": str(start),
            "End": str(start + max(len(ref), 1) - 1),
            "Ref": ref,
            "Alt": alt,
            "1000g2015aug_all": "0.01",
            "1000g2015aug_afr": "0.02",
            "ExAC_ALL": "0.003",
            "esp6500siv2_all": "0.004",
            "avsnp147": "rs100",
            "cosmic70": "ID=COSM1,COSM2;OCCURENCE=1(kidney),2(thyroid)",
            "CLNSIG": "Benign",
            "CLNDN": "dis",
            "CLNDISDB": "MedGen:C1",
            "MutationAssessor_rankscore": "0.5",
            "gnomAD_genome_ALL": "3.086e-05",
            "gnomAD_exome_AFR": "1e-05",
            "gnomad41_genome_AF": "0.9",
            "gnomad41_exome_AF_afr": "2e-05",
            "gnomad41_exome_AF_remaining": "3e-05",
            "AAChange.refGene": "G:NM_1:exon1:c.1A>T:p.X1Y",
            "AAChange.refGeneWithVer": "G:NM_1.1:exon1:c.1A>T:p.X1Y",
            "GeneDetail.refGeneWithVer": "NM_1.1:exon1:r.spl",
            "rmsk": "Name=L1,Alu",
        }
    )
    cols.update(extra)
    return cols


def record(chrom: str, pos: int, ref: str, alt: str, **fields: object) -> dict:
    """Return a full annotation record."""
    rec = {
        "allele_key": f"{chrom}:{pos}:{ref}:{alt}",
        "chrom": chrom,
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "symbolic": False,
        "annovar": {},
        "intervar": " InterVar: Benign PVS1=0 ",
        "intervar_class": "Benign",
        "vep_variant": {},
        "vep_transcripts": [],
        "clinvar": "Likely benign",
        "clinvar_vid": "123",
        "taiwan_biobank": "0.00750",
        "hgmd_class": "DM",
        "genes": [],
        "panel": "P1,P2",
        "phenotypes": "Pheno",
        "inheritance": "AD",
        "transcript_source": "default_transcript",
        "consequence": "missense_variant,splice_region_variant",
        "gene": "GENEA",
        "transcript": "NM_000001.3",
        "exon": "4/10",
        "hgvsc": "NM_000001.3:c.647A>T",
        "hgvsp": "NP_000001.1:p.Leu216%3D",
        "am_class": "",
        "am_pathogenicity": "",
        "maxentscan": {},
    }
    rec.update(fields)
    return rec


def tx(symbol: str, consequence: str) -> dict:
    return {
        "feature": "NM_1.1",
        "feature_type": "Transcript",
        "symbol": symbol,
        "consequence": consequence,
    }


def write_inputs(tmp_path: Path) -> tuple[Path, Path]:
    """Write the synthetic merged VCF and records; return their paths."""
    vcf = tmp_path / "merged.vcf"
    vcf.write_text(VCF_TEXT)
    pysam.tabix_index(str(vcf), preset="vcf", force=True)
    records = [
        record(
            "chr2",
            100,
            "G",
            "A",
            annovar=annovar("chr2", 100, "G", "A"),
            vep_transcripts=[tx("GENEA", "missense_variant"), tx("GENEB", "")],
        ),
        record(
            "chr2",
            200,
            "CT",
            "C",
            annovar=annovar("chr2", 201, "T", "-"),
            transcript="",
        ),
        record(
            "chr2",
            300,
            "A",
            "<DEL>",
            symbolic=True,
            vep_transcripts=[
                tx("GENEA", "intron_variant,upstream_gene_variant"),
                tx("GENEA", "feature_truncation"),
            ],
        ),
        record("chr10", 50, "T", "TG"),
    ]
    path = tmp_path / "records.jsonl.gz"
    with gzip.open(path, "wt") as handle:
        for rec in records:
            handle.write(json.dumps(rec) + "\n")
    return vcf.with_name("merged.vcf.gz"), path


def read_table(path: Path) -> tuple[list[str], list[list[str]]]:
    lines = gzip.open(path, "rt").read().split("\n")
    assert lines[0] == "#keywords="
    return lines[1:3], [x.split("\t") for x in lines[3:] if x]


def test_header_resource_is_valid() -> None:
    assert len(HEADER) == 365
    assert HEADER[359] == "Genotype"


def test_rows_sort_and_genotype(tmp_path: Path) -> None:
    merged, records = write_inputs(tmp_path)
    out = tmp_path / "s.paralog.tsv.gz"
    n = write_mviewer_table(records, merged, "S1", "female", out)
    comment, rows = read_table(out)
    assert comment[0] == "#sex=Female"
    assert gzip.open(out, "rt").read().split("\n")[2].split("\t") == HEADER
    # chr2:100 x2 (hap2 before hap10 in natural order), chr2:200 hap2,
    # chr2:300 hap2 (symbolic), chr10:50 hap10; GT 0 and . give no row.
    key = [(r[0], r[COL["Genotype"]]) for r in rows]
    assert n == len(rows) == 5
    assert key == [
        ("2", "geneA_hap2"),
        ("2", "geneA_hap10"),
        ("2", "geneA_hap2"),
        ("2", "geneA_hap2"),
        ("10", "geneA_hap10"),
    ]
    assert all(len(r) == 365 for r in rows)


def test_explicit_allele_columns(tmp_path: Path) -> None:
    merged, records = write_inputs(tmp_path)
    out = tmp_path / "t.gz"
    write_mviewer_table(records, merged, "S1", "Male", out)
    _, rows = read_table(out)
    snv = rows[1]  # chr2:100, geneA_hap10
    get = lambda name: snv[COL[name]]  # noqa: E731
    assert [get(c) for c in ("Chr", "Start", "End", "Ref", "Alt")] == [
        "2", "100", "100", "G", "A",
    ]  # fmt: skip
    assert get("Func.refgene") == "exonic"
    # Otherinfo: zygosity, QUAL, depth, then the haplotype's own VCF cells.
    assert [get(f"Otherinfo{i}") for i in range(1, 14)] == [
        "hom", "0", ".", "chr2", "100", ".", "G", "A", ".", "PASS", ".",
        "GT:DP", "1:12",
    ]  # fmt: skip
    other = rows[0]  # same allele, other haplotype: its own sample string
    assert other[COL["Otherinfo13"]] == "1:8"
    # Legacy-named columns follow the ANNOVAR protocol columns.
    assert get("1000G_ALL") == "0.01"
    assert get("1000G_AFR") == "0.02"
    assert get("ExAC_Freq") == "0.003"
    assert get("ESP6500si_ALL") == "0.004"
    assert get("dbSNP") == get("dbSNP ID") == "rs100"
    assert get("COSMIC_ID") == "COSM1;COSM2"
    assert get("COSMIC_DIS") == "1(kidney);2(thyroid)"
    assert [get(c) for c in ("ClinVar_SIG", "ClinVar_DIS", "ClinVar_DB")] == [
        "Benign", "dis", "MedGen:C1",
    ]  # fmt: skip
    assert get("MutationAssessor_score_rankscore") == "0.5"
    # gnomAD legacy columns come from gnomad41, not from the v2 columns.
    assert get("gnomAD_genome_ALL") == "0.9"
    assert get("gnomAD_exome_AFR") == "2e-05"
    assert get("gnomAD_exome_OTH") == "3e-05"
    assert get("gnomAD_exome_SAS") == "."
    # Versioned AAChange and GeneDetail, joined by '/'.
    assert get("AAChange.refgene") == (
        "G:NM_1.1:exon1:c.1A>T:p.X1Y/NM_1.1:exon1:r.spl"
    )
    assert get("ClinVar_ID") == get("GWAS_DIS") == get("CG46") == ""
    # Rebuilt and lookup columns.
    assert get("Variant") == "G>A/A"
    assert get("Coordinate") == "100"
    assert get("Type") == "snv"
    assert get("Gene") == "GENEA;GENEB"
    assert get("Consequence") == "missense_variant;splice_region_variant"
    assert get("Panel") == "P1,P2"
    assert get("HGMD Variant Class") == "DM"
    assert get("ClinVar") == "Likely benign"
    assert get("ClinVar_VID") == "123"
    assert get("Taiwan Biobank") == "0.00750"
    assert get("InterVar") == " InterVar: Benign PVS1=0 "
    assert get("Interpretation Result") == "Benign"
    # Max over 1000G, ExAC, ESP, gnomAD and Taiwan Biobank; the gnomAD
    # columns carry the gnomad41 values (gnomad41_genome_AF is 0.9).
    assert get("Max Allele Freq") == "0.9"
    assert [get(c) for c in ("Transcript", "Exon", "Nucleotide")] == [
        "NM_000001.3", "exon4", "c.647A>T",
    ]  # fmt: skip
    assert get("AAChange") == "p.Leu216="
    assert get("Severity Score") == ""
    assert get("Methylated loci") == get("am_class") == ""


def test_indel_uses_annovar_start_and_vcf_id(tmp_path: Path) -> None:
    merged, records = write_inputs(tmp_path)
    out = tmp_path / "u.gz"
    write_mviewer_table(records, merged, "S1", "Male", out)
    _, rows = read_table(out)
    dele = rows[2]  # chr2:200, only geneA_hap2 carries it
    assert [dele[COL[c]] for c in ("Start", "Ref", "Alt")] == ["201", "T", "-"]
    assert dele[COL["Coordinate"]] == "201"
    assert dele[COL["Type"]] == "deletion"
    assert dele[COL["Variant"]] == "CT>C/C"
    assert dele[COL["Otherinfo2"]] == "50"
    assert dele[COL["Otherinfo6"]] == "rsX"
    assert dele[COL["Transcript"]] == dele[COL["Exon"]] == ""


def test_symbolic_row_keeps_only_allele_and_vep_fields(tmp_path: Path) -> None:
    merged, records = write_inputs(tmp_path)
    out = tmp_path / "v.gz"
    write_mviewer_table(records, merged, "S1", "Male", out)
    _, rows = read_table(out)
    sym = rows[3]
    filled = {HEADER[i]: v for i, v in enumerate(sym) if v}
    assert filled == {
        "Chr": "2",
        "Start": "300",
        "End": "350",
        "Ref": "A",
        "Alt": "<DEL>",
        "Consequence": "feature_truncation;intron_variant;"
        "upstream_gene_variant",
        "Gene": "GENEA",
        "Genotype": "geneA_hap2",
    }


def test_record_without_annovar_row_uses_vcf_values(tmp_path: Path) -> None:
    merged, records = write_inputs(tmp_path)
    out = tmp_path / "w.gz"
    write_mviewer_table(records, merged, "S1", "Male", out)
    _, rows = read_table(out)
    row = rows[4]  # chr10:50 T>TG, ANNOVAR dropped it
    assert [row[COL[c]] for c in ("Chr", "Start", "End", "Ref", "Alt")] == [
        "10", "50", "50", "T", "TG",
    ]  # fmt: skip
    assert row[COL["Func.refgene"]] == ""
    assert row[COL["Type"]] == "insertion"
    assert row[COL["dbSNP"]] == "."


def test_missing_record_and_bad_sex_raise(tmp_path: Path) -> None:
    merged, records = write_inputs(tmp_path)
    with pytest.raises(ValueError, match="sex"):
        write_mviewer_table(records, merged, "S1", "x", tmp_path / "a.gz")
    write_mviewer_table(records, merged, "S1", "unknown", tmp_path / "u.gz")
    with gzip.open(tmp_path / "u.gz", "rt") as fh:
        assert "#sex=Unknown" in fh.read(60)
    empty = tmp_path / "empty.jsonl.gz"
    with gzip.open(empty, "wt"):
        pass
    with pytest.raises(ValueError, match="no annotation record"):
        write_mviewer_table(empty, merged, "S1", "Male", tmp_path / "b.gz")


def test_max_allele_freq_and_formats() -> None:
    assert max_allele_freq({}) == "0"
    cells = {"1000G_ALL": ".", "gnomAD_exome_ALL": "3.086e-05", "x_AF": "0.9"}
    assert max_allele_freq(cells) == "3.086E-05"
    assert max_allele_freq({"ExAC_AFR": "1", "Taiwan Biobank": "0.5"}) == "1"
    assert net_float(0.0075) == "0.0075"
    assert variant_type("A", "AT") == "insertion"
    assert variant_type("AT", "A") == "deletion"
    assert variant_type("AT", "GC") == "other"


def test_aachange_refgene_rules() -> None:
    plain = {"aachange.refgene": "plain"}
    assert _aachange_refgene(plain) == "plain"
    assert _aachange_refgene({"aachange.refgenewithver": "."}) == "."
    detail = {"aachange.refgenewithver": ".", "genedetail.refgenewithver": "d"}
    assert _aachange_refgene(detail) == "d"
    both = {"aachange.refgenewithver": "a", "genedetail.refgenewithver": "d"}
    assert _aachange_refgene(both) == "a/d"


def test_compare_tables_documented_deviations(tmp_path: Path) -> None:
    merged, records = write_inputs(tmp_path)
    base = tmp_path / "base.gz"
    write_mviewer_table(records, merged, "S1", "Male", base)
    cells = [x.split("\t") for x in gzip.open(base, "rt").read().split("\n")]
    for row in cells[3:8]:
        for name in ("Otherinfo6", "Otherinfo13", "dbSNP", "dbSNP ID"):
            row[COL[name]] = "changed"
        if row[COL["Alt"]].startswith("<"):
            row[COL["Func.refgene"]] = "changed"
    new = tmp_path / "new.gz"
    with gzip.open(new, "wt") as handle:
        handle.write("\n".join("\t".join(c) for c in cells))
    diff = compare_tables(base, new)
    assert not any(v for k, v in diff.items() if not k.startswith("__"))
    # A symbolic row ignores only the two refgene columns, not others.
    assert "dbSNP" not in diff and "Func.refgene" in diff
    strict = compare_tables(
        base, new, ignore_columns=(), symbolic_ignore_columns=()
    )
    assert strict["Func.refgene"] == 1 and strict["dbSNP"] > 0


def test_compare_tables(tmp_path: Path) -> None:
    merged, records = write_inputs(tmp_path)
    base = tmp_path / "base.gz"
    write_mviewer_table(records, merged, "S1", "Male", base)
    new = tmp_path / "new.gz"
    lines = gzip.open(base, "rt").read().split("\n")
    cells = [x.split("\t") for x in lines]
    first = 3
    cells[first][COL["Max Allele Freq"]] = "999"  # ignored by default
    cells[first][COL["rmsk"]] = "Name=Alu,L1"  # same names, other order
    cells[first][COL["CADD_phred"]] = "5"  # real difference
    cells[first + 1][COL["Gene.refgene"]] = "OTHER"
    del cells[first + 2]  # one baseline row has no partner
    with gzip.open(new, "wt") as handle:
        handle.write("\n".join("\t".join(c) for c in cells))
    diff = compare_tables(base, new)
    assert diff["CADD_phred"] == 1
    assert diff["Gene.refgene"] == 1
    assert diff["rmsk"] == 0
    assert "Max Allele Freq" not in diff
    assert all(c not in diff for c in DEFAULT_IGNORED_COLUMNS)
    assert diff["__unmatched_baseline__"] == 1
    assert diff["__unmatched_new__"] == 0
    assert diff["__rows_baseline__"] == 5
    assert diff["__rows_new__"] == 4
    plain = compare_tables(base, new, sort_columns=())
    assert plain["rmsk"] == 1
    assert compare_tables(base, base, ignore_columns=())["Severity Score"] == 0
