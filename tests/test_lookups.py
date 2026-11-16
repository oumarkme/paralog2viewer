# SPDX-License-Identifier: GPL-3.0-or-later
"""Tests for paralog2viewer.lookups with tiny synthetic tables."""

from __future__ import annotations

from pathlib import Path

import pysam
import pytest

from paralog2viewer.datalayout import DataLayout
from paralog2viewer.lookups import (
    CLINVAR_COLUMNS,
    LookupSchemaError,
    LookupTables,
    VariantTable,
    load_default_transcripts,
    minimal_allele,
    parse_omim_phenotypes,
)

CLINVAR = """\
Chr\tPos\tRef\tAlt\tClinVar\tClinVar_VID
6\t100\tG\tA\tPathogenic\t111
6\t200\tCT\tC\tLikely pathogenic\t222
6\t301\tT\t-\tBenign\t333
6\t400\t-\tGG\tUncertain significance\t444
"""
TAIWAN = """\
Chr\tPos\tRef\tFreq.
chr6\t100\tG\tG:0.9925 A:0.0075 *:0.0010
chr6\t200\tCT\tCT:0.9000 C:0.1000
chr6\t300\tA\tA:0.9990 *:0.0010
chr6\t400\tAG\tAG:0.5 A:0.5
"""
HGMD = """\
##fileformat=VCFv4.2
##contig=<ID=chr6>
##INFO=<ID=CLASS,Number=1,Type=String,Description="Class">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr6\t100\tCM1\tG\tA\t.\t.\tCLASS=DM
chr6\t200\tCM2\tCT\tC\t.\t.\tCLASS=DM?
"""
PANELS = """\
Gene\tPanel
CYP21A2\tPanel B
CYP21A2\tPanel A
TNXB\tPanel A
"""
PHENOS = (
    "Adrenal hyperplasia, 201910 (3), Autosomal recessive; "
    "Other, 100000 (3), Autosomal dominant, Autosomal recessive"
)
XL = "A, 1 (3), X-linked recessive; B, 2 (3), X-linked dominant"
GENEMAP2 = f"""\
# Copyright
# Chromosome\tApproved Gene Symbol\tGene Symbols\tPhenotypes
chr6\tCYP21A2\tCYP21A2, CYP21\t{PHENOS}
chr6\tNOPHENO\tNOPHENO\t
chrX\tIKBKG\tIKBKG, NEMO\t{XL}
chr1\tGBA1\tGBA1, GBA\tG, 3 (3), Multifactorial, Autosomal recessive
"""
DEFAULTS = """\
Gene\tTranscript
CYP21A2\tNM_000500.9
TNXB\tNM_019105.8\tENST00000375244.7
"""


def _index(path: Path, **kwargs: object) -> None:
    """Compress and tabix-index a plain-text table (writes path + .gz)."""
    pysam.tabix_index(str(path), force=True, **kwargs)  # type: ignore[arg-type]


def make_profile(profile: Path) -> None:
    """Write all profile tables of the data layout into ``profile``."""
    profile.mkdir(parents=True)
    tsv = {"seq_col": 0, "start_col": 1, "end_col": 1, "line_skip": 1}
    (profile / "ClinVar.tsv").write_text(CLINVAR)
    _index(profile / "ClinVar.tsv", **tsv)
    (profile / "TaiwanBioBank_GRCh38.freq-only.tsv").write_text(TAIWAN)
    _index(profile / "TaiwanBioBank_GRCh38.freq-only.tsv", **tsv)
    (profile / "HGMD_Pro_hg38.vcf").write_text(HGMD)
    _index(profile / "HGMD_Pro_hg38.vcf", preset="vcf")
    (profile / "PanelList.tsv").write_text(PANELS)
    (profile / "genemap2.txt").write_text(GENEMAP2)
    (profile / "DefaultTranscript.tsv").write_text(DEFAULTS)


@pytest.fixture
def tables(tmp_path: Path):
    """Open lookup tables on the synthetic profile."""
    make_profile(tmp_path / "profile")
    with LookupTables(DataLayout(tmp_path)) as opened:
        yield opened


def test_minimal_allele_vcf_and_annovar_forms_agree() -> None:
    assert minimal_allele(300, "AT", "A") == minimal_allele(301, "T", "-")
    assert minimal_allele(400, "A", "AGG") == minimal_allele(400, "-", "GG")
    assert minimal_allele(100, "G", "A") == (100, "G", "A")
    assert minimal_allele(10, "ACGT", "AGGT") == (11, "C", "G")


def test_variant_fields_snv_all_tables(tables: LookupTables) -> None:
    got = tables.variant_fields("chr6", 100, "G", "A")
    assert got == {
        "clinvar": "Pathogenic",
        "clinvar_vid": "111",
        "taiwan_biobank": "0.00750",
        "hgmd_class": "DM",
    }


def test_variant_fields_vcf_and_annovar_style_indels(
    tables: LookupTables,
) -> None:
    assert tables.variant_fields("chr6", 200, "CT", "C")["clinvar"] == (
        "Likely pathogenic"
    )
    assert tables.variant_fields("chr6", 300, "AT", "A")["clinvar_vid"] == "333"
    assert (
        tables.variant_fields("chr6", 400, "A", "AGG")["clinvar_vid"] == "444"
    )
    assert tables.variant_fields("chr6", 200, "CT", "C")["hgmd_class"] == "DM?"


def test_taiwan_biobank_exact_site_rule(tables: LookupTables) -> None:
    taiwan = tables.taiwan
    assert taiwan.query("chr6", 200, "CT", "C") == "0.10000"
    # ALT absent from an existing row: the whole Freq. text, as in the
    # baseline.
    assert taiwan.query("chr6", 100, "G", "T") == "G:0.9925 A:0.0075 *:0.0010"
    # The star allele is plain text of the row; a REF mismatch finds no row.
    assert taiwan.query("chr6", 300, "A", "T") == "A:0.9990 *:0.0010"
    assert taiwan.query("chr6", 400, "A", "T") == ""
    # Differently anchored query (ANNOVAR-style deletion): no normalisation,
    # so no match; the VCF-style query matches.
    assert taiwan.query("chr6", 401, "G", "-") == ""
    assert taiwan.query("chr6", 400, "AG", "A") == "0.50000"


def test_variant_fields_same_position_other_allele_is_empty(
    tables: LookupTables,
) -> None:
    got = tables.variant_fields("chr6", 100, "G", "T")
    assert got["clinvar"] == got["clinvar_vid"] == got["hgmd_class"] == ""
    # Taiwan Biobank reports the whole Freq. text for an absent ALT.
    assert got["taiwan_biobank"].startswith("G:0.9925")


def test_variant_fields_unknown_contig_is_empty(tables: LookupTables) -> None:
    assert set(tables.variant_fields("chr9", 100, "G", "A").values()) == {""}


def test_clinvar_ignores_non_grch38_rows(tmp_path: Path) -> None:
    table = tmp_path / "ClinVar.tsv"
    table.write_text(
        "Chromosome\tPositionVCF\tReferenceAlleleVCF\tAlternateAlleleVCF"
        "\tClinicalSignificance\tAssembly\tVariationID\n"
        "6\t100\tG\tA\tBenign\tGRCh37\t1\n"
        "6\t100\tG\tA\tPathogenic\tGRCh38\t2\n"
        "6\t500\tG\tA\tBenign\tGRCh37\t3\n"
    )
    _index(table, seq_col=0, start_col=1, end_col=1, line_skip=1)
    clinvar = VariantTable(tmp_path / "ClinVar.tsv.gz", CLINVAR_COLUMNS)
    assert clinvar.query("chr6", 100, "G", "A") == [
        {"clinvar": "Pathogenic", "clinvar_vid": "2"}
    ]
    assert clinvar.query("chr6", 500, "G", "A") == []
    clinvar.close()


def test_gene_info(tables: LookupTables) -> None:
    info = tables.genes.gene_info(["CYP21A2", "TNXB", "UNKNOWN"])
    # Sorted within a gene; not deduplicated across genes.
    assert info["panel"] == "Panel A,Panel B,Panel A"
    assert info["phenotypes"] == PHENOS
    assert info["inheritance"] == "AD, AR"
    assert tables.genes.gene_info([]) == {
        "panel": "",
        "phenotypes": "",
        "inheritance": "",
    }


def test_gene_info_alias_x_linked_and_joins(tables: LookupTables) -> None:
    # GBA is only an alias in genemap2: it falls back to the GBA1 row.
    gba = tables.genes.gene_info(["GBA"])
    assert gba["phenotypes"] == "G, 3 (3), Multifactorial, Autosomal recessive"
    assert gba["inheritance"] == "AR, Multi"
    # XL accompanies XLR/XLD; order XLR, XLD, XL; genes join with ", ".
    both = tables.genes.gene_info(["IKBKG", "GBA"])
    assert both["inheritance"] == "AR, Multi, XLR, XLD, XL"
    assert both["phenotypes"] == (
        f"{XL}, G, 3 (3), Multifactorial, Autosomal recessive"
    )
    # An approved symbol with no phenotype is not replaced by an alias.
    assert tables.genes.gene_info(["NOPHENO"])["phenotypes"] == ""


def test_default_transcripts_ignore_version_and_header(
    tmp_path: Path,
) -> None:
    path = tmp_path / "d.tsv"
    path.write_text(DEFAULTS)
    assert load_default_transcripts(path) == {
        "CYP21A2": {"NM_000500"},
        "TNXB": {"NM_019105", "ENST00000375244"},
    }


def test_parse_omim_phenotypes_unparsed_entry_kept() -> None:
    out = parse_omim_phenotypes("Weird text; Name, 123456 (3)")
    assert out[0]["name"] == "Weird text"
    assert out[1] == {
        "text": "Name, 123456 (3)",
        "name": "Name",
        "mim": "123456",
        "mapping_key": "3",
        "inheritance": "",
    }


def test_schema_error_when_alleles_missing(tmp_path: Path) -> None:
    plain = tmp_path / "bad.tsv"
    plain.write_text("Chr\tPos\tValue\n1\t5\tx\n")
    _index(plain, seq_col=0, start_col=1, end_col=1, line_skip=1)
    with pytest.raises(LookupSchemaError):
        VariantTable(
            tmp_path / "bad.tsv.gz", {"ref": ("ref",), "alt": ("alt",)}
        )


def test_bare_chrom_strips_only_a_chr_prefix() -> None:
    from paralog2viewer.datalayout import bare_chrom

    assert [bare_chrom(c) for c in ("chr1", "1", "chrX", "chr", "MT")] == [
        "1",
        "1",
        "X",
        "chr",
        "MT",
    ]
