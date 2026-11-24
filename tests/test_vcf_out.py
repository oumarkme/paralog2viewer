# SPDX-License-Identifier: GPL-3.0-or-later
"""Tests for the annotated VCF writer (synthetic data only)."""

import gzip
import json
import shutil
import subprocess
from pathlib import Path

import pysam
import pytest

from paralog2viewer.vcf_out import write_annotated_vcf

VCF_TEXT = """##fileformat=VCFv4.2
##contig=<ID=chr2,length=1000>
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHP1\tHP2
chr2\t100\t.\tG\tA\t.\tPASS\t.\tGT\t1\t0
chr2\t200\t.\tCT\tC\t.\tPASS\t.\tGT\t0\t1
chr2\t300\t.\tA\t<DEL>\t.\tPASS\tEND=350;SVTYPE=DEL\tGT\t1\t.
"""
HGMD_SENTINEL = "HGMD_SENTINEL_DM"
TWB_SENTINEL = "TWB_SENTINEL_0.0123"


def tx(feature: str, consequence: str, **plugins: str) -> dict:
    return {
        "feature": feature,
        "feature_type": "Transcript",
        "gene_id": "ENSG1",
        "symbol": "GENEA",
        "consequence": consequence,
        "impact": "MODERATE",
        "biotype": "protein_coding",
        "canonical": "YES",
        "mane_select": "",
        "mane_plus_clinical": "",
        "exon": "4/10",
        "intron": "",
        "hgvsc": f"{feature}:c.1A>T",
        "hgvsp": "NP_1.1:p.Leu2%3D",
        "cdna_position": "10",
        "cds_position": "1",
        "protein_position": "1",
        "amino_acids": "L/F",
        "codons": "Cta/Tta",
        "plugins": plugins,
    }


def record(pos: int, ref: str, alt: str, **fields: object) -> dict:
    rec = {
        "allele_key": f"chr2:{pos}:{ref}:{alt}",
        "chrom": "chr2",
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "symbolic": False,
        "annovar": {},
        "intervar": "",
        "intervar_class": "",
        "vep_transcripts": [],
        "clinvar": "",
        "taiwan_biobank": TWB_SENTINEL,
        "hgmd_class": HGMD_SENTINEL,
    }
    rec.update(fields)
    return rec


@pytest.fixture
def inputs(tmp_path: Path) -> tuple[Path, Path]:
    vcf = tmp_path / "merged.vcf"
    vcf.write_text(VCF_TEXT)
    pysam.tabix_index(str(vcf), preset="vcf", force=True)
    records = [
        record(
            100,
            "G",
            "A",
            intervar_class="Likely pathogenic",
            clinvar="Pathogenic/Likely pathogenic",
            annovar={"Gene.refGene": "GENEA", "CADD_phred": "27.6"},
            vep_transcripts=[
                tx(
                    "NM_1.1",
                    "missense_variant,splice_region_variant",
                    MaxEntScan_alt="5.5",
                    am_class="likely_pathogenic",
                    am_pathogenicity="0.97",
                ),
                tx("NM_2.1", "intron_variant"),
            ],
        ),
        record(200, "CT", "C", annovar={"CLNSIG": "Benign", "avsnp147": "."}),
        record(300, "A", "<DEL>", symbolic=True),
    ]
    path = tmp_path / "records.jsonl.gz"
    with gzip.open(path, "wt") as handle:
        for rec in records:
            handle.write(json.dumps(rec) + "\n")
    return vcf.with_name("merged.vcf.gz"), path


def test_default_mode_info_fields(
    inputs: tuple[Path, Path], tmp_path: Path
) -> None:
    merged, records = inputs
    out = write_annotated_vcf(merged, records, tmp_path / "s.paralog.vcf.gz")
    assert Path(str(out) + ".tbi").exists()
    with pysam.VariantFile(str(out)) as vcf:
        header = vcf.header
        assert list(header.samples) == ["HP1", "HP2"]
        description = header.info["CSQ"].description
        subfields = description.split("Format: ")[1].split("|")
        assert {"MaxEntScan_alt", "am_class", "am_pathogenicity"} <= set(
            subfields
        )
        recs = list(vcf)
        assert [r.pos for r in recs] == [100, 200, 300]
        assert [r.samples["HP1"]["GT"] for r in recs] == [(1,), (0,), (1,)]
        snv, dele, sym = recs
        entries = [
            dict(zip(subfields, e.split("|"), strict=True))
            for e in snv.info["CSQ"]
        ]
        assert len(entries) == 2
        assert entries[0]["Allele"] == "A"
        assert (
            entries[0]["Consequence"]
            == "missense_variant&splice_region_variant"
        )
        assert entries[0]["MaxEntScan_alt"] == "5.5"
        assert entries[0]["am_class"] == "likely_pathogenic"
        assert entries[0]["am_pathogenicity"] == "0.97"
        assert entries[0]["HGVSp"] == "NP_1.1:p.Leu2%3D"
        assert entries[1]["MaxEntScan_alt"] == ""
        assert snv.info["INTERVAR"] == ("Likely_pathogenic",)
        assert snv.info["CLNSIG"] == ("Pathogenic/Likely_pathogenic",)
        assert dele.info["CLNSIG"] == ("Benign",)  # ANNOVAR fallback
        assert "INTERVAR" not in dele.info
        assert "CSQ" not in sym.info
        assert sym.info["SVTYPE"] == "DEL"
        assert "ANN_CADD_phred" not in header.info
    text = gzip.open(out, "rt").read()
    assert HGMD_SENTINEL not in text
    assert TWB_SENTINEL not in text


def test_full_mode_adds_annovar_but_never_licensed_values(
    inputs: tuple[Path, Path], tmp_path: Path
) -> None:
    merged, records = inputs
    out = write_annotated_vcf(
        merged, records, tmp_path / "f.vcf.gz", mode="full"
    )
    with pysam.VariantFile(str(out)) as vcf:
        snv, dele, _ = list(vcf)
        assert snv.info["ANN_CADD_phred"] == ("27.6",)
        assert snv.info["ANN_Gene_refGene"] == ("GENEA",)
        assert "ANN_avsnp147" not in dele.info  # "." counts as missing
    text = gzip.open(out, "rt").read()
    assert HGMD_SENTINEL not in text
    assert TWB_SENTINEL not in text
    assert "Taiwan" not in text
    assert "HGMD" not in text


def test_bad_mode_and_missing_record(
    inputs: tuple[Path, Path], tmp_path: Path
) -> None:
    merged, records = inputs
    with pytest.raises(ValueError, match="mode"):
        write_annotated_vcf(merged, records, tmp_path / "x.vcf.gz", "all")
    empty = tmp_path / "empty.jsonl.gz"
    with gzip.open(empty, "wt"):
        pass
    with pytest.raises(ValueError, match="no annotation record"):
        write_annotated_vcf(merged, empty, tmp_path / "y.vcf.gz")


def test_bcftools_reads_the_output(
    inputs: tuple[Path, Path], tmp_path: Path
) -> None:
    bcftools = shutil.which("bcftools")
    if bcftools is None:
        pytest.skip("bcftools not available")
    merged, records = inputs
    out = write_annotated_vcf(merged, records, tmp_path / "b.vcf.gz")
    result = subprocess.run(
        [bcftools, "query", "-f", "%POS\\t%INFO/INTERVAR\\n", str(out)],
        capture_output=True,
        text=True,
        check=True,
    )
    assert result.stdout.splitlines() == [
        "100\tLikely_pathogenic",
        "200\t.",
        "300\t.",
    ]
    region = subprocess.run(
        [bcftools, "view", "-H", "-r", "chr2:200-200", str(out)],
        capture_output=True,
        text=True,
        check=True,
    )
    assert region.stdout.startswith("chr2\t200\t")
