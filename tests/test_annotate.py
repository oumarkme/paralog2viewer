# SPDX-License-Identifier: GPL-3.0-or-later
"""Tests for paralog2viewer.annotate with mocked VEP, ANNOVAR and InterVar."""

from __future__ import annotations

import json
import threading
from pathlib import Path

import pysam
import pytest

from paralog2viewer import annotate, proc
from paralog2viewer.datalayout import DataLayout
from test_lookups import make_profile

MERGED = """\
##fileformat=VCFv4.2
##contig=<ID=chr6,length=170805979>
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Length">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHP1\tHP2
chr6\t100\t.\tG\tA,T\t.\tPASS\t.\tGT\t1\t2
chr6\t100\t.\tG\tA\t.\tPASS\t.\tGT\t1\t0
chr6\t200\t.\tCT\tC\t.\tPASS\t.\tGT\t0\t1
chr6\t500\t.\tA\t<DEL>\t.\tPASS\tEND=600;SVTYPE=DEL;SVLEN=-100\tGT\t1\t.
chr6\t700\t.\tC\tC[chr6:900[\t.\tPASS\tSVTYPE=BND\tGT\t1\t.
"""
KEY_SNV = "chr6:100:G:A"
KEY_SNV2 = "chr6:100:G:T"
KEY_DEL = "chr6:200:CT:C"
KEY_SYM = "chr6:500:A:<DEL>"
KEY_BND = "chr6:700:C:C[chr6:900["
ALL_KEYS = [KEY_SNV, KEY_SNV2, KEY_DEL, KEY_SYM, KEY_BND]
VEP_HEADER = (
    "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\t"
    "Consequence\tcDNA_position\tCDS_position\tProtein_position\t"
    "Amino_acids\tCodons\tExisting_variation\tExtra"
)


def write_merged(path: Path, text: str = MERGED) -> Path:
    """Write and index a merged VCF; return the .vcf.gz path."""
    path.write_text(text)
    pysam.tabix_index(str(path), preset="vcf", force=True)
    return Path(str(path) + ".gz")


def vcf_ids(path: Path) -> list[str]:
    """Return the ID column of a plain VCF."""
    return [
        line.split("\t")[2]
        for line in path.read_text().splitlines()
        if not line.startswith("#")
    ]


class FakeAnnotators:
    """Replacement for proc.run that writes plausible tool outputs."""

    def __init__(self, echo_id: bool = True, fail: str = "") -> None:
        self.echo_id = echo_id
        self.fail = fail
        self.calls: list[list[str]] = []
        self.order: list[str] = []
        self.barrier = threading.Barrier(2, timeout=20)

    def __call__(self, cmd, *, cwd=None, env=None, stdout=None, stderr=None):
        argv = [str(c) for c in cmd]
        self.calls.append(argv)
        if argv[0] == "vep":
            self._vep(argv)
        elif argv[0] == "perl":
            self._annovar(argv)
        elif argv[0] == "python3":
            self._intervar(argv)
        else:  # pragma: no cover
            raise AssertionError(argv)

    @staticmethod
    def _sites(path: str) -> list[list[str]]:
        return [
            line.split("\t")
            for line in Path(path).read_text().splitlines()
            if not line.startswith("#")
        ]

    def _vep(self, argv: list[str]) -> None:
        src = argv[argv.index("--input_file") + 1]
        dst = Path(argv[argv.index("--output_file") + 1])
        small = Path(src).name == "small.vcf"
        if small:
            self.order.append("vep_small")
            self.barrier.wait()
        lines = ["## fake VEP", VEP_HEADER]
        for _chrom, pos, ident, ref, alt, *_ in self._sites(src):
            shown = ident if self.echo_id else "."
            start = int(pos)
            vep_alt = alt
            if small and ref[0] == alt[0]:
                ref, vep_alt, start = ref[1:], alt[1:], start + 1
            loc = (
                f"6:{start}"
                if len(ref) <= 1
                else f"6:{start}-{start + len(ref) - 1}"
            )
            if small and not ref:
                loc = f"6:{start}-{start - 1}"
            vep_alt = vep_alt or "-"
            if not small:
                lines.append(
                    f"{shown}\t{loc}\t{alt}\t-\t-\t-\tintergenic_variant"
                    "\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER"
                )
                continue
            common = "SYMBOL=CYP21A2;IMPACT=MODERATE;EXON=2/10"
            lines.append(
                f"{shown}\t{loc}\t{vep_alt}\tENSG1\tENST00000000001\tTranscript"
                f"\tmissense_variant\t1\t1\t1\tA/V\tGcc/Gtc\trs1\t{common};"
                "CANONICAL=YES;MANE_SELECT=NM_000500.9;"
                "HGVSc=ENST00000000001.1:c.1A>G;am_class=likely_pathogenic"
            )
            lines.append(
                f"{shown}\t{loc}\t{vep_alt}\tENSG1\tNM_000500.9\tTranscript"
                f"\tmissense_variant\t1\t1\t1\tA/V\tGcc/Gtc\trs1\t{common};"
                "MANE_SELECT=NM_000500.9;HGVSc=NM_000500.9:c.1A>G;"
                "HGVSp=NP_1.1:p.Ala1Val;am_pathogenicity=0.9;"
                "MaxEntScan_diff=0.1;VARIANT_CLASS=SNV;gnomADg_AF=0.01"
            )
        dst.write_text("\n".join(lines) + "\n")

    def _annovar(self, argv: list[str]) -> None:
        self.order.append("annovar")
        self.barrier.wait()
        if self.fail == "annovar":
            raise proc.ToolError("annovar failed")
        prefix = argv[argv.index("-out") + 1]
        header = ["Chr", "Start", "End", "Ref", "Alt", "Func.refGene",
                  "Gene.refGene", "ClinVar_SIG"]  # fmt: skip
        header += [f"Otherinfo{i}" for i in range(1, 12)]
        rows, avin = ["\t".join(header)], []
        for chrom, pos, ident, ref, alt, qual, flt, info in self._sites(
            argv[2]
        ):
            start, a_ref, a_alt = int(pos), ref, alt
            if ref[0] == alt[0]:
                a_ref, a_alt, start = ref[1:] or "-", alt[1:] or "-", start + 1
            sig = "Pathogenic" if alt == "A" else "."
            cells = ["6", str(start), str(start + len(a_ref) - 1), a_ref,
                     a_alt, "exonic", "CYP21A2(dist=0);TNXB", sig,
                     "het", "30", "10", chrom, pos, ident, ref, alt, qual,
                     flt, info]  # fmt: skip
            rows.append("\t".join(cells))
            avin.append("\t".join(cells[:5]))
        Path(prefix + ".hg38_multianno.txt").write_text("\n".join(rows) + "\n")
        Path(prefix + ".avinput").write_text("\n".join(avin) + "\n")

    def _intervar(self, argv: list[str]) -> None:
        assert "annovar" in self.order, "InterVar must run after ANNOVAR"
        self.order.append("intervar")
        out = next(a for a in argv if a.startswith("--output=")).split("=")[1]
        multianno = Path(out + ".hg38_multianno.txt").read_text().splitlines()
        header = "#Chr\tStart\tEnd\tRef\tAlt\tRef.Gene"
        lines = [header + "\t InterVar: InterVar and Evidence "]
        for row in multianno[1:]:
            c = row.split("\t")
            verdict = "Likely pathogenic" if c[4] == "A" else "Benign"
            lines.append(
                "\t".join(c[:5])
                + f"\tCYP21A2\t InterVar: {verdict} PVS1=0 PS=[0, 0]"
            )
        Path(out + ".hg38_multianno.txt.intervar").write_text(
            "\n".join(lines) + "\n"
        )


@pytest.fixture
def layout(tmp_path: Path) -> DataLayout:
    """A data directory with a manifest and the synthetic profile."""
    root = tmp_path / "data"
    root.mkdir()
    (root / "manifest.json").write_text(json.dumps({"build": "GRCh38"}))
    make_profile(root / "profile")
    return DataLayout(root)


@pytest.fixture
def fake(monkeypatch: pytest.MonkeyPatch) -> FakeAnnotators:
    """Patch proc.run with the fake annotators."""
    tools = FakeAnnotators()
    monkeypatch.setattr(proc, "run", tools)
    return tools


def run_all(
    tmp_path: Path,
    layout: DataLayout,
    name: str = "w",
    merged: Path | None = None,
) -> tuple[Path, Path]:
    """Run annotate_alleles; return the records path and the merged VCF."""
    merged = merged or write_merged(tmp_path / "merged.vcf")
    out = annotate.annotate_alleles(
        merged,
        tmp_path / "ref.fa",
        layout,
        tmp_path / name,
        tmp_path / "cache" / "annot.sqlite",
        threads=2,
    )
    return out, merged


def test_sites_only_vcf_splits_small_and_symbolic(tmp_path: Path) -> None:
    merged = write_merged(tmp_path / "merged.vcf")
    sites = annotate.sites_only_vcf(merged, tmp_path / "s")
    assert vcf_ids(sites.small) == [KEY_SNV, KEY_SNV2, KEY_DEL]
    assert vcf_ids(sites.symbolic) == [KEY_SYM, KEY_BND]
    text = sites.symbolic.read_text()
    assert "##FORMAT" not in text and "HP1" not in text
    assert "\tEND=600;SVTYPE=DEL;SVLEN=-100" in text
    only = annotate.sites_only_vcf(merged, tmp_path / "o", only={KEY_DEL})
    assert vcf_ids(only.small) == [KEY_DEL] and only.symbolic is None


def test_annotate_alleles_end_to_end(
    tmp_path: Path, layout: DataLayout, fake: FakeAnnotators
) -> None:
    out, _ = run_all(tmp_path, layout)
    records = {r["allele_key"]: r for r in annotate.read_records(out)}
    assert list(records) == ALL_KEYS
    assert all(tuple(r) == annotate.RECORD_FIELDS for r in records.values())

    snv = records[KEY_SNV]
    assert snv["annovar"]["Chr"] == "6" and "Otherinfo1" not in snv["annovar"]
    assert snv["intervar_class"] == "Likely pathogenic"
    assert snv["intervar"].startswith(" InterVar: Likely pathogenic PVS1=0")
    assert (snv["clinvar"], snv["clinvar_vid"]) == ("Pathogenic", "111")
    assert snv["taiwan_biobank"] == "0.00750" and snv["hgmd_class"] == "DM"
    # DefaultTranscript.tsv names NM_000500, so the RefSeq row wins over
    # the canonical Ensembl row.
    assert snv["transcript_source"] == "default_transcript"
    assert snv["transcript"] == "NM_000500.9"
    assert snv["hgvsc"] == "NM_000500.9:c.1A>G" and snv["exon"] == "2/10"
    assert snv["gene"] == "CYP21A2" and snv["consequence"] == "missense_variant"
    assert snv["am_pathogenicity"] == "0.9" and snv["am_class"] == ""
    assert snv["maxentscan"] == {"MaxEntScan_diff": "0.1"}
    assert snv["vep_variant"]["VARIANT_CLASS"] == "SNV"
    assert snv["vep_variant"]["gnomADg_AF"] == "0.01"
    assert len(snv["vep_transcripts"]) == 2
    # The fake ANNOVAR names CYP21A2 and TNXB, so both genes are looked up.
    assert snv["genes"] == ["CYP21A2", "TNXB"]
    assert snv["panel"] == "Panel A,Panel B,Panel A"
    assert snv["inheritance"] == "AD, AR"

    # Indel: joined through the original VCF columns despite ANNOVAR's
    # rewritten coordinates, and through the VEP trimming of the allele.
    assert records[KEY_DEL]["annovar"]["Start"] == "201"
    assert records[KEY_DEL]["clinvar"] == "Likely pathogenic"
    assert records[KEY_DEL]["transcript"] == "NM_000500.9"

    sym = records[KEY_SYM]
    assert sym["annovar"] == {} and sym["intervar"] == ""
    assert (
        sym["symbolic"] is True and sym["consequence"] == "intergenic_variant"
    )
    assert sym["transcript_source"] == "first_row" and sym["gene"] == ""

    assert annotate.pathogenic_allele_keys(out) == {KEY_SNV, KEY_DEL}


def test_commands_use_layout_paths_and_run_concurrently(
    tmp_path: Path, layout: DataLayout, fake: FakeAnnotators
) -> None:
    run_all(tmp_path, layout)  # the barrier fails if VEP/ANNOVAR are serial
    vep = next(c for c in fake.calls if c[0] == "vep" and "--hgvs" in c)
    joined = " ".join(vep)
    for option in ("--merged", "--use_given_ref", "--offline", "--tab",
                   "--cache_version 116", "--fork 2",
                   f"--dir {layout.vep_cache}",
                   f"--dir_plugins {layout.vep_plugins}",
                   f"MaxEntScan,{layout.maxentscan}",
                   f"AlphaMissense,file={layout.alphamissense}"):  # fmt: skip
        assert option in joined
    annovar = next(c for c in fake.calls if c[0] == "perl")
    assert str(layout.annovar_humandb) in annovar and "-vcfinput" in annovar
    intervar = next(c for c in fake.calls if c[0] == "python3")
    assert f"--database_intervar={layout.intervardb}" in intervar
    assert "--skip_annovar" in intervar
    symbolic = [c for c in fake.calls if c[0] == "vep" and "--hgvs" not in c]
    assert len(symbolic) == 1 and "--symbol" in symbolic[0]
    assert fake.order.index("annovar") < fake.order.index("intervar")


def test_join_falls_back_to_location_and_allele(
    tmp_path: Path, layout: DataLayout, fake: FakeAnnotators
) -> None:
    fake.echo_id = False
    out, _ = run_all(tmp_path, layout)
    records = {r["allele_key"]: r for r in annotate.read_records(out)}
    assert records[KEY_SNV]["transcript"] == "NM_000500.9"
    assert records[KEY_DEL]["transcript"] == "NM_000500.9"
    assert records[KEY_SYM]["consequence"] == "intergenic_variant"


def test_cache_hits_skip_annotators_and_manifest_change_invalidates(
    tmp_path: Path, layout: DataLayout, fake: FakeAnnotators
) -> None:
    first, merged = run_all(tmp_path, layout, "w1")
    assert fake.calls
    fake.calls.clear()
    second, _ = run_all(tmp_path, layout, "w2")
    assert not fake.calls
    assert list(annotate.read_records(first)) == list(
        annotate.read_records(second)
    )

    # One new allele: only it reaches the annotators.
    new = "chr6\t300\t.\tA\tG\t.\tPASS\t.\tGT\t1\t1\n"
    extra = MERGED.replace("chr6\t500", new + "chr6\t500")
    grown = write_merged(tmp_path / "grown.vcf", extra)
    fake.barrier = threading.Barrier(2, timeout=20)
    run_all(tmp_path, layout, "w3", merged=grown)
    assert vcf_ids(tmp_path / "w3" / "annotation" / "small.vcf") == [
        "chr6:300:A:G"
    ]

    fake.calls.clear()
    fake.barrier = threading.Barrier(2, timeout=20)
    (layout.root / "manifest.json").write_text(json.dumps({"build": "new"}))
    run_all(tmp_path, layout, "w4")
    assert fake.calls


def test_tool_failure_raises_and_nothing_is_cached(
    tmp_path: Path, layout: DataLayout, fake: FakeAnnotators
) -> None:
    fake.fail = "annovar"
    with pytest.raises(proc.ToolError):
        run_all(tmp_path, layout)
    cache = annotate.RecordCache(
        tmp_path / "cache" / "annot.sqlite",
        annotate.manifest_checksum(layout),
    )
    assert cache.get_many(ALL_KEYS) == {}


def _row(ident: str, location: str, allele: str) -> dict:
    return {
        annotate.VEP_KEY_COLUMN: ident,
        "Location": location,
        "Allele": allele,
        "_extra": {},
    }


def test_group_vep_rows_by_id_and_by_trimmed_location() -> None:
    alleles = [
        annotate.Allele("chr6", 100, "G", "A", False),
        annotate.Allele("chr6", 200, "CT", "C", False),
        annotate.Allele("chr6", 300, "A", "AGG", False),
    ]
    rows = [
        _row(alleles[0].key, "6:100", "A"),
        _row(".", "6:201", "-"),
        _row("x", "6:301-300", "GG"),  # insertion, start after end
    ]
    got = annotate.group_vep_rows(rows, alleles)
    assert [len(got[a.key]) for a in alleles] == [1, 1, 1]
    flipped = annotate.group_vep_rows([_row("x", "6:300-301", "GG")], alleles)
    assert list(flipped) == [alleles[2].key]
    assert annotate.group_vep_rows([_row("x", "6:999", "A")], alleles) == {}


def _tx(feature: str, **kwargs: str) -> dict:
    base = {"feature": feature, "feature_type": "Transcript", "mane_select": "",
            "canonical": ""}  # fmt: skip
    return base | kwargs


def test_choose_transcript_order() -> None:
    default = {"NM_2"}
    regulatory = {"feature": "ENSR1", "feature_type": "RegulatoryFeature",
                  "mane_select": "", "canonical": ""}  # fmt: skip
    rows = [
        regulatory,
        _tx("ENST1", canonical="YES"),
        _tx("NM_1.3", mane_select="NM_1.3"),
        _tx("NM_2.1"),
    ]
    assert annotate.choose_transcript(rows, default)[0]["feature"] == "NM_2.1"
    pick, source = annotate.choose_transcript(rows, set())
    assert (pick["feature"], source) == ("NM_1.3", "mane_select")
    pick, source = annotate.choose_transcript(
        [_tx("ENST1", canonical="YES"), _tx("NM_5.1", canonical="YES")], set()
    )
    assert (pick["feature"], source) == ("NM_5.1", "canonical_refseq")
    pick, source = annotate.choose_transcript(
        [_tx("ENST9"), _tx("ENST1", canonical="YES")], set()
    )
    assert (pick["feature"], source) == ("ENST1", "canonical")
    assert annotate.choose_transcript([_tx("ENST9")], set())[1] == (
        "first_transcript"
    )
    assert annotate.choose_transcript([regulatory], set())[1] == "first_row"
    assert annotate.choose_transcript([], set()) == (None, "none")


def test_annovar_genes_use_refgene_and_knowngene_only() -> None:
    annovar = {
        "Gene.refGene": "GENEA(dist=5),GENEB;GENEA",
        "Gene.ensGene": "ENSG00000000001;GENEE",
        "Gene.knownGene": "ENST00000000001.2;GENEK;GENEB",
    }
    assert annotate._annovar_genes(annovar) == ["GENEA", "GENEB", "GENEK"]
    assert annotate._annovar_genes({"Gene.refGene": "."}) == []
    assert annotate._annovar_genes({}) == []


def test_read_intervar_keeps_first_gene_line(tmp_path: Path) -> None:
    path = tmp_path / "x.intervar"
    path.write_text(
        "#Chr\tStart\tEnd\tRef\tAlt\tRef.Gene\t InterVar: InterVar and"
        " Evidence \n"
        "1\t5\t5\t-\tTT\tGENEA\t first\n"
        "1\t5\t5\t-\tTT\tGENEB\t second\n"
    )
    assert annotate.read_intervar(path) == {
        ("1", "5", "5", "-", "TT"): " first"
    }


def test_read_multianno_requires_otherinfo(tmp_path: Path) -> None:
    path = tmp_path / "m.txt"
    path.write_text("Chr\tStart\n6\t1\n")
    with pytest.raises(ValueError, match="Otherinfo"):
        annotate.read_multianno(path)


def test_is_pathogenic_is_exact() -> None:
    def record(sig: str = "", table: str = "", iv: str = "") -> dict:
        return {"annovar": {"ClinVar_SIG": sig}, "clinvar": table,
                "intervar_class": iv}  # fmt: skip

    assert annotate.is_pathogenic(record(sig="Likely_pathogenic"))
    assert annotate.is_pathogenic(record(sig="Pathogenic/Likely_pathogenic"))
    assert annotate.is_pathogenic(record(table="Likely pathogenic"))
    assert annotate.is_pathogenic(record(iv="Pathogenic"))
    assert annotate.is_pathogenic({**record(), "hgmd_class": "DM"})
    assert annotate.is_pathogenic({**record(), "hgmd_class": "DM?"})
    assert not annotate.is_pathogenic({**record(), "hgmd_class": "DP"})
    assert not annotate.is_pathogenic(
        record(sig="Conflicting_interpretations_of_pathogenicity")
    )
    conflict = "Conflicting classifications of pathogenicity[P(27)|LP(2)]"
    assert not annotate.is_pathogenic(
        record(sig="Pathogenic/Likely_pathogenic", table=conflict)
    )
    assert annotate.is_pathogenic(
        record(table=conflict, iv="Likely pathogenic")
    )
    assert not annotate.is_pathogenic(record(sig="Benign", iv="Likely benign"))
    assert not annotate.is_pathogenic(record(iv="Uncertain significance"))


def test_cache_is_not_reused_across_tool_versions(
    tmp_path, monkeypatch
) -> None:
    cache = annotate.RecordCache(tmp_path / "c.sqlite", "sum")
    cache.put_many([{"allele_key": "chr1:1:A:C"}])
    assert set(cache.get_many(["chr1:1:A:C"])) == {"chr1:1:A:C"}
    monkeypatch.setattr(annotate, "CACHE_SCHEMA", "2+99.0.0")
    assert cache.get_many(["chr1:1:A:C"]) == {}
    cache.close()
