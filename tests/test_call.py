# SPDX-License-Identifier: GPL-3.0-or-later
"""Tests for per-haplotype calling and call combination."""

import logging
import random
import shutil
from pathlib import Path

import pysam
import pytest

from paralog2viewer.call import (
    GPU_FLAGS,
    clair3_command,
    combine_calls,
    filter_clair3_vcf,
    sawfish_commands,
)

HAP = "hapA"
NEEDS_BCFTOOLS = pytest.mark.skipif(
    shutil.which("bcftools") is None, reason="bcftools not on PATH"
)


def make_reference(path: Path) -> str:
    """Write a 600 bp chr1 with a C-AAAAAA-G homopolymer at 1-based 101."""
    rng = random.Random(1)
    seq = [rng.choice("ACGT") for _ in range(600)]
    seq[99:107] = list("CAAAAAAG")
    text = "".join(seq)
    path.write_text(f">chr1\n{text}\n")
    pysam.faidx(str(path))
    return text


def write_vcf(path: Path, header_extra: str, rows: list[str]) -> Path:
    """Write and bgzip a one-sample VCF named HAP."""
    lines = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=chr1,length=600>",
        '##FILTER=<ID=PASS,Description="p">',
        '##FILTER=<ID=LowQual,Description="l">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="g">',
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="g">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="d">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="a">',
        '##FORMAT=<ID=AF,Number=1,Type=Float,Description="a">',
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="t">',
        '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="l">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="e">',
        header_extra,
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{HAP}",
        *rows,
    ]
    path.write_text("\n".join(x for x in lines if x) + "\n")
    pysam.tabix_index(str(path), preset="vcf", force=True)
    return Path(f"{path}.gz")


def test_clair3_command_has_main_nf_options() -> None:
    cmd = clair3_command(
        Path("r.bam"), Path("ref.fa"), Path("r.bed"), "hapA",
        Path("out"), threads=4,
    )  # fmt: skip
    for option in (
        "--platform=ont", "--haploid_precise", "--snp_min_af=0.10",
        "--indel_min_af=0.13", "--var_pct_full=0.2", "--ref_pct_full=0.2",
        "--var_pct_phasing=0.8", "--no_phasing_for_fa",
        "--chunk_size=250000000", "--sample_name=hapA", "--threads=4",
        "--bed_fn=r.bed", "--bam_fn=r.bam", "--ref_fn=ref.fa",
    ):  # fmt: skip
        assert option in cmd


def test_gpu_flag_is_not_guessed(caplog: pytest.LogCaptureFixture) -> None:
    args = (Path("r.bam"), Path("f.fa"), Path("r.bed"), "h", Path("o"))
    base = clair3_command(*args, threads=2)
    with caplog.at_level(logging.WARNING, logger="paralog2viewer"):
        gpu = clair3_command(*args, threads=2, gpu=True)
    assert GPU_FLAGS == ()
    assert gpu == base
    assert "GPU flag is not verified" in caplog.text


def test_sawfish_commands_match_main_nf() -> None:
    discover, joint = sawfish_commands(
        Path("r.bam"), Path("f.fa"), Path("cn.bed"), Path("d"), Path("j"),
        threads=3,
    )  # fmt: skip
    assert discover[:2] == ["sawfish", "discover"]
    assert "--disable-cnv" in discover
    assert discover[discover.index("--cov-regex") + 1] == "."
    assert joint[:2] == ["sawfish", "joint-call"]
    assert "--report-supporting-reads" in joint
    assert "--treat-single-copy-as-haploid" in joint


def test_filter_clair3_keeps_only_records_in_bed(tmp_path: Path) -> None:
    make_reference(tmp_path / "ref.fa")
    bed = tmp_path / "region.bed"
    bed.write_text("chr1\t50\t450\n")
    src = write_vcf(
        tmp_path / "in.vcf",
        "",
        [
            "chr1\t60\t.\tA\tC\t10\tPASS\t.\tGT\t1",
            "chr1\t500\t.\tA\tC\t10\tPASS\t.\tGT\t1",
        ],
    )
    # Records need real REF only for the caller; filtering is by position.
    out = tmp_path / "out.vcf.gz"
    n_out = filter_clair3_vcf(src, out, bed, HAP, tmp_path / "x.tsv")
    assert n_out == 1
    with pysam.VariantFile(str(out)) as vcf:
        assert [r.pos for r in vcf] == [60]
    assert "chr1\t500" in (tmp_path / "x.tsv").read_text()
    assert Path(f"{out}.tbi").exists()
    with pytest.raises(ValueError, match="sample identity"):
        filter_clair3_vcf(src, out, bed, "other", tmp_path / "y.tsv")


@NEEDS_BCFTOOLS
def test_combine_calls(tmp_path: Path) -> None:
    ref = tmp_path / "ref.fa"
    seq = make_reference(ref)
    hap_dir = tmp_path / HAP
    hap_dir.mkdir()
    (hap_dir / "region.bed").write_text("chr1\t50\t450\n")
    snv_ref = seq[149]
    snv_alt = "A" if snv_ref != "A" else "C"
    out_ref = seq[459]
    clair = write_vcf(
        tmp_path / "clair.vcf",
        "",
        [
            # right-aligned 1-bp deletion in the homopolymer
            "chr1\t105\t.\tAA\tA\t30\tPASS\t.\tGT:GQ:DP:AD:AF\t1:40:20:2,18:0.9",
            f"chr1\t150\t.\t{snv_ref}\t{snv_alt}\t30\tPASS\t.\tGT:GQ:DP:AD:AF"
            "\t1:35:10:1,9:0.9",
            f"chr1\t160\t.\t{seq[159]}\tT\t30\tPASS\t.\tGT\t0",
            f"chr1\t170\t.\t{seq[169]}\tT\t30\tLowQual\t.\tGT\t1",
            f"chr1\t460\t.\t{out_ref}\tT\t30\tPASS\t.\tGT\t1",
        ],
    )
    sawfish = write_vcf(
        tmp_path / "saw.vcf",
        "",
        [
            f"chr1\t200\t.\t{seq[199]}\t<DEL>\t.\tPASS\t"
            "SVTYPE=DEL;SVLEN=-100;END=300\tGT:GQ:AD\t1:50:5,7",
            f"chr1\t250\t.\t{seq[249]}\t{seq[249]}[chr1:400[\t.\tPASS\t"
            "SVTYPE=BND\tGT:GQ:AD\t1:50:4,6",
            f"chr1\t260\t.\t{seq[259]}\t<DUP>\t.\tLowQual\t"
            "SVTYPE=DUP;SVLEN=10;END=270\tGT:GQ:AD\t1:5:4,6",
        ],
    )
    counts = combine_calls(hap_dir, ref, clair3_vcf=clair, sawfish_vcf=sawfish)
    assert counts["clair3:non_alt_or_nonhaploid"] == 1
    assert counts["clair3:filtered"] == 1
    assert counts["clair3:outside_bounds"] == 1
    assert counts["sawfish:filtered"] == 1
    assert counts["alleles_written"] == 4
    with pysam.VariantFile(str(hap_dir / "calls.vcf.gz")) as vcf:
        recs = list(vcf)
    by_alt = {r.alts[0]: r for r in recs}
    # Deletion left-aligned onto the C anchor before the homopolymer.
    deletion = by_alt["C"]
    assert (deletion.pos, deletion.ref) == (100, "CA")
    assert deletion.info["SRC"] == ("clair3",)
    sample = deletion.samples[HAP]
    assert sample["GT"] == (1,)
    assert (sample["DP"], sample["GQ"], sample["AD"]) == (20, 40, (2, 18))
    assert sample["AF"] == pytest.approx(0.9)
    # Symbolic and breakend alleles keep their true ALT and coordinates.
    sv = by_alt["<DEL>"]
    assert (sv.pos, sv.stop, sv.info["SVTYPE"]) == (200, 300, "DEL")
    assert sv.samples[HAP]["SR"] == 7
    assert sv.info["SRC"] == ("sawfish",)
    bnd = by_alt[f"{seq[249]}[chr1:400["]
    assert (bnd.pos, bnd.samples[HAP]["SR"]) == (250, 6)
    assert [r.pos for r in recs] == sorted(r.pos for r in recs)
    for rec in recs:
        assert seq[rec.start : rec.start + len(rec.ref)] == rec.ref
    assert Path(f"{hap_dir}/calls.vcf.gz.tbi").exists()


@NEEDS_BCFTOOLS
def test_combine_calls_rejects_ref_mismatch(tmp_path: Path) -> None:
    ref = tmp_path / "ref.fa"
    seq = make_reference(ref)
    hap_dir = tmp_path / HAP
    hap_dir.mkdir()
    (hap_dir / "region.bed").write_text("chr1\t50\t450\n")
    wrong = "A" if seq[149] != "A" else "C"
    clair = write_vcf(
        tmp_path / "c.vcf",
        "",
        [f"chr1\t150\t.\t{wrong}\tT\t30\tPASS\t.\tGT\t1"],
    )
    saw = write_vcf(tmp_path / "s.vcf", "", [])
    with pytest.raises(ValueError, match="REF differs"):
        combine_calls(hap_dir, ref, clair3_vcf=clair, sawfish_vcf=saw)
