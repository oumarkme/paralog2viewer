# SPDX-License-Identifier: GPL-3.0-or-later
"""Tests for merging per-haplotype calls into one multi-sample VCF."""

import random
from pathlib import Path

import pysam
import pytest

from paralog2viewer.call import FORMAT_LINES, INFO_LINES
from paralog2viewer.merge import merge_haplotype_calls


def make_reference(path: Path) -> dict[str, str]:
    """Write chr1 (600 bp) and chr2 (300 bp) and return their sequences."""
    rng = random.Random(2)
    seqs = {
        "chr1": "".join(rng.choice("ACGT") for _ in range(600)),
        "chr2": "".join(rng.choice("ACGT") for _ in range(300)),
    }
    path.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
    pysam.faidx(str(path))
    return seqs


def make_hap(
    root: Path,
    name: str,
    chrom: str,
    region: tuple[int, int],
    reads: list[tuple[int, int]],
    rows: list[str],
) -> Path:
    """Create a haplotype directory with region.bed, reads.bam, calls."""
    hap = root / name
    hap.mkdir(parents=True)
    (hap / "region.bed").write_text(f"{chrom}\t{region[0]}\t{region[1]}\n")
    header = {"SQ": [{"SN": "chr1", "LN": 600}, {"SN": "chr2", "LN": 300}]}
    with pysam.AlignmentFile(str(hap / "reads.bam"), "wb", header=header) as f:
        for i, (start, length) in enumerate(reads):
            read = pysam.AlignedSegment(f.header)
            read.query_name = f"r{i}"
            read.query_sequence = "A" * length
            read.flag = 0
            read.reference_id = f.get_tid(chrom)
            read.reference_start = start
            read.mapping_quality = 60
            read.cigar = [(0, length)]
            f.write(read)
    pysam.index(str(hap / "reads.bam"))
    lines = [
        "##fileformat=VCFv4.2",
        f"##contig=<ID={chrom},length=600>",
        '##FILTER=<ID=PASS,Description="p">',
        *INFO_LINES,
        *FORMAT_LINES,
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{name}",
        *rows,
    ]
    (hap / "calls.vcf").write_text("\n".join(lines) + "\n")
    pysam.tabix_index(str(hap / "calls.vcf"), preset="vcf", force=True)
    return hap


@pytest.fixture
def merged(tmp_path: Path):
    """Two haplotypes on chr1; hapA calls an SNV and a deletion allele."""
    ref = tmp_path / "ref.fa"
    seqs = make_reference(ref)
    s = seqs["chr1"]
    alt = "A" if s[99] != "A" else "C"
    row_snv = (
        f"chr1\t100\t.\t{s[99]}\t{alt}\t.\tPASS\tSRC=clair3\t"
        "GT:DP:AD:AF:GQ\t1:20:2,18:0.9:40"
    )
    row_del = (
        f"chr1\t300\t.\t{s[299]}\t<DEL>\t.\tPASS\t"
        "SRC=sawfish;SVTYPE=DEL;SVLEN=-50;END=350\tGT:AD:SR\t1:3,5:5"
    )
    # hapA reads: 50-200 and 280-380. hapB reads: 50-150 only.
    a = make_hap(
        tmp_path / "w" / "hap", "hapA", "chr1", (50, 380),
        [(50, 150), (280, 100)], [row_snv, row_del],
    )  # fmt: skip
    b = make_hap(
        tmp_path / "w" / "hap", "hapB", "chr1", (50, 150),
        [(50, 100)], [],
    )  # fmt: skip
    out = merge_haplotype_calls([a, b], ref, tmp_path / "w")
    return out, seqs, ref


def test_genotypes_one_zero_missing(merged) -> None:
    out, _, _ = merged
    with pysam.VariantFile(str(out)) as vcf:
        assert list(vcf.header.samples) == ["hapA", "hapB"]
        snv, deletion = list(vcf)
    # hapA carries both; hapB covers the SNV site (0), not the DEL site (.).
    assert snv.samples["hapA"]["GT"] == (1,)
    assert snv.samples["hapB"]["GT"] == (0,)
    assert deletion.samples["hapA"]["GT"] == (1,)
    assert deletion.samples["hapB"]["GT"] == (None,)


def test_fields_and_info(merged) -> None:
    out, _, _ = merged
    with pysam.VariantFile(str(out)) as vcf:
        snv, deletion = list(vcf)
    a = snv.samples["hapA"]
    assert (a["DP"], a["AD"], a["GQ"]) == (20, (2, 18), 40)
    assert a["AF"] == pytest.approx(0.9)
    assert snv.samples["hapB"]["DP"] is None
    assert snv.info["SRC"] == ("clair3",)
    assert snv.info["REGION"] == ("chr1:51-380",)
    assert deletion.info["SRC"] == ("sawfish",)
    assert deletion.samples["hapA"]["SR"] == 5
    assert (deletion.alts, deletion.stop) == (("<DEL>",), 350)


def test_header_contigs_and_ref_match_fasta(merged) -> None:
    out, seqs, _ = merged
    with pysam.VariantFile(str(out)) as vcf:
        lengths = {c: vcf.header.contigs[c].length for c in vcf.header.contigs}
        assert lengths == {"chr1": 600, "chr2": 300}
        for rec in vcf:
            seq = seqs[rec.chrom]
            assert seq[rec.start : rec.start + len(rec.ref)] == rec.ref
    assert Path(f"{out}.tbi").exists()


def test_shared_allele_is_one_record_with_both_genotypes(
    tmp_path: Path,
) -> None:
    ref = tmp_path / "ref.fa"
    s = make_reference(ref)["chr1"]
    alt = "A" if s[99] != "A" else "C"
    row = f"chr1\t100\t.\t{s[99]}\t{alt}\t.\tPASS\tSRC=clair3\tGT:DP\t1:10"
    h1 = make_hap(tmp_path / "h", "hap1", "chr1", (50, 150), [(50, 100)], [row])
    h2 = make_hap(tmp_path / "h", "hap2", "chr1", (60, 160), [(60, 100)], [row])
    out = merge_haplotype_calls([h1, h2], ref, tmp_path / "w")
    with pysam.VariantFile(str(out)) as vcf:
        (rec,) = list(vcf)
    assert [rec.samples[h]["GT"] for h in ("hap1", "hap2")] == [(1,), (1,)]
    assert rec.info["REGION"] == ("chr1:51-150", "chr1:61-160")


def test_min_depth_and_empty_bam(tmp_path: Path) -> None:
    ref = tmp_path / "ref.fa"
    s = make_reference(ref)["chr1"]
    alt = "A" if s[99] != "A" else "C"
    row = f"chr1\t100\t.\t{s[99]}\t{alt}\t.\tPASS\tSRC=clair3\tGT\t1"
    h1 = make_hap(tmp_path / "h", "hap1", "chr1", (50, 150), [(50, 100)], [row])
    h2 = make_hap(tmp_path / "h", "hap2", "chr1", (50, 150), [(50, 100)], [])
    h3 = make_hap(tmp_path / "h", "hap3", "chr1", (50, 150), [], [])
    kw = {"reference": ref}
    out = merge_haplotype_calls([h1, h2, h3], work_dir=tmp_path / "w1", **kw)
    with pysam.VariantFile(str(out)) as vcf:
        (rec,) = list(vcf)
    assert rec.samples["hap2"]["GT"] == (0,)
    assert rec.samples["hap3"]["GT"] == (None,)
    out = merge_haplotype_calls(
        [h1, h2, h3], work_dir=tmp_path / "w2", min_depth=2, **kw
    )
    with pysam.VariantFile(str(out)) as vcf:
        (rec,) = list(vcf)
    assert rec.samples["hap2"]["GT"] == (None,)


def test_duplicate_hp_names_rejected(tmp_path: Path) -> None:
    ref = tmp_path / "ref.fa"
    make_reference(ref)
    a = make_hap(tmp_path / "x", "hap1", "chr1", (0, 100), [(0, 50)], [])
    b = make_hap(tmp_path / "y", "hap1", "chr1", (0, 100), [(0, 50)], [])
    with pytest.raises(ValueError, match="not unique"):
        merge_haplotype_calls([a, b], ref, tmp_path / "w")


def test_ref_mismatch_rejected(tmp_path: Path) -> None:
    ref = tmp_path / "ref.fa"
    s = make_reference(ref)["chr1"]
    wrong = "A" if s[99] != "A" else "C"
    row = f"chr1\t100\t.\t{wrong}\tT\t.\tPASS\tSRC=clair3\tGT\t1"
    h = make_hap(tmp_path / "h", "hap1", "chr1", (50, 150), [(50, 100)], [row])
    with pytest.raises(ValueError, match="REF differs"):
        merge_haplotype_calls([h], ref, tmp_path / "w")


def test_symbolic_alleles_keep_end_and_split_by_end(tmp_path: Path) -> None:
    """Same POS/ALT with different END stays two records, each with END."""
    ref = tmp_path / "ref.fa"
    s = make_reference(ref)["chr1"]
    row = (
        f"chr1\t300\t.\t{s[299]}\t<DUP>\t.\tPASS\t"
        "SRC=sawfish;SVTYPE=DUP;END={end}\tGT\t1"
    )
    h1 = make_hap(
        tmp_path / "h", "hap1", "chr1", (50, 450), [(280, 100)],
        [row.format(end=350)],
    )  # fmt: skip
    h2 = make_hap(
        tmp_path / "h", "hap2", "chr1", (50, 450), [(280, 100)],
        [row.format(end=400)],
    )  # fmt: skip
    out = merge_haplotype_calls([h1, h2], ref, tmp_path / "w")
    with pysam.VariantFile(str(out)) as vcf:
        recs = list(vcf)
    assert [r.stop for r in recs] == [350, 400]
    assert [str(r).split("\t")[7].count("END=") for r in recs] == [1, 1]
    assert "END=350" in str(recs[0]) and "END=400" in str(recs[1])
    assert [r.samples["hap1"]["GT"] for r in recs] == [(1,), (0,)]
    assert [r.samples["hap2"]["GT"] for r in recs] == [(0,), (1,)]
