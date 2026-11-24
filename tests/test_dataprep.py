# SPDX-License-Identifier: GPL-3.0-or-later
"""Tests for the annotation data builder, on tiny synthetic tables."""

import dataclasses
import gzip
import json
from pathlib import Path

import pysam
import pytest

from paralog2viewer import dataprep
from paralog2viewer.datalayout import DataLayout
from paralog2viewer.dataprep import (
    DataSources,
    RegionSet,
    build_data,
    build_regions,
    subset_annovar_table,
    subset_tabix,
    subset_table,
    uncovered_intervals,
    verify_data,
    write_annovar_idx,
)
from paralog2viewer.split import load_native_regions

LENGTHS = {"chr1": 5000, "chr2": 3000}


def _bgzip_index(text: str, dst: Path, **kwargs) -> None:
    plain = dst.with_suffix(".plain")
    plain.write_text(text)
    pysam.tabix_compress(str(plain), str(dst), force=True)
    pysam.tabix_index(str(dst), force=True, **kwargs)
    plain.unlink()


def _annovar_rows(chroms=("1", "2"), n=60, step=37) -> list[str]:
    return [
        f"{c}\t{100 + i * step}\t{100 + i * step}\tA\tG\t0.{i}\n"
        for c in chroms
        for i in range(n)
    ]


def _read_gz_lines(path: Path) -> list[str]:
    with gzip.open(path, "rt") as handle:
        return handle.read().splitlines()


# ---------------------------------------------------------------- regions


def test_regions_from_json_pad_clip_merge(tmp_path: Path) -> None:
    data = {
        # 1-based inclusive 1001-1100; native locus of GBA is added too.
        "GBA": {"phase_region": "38:chr1:1001-1100", "final_haplotypes": {}},
        "CFHclust": {"phase_region": "38:chr1:1500-1600,38:chr2:100-200"},
        "none": {"phase_region": None, "final_haplotypes": None},
    }
    js = tmp_path / "a.json"
    js.write_text(json.dumps(data))
    assert "GBA" in load_native_regions()
    regions = build_regions(js, 100, {"chr1": 300_000_000, "chr2": 3000})
    # 1000-1100 padded to 900-1200; 1499-1600 padded to 1399-1700;
    # chr2 99-200 padded to 0-300 (clipped at 0).
    assert ("chr1", 900, 1200) in regions
    assert ("chr1", 1399, 1700) in regions
    assert ("chr2", 0, 300) in regions
    # Native GBA locus (155210380-155219657, 1-based) is padded and merged.
    assert ("chr1", 155210279, 155219757) in regions
    assert regions == sorted(regions)


def test_regions_from_bed_merges_touching_and_clips(tmp_path: Path) -> None:
    bed = tmp_path / "r.bed"
    bed.write_text("# c\nchr1\t100\t200\nchr1\t300\t400\nchr2\t2900\t2990\n")
    regions = build_regions(bed, 50, LENGTHS)
    # 50-250 and 250-450 touch, so they merge; chr2 end is clipped to 3000.
    assert regions == [("chr1", 50, 450), ("chr2", 2850, 3000)]


def test_regions_reject_bad_input(tmp_path: Path) -> None:
    bed = tmp_path / "r.bed"
    bed.write_text("chrZ\t1\t10\n")
    with pytest.raises(ValueError, match="not in the reference"):
        build_regions(bed, 0, LENGTHS)
    bed.write_text("chr1\t1\t9999\n")
    with pytest.raises(ValueError, match="exceeds contig"):
        build_regions(bed, 0, LENGTHS)
    js = tmp_path / "a.json"
    js.write_text("{}")
    with pytest.raises(ValueError, match="mixed|either"):
        build_regions([js, bed], 0, LENGTHS)
    with pytest.raises(ValueError, match="no regions"):
        build_regions(js, 0, LENGTHS)
    js.write_text(json.dumps({"x": {"phase_region": "37:chr1:1-5"}}))
    with pytest.raises(ValueError, match="phase_region"):
        build_regions(js, 0, LENGTHS)


def test_uncovered_intervals_reports_offending_parts() -> None:
    regions = [("chr1", 100, 200), ("chr1", 300, 400)]
    calling = [
        ("chr1", 120, 180),  # inside
        ("chr1", 150, 350),  # gap 200-300
        ("chr1", 50, 120),  # starts before
        ("chr2", 0, 10),  # unknown contig
    ]
    assert uncovered_intervals(regions, calling) == [
        ("chr1", 200, 300),
        ("chr1", 50, 100),
        ("chr2", 0, 10),
    ]


# ------------------------------------------------------------ tabular data


def test_subset_table_strips_chr_for_matching_and_sorts(tmp_path: Path) -> None:
    src = tmp_path / "clinvar.tsv.gz"
    rows = ["Chromosome\tPositionVCF\tRef\n"]
    rows += [f"1\t{p}\tA\n" for p in (50, 150, 250, 350)]
    rows += [f"2\t{p}\tA\n" for p in (10, 20)]
    with gzip.open(src, "wt") as handle:
        handle.writelines(rows)
    dst = tmp_path / "out.tsv.gz"
    # Region names carry 'chr'; the table does not.
    kept = subset_table(src, dst, RegionSet([("chr1", 100, 300)]))
    assert kept == 2
    assert _read_gz_lines(dst) == [
        "Chromosome\tPositionVCF\tRef",
        "1\t150\tA",
        "1\t250\tA",
    ]
    with pysam.TabixFile(str(dst)) as tbx:
        assert [r.split("\t")[1] for r in tbx.fetch("1", 0, 1000)] == [
            "150",
            "250",
        ]


def test_subset_table_keeps_chr_prefix_and_boundaries(tmp_path: Path) -> None:
    src = tmp_path / "tb.tsv.gz"
    with gzip.open(src, "wt") as handle:
        handle.write("Chr\tPos\tRef\tFreq.\n")
        handle.writelines(f"chr1\t{p}\tA\tA:1\n" for p in (100, 101, 200, 201))
    dst = tmp_path / "out.tsv.gz"
    # BED 100-200 is 1-based 101..200: 100 and 201 fall outside.
    assert subset_table(src, dst, RegionSet([("chr1", 100, 200)])) == 2
    assert [x.split("\t")[1] for x in _read_gz_lines(dst)[1:]] == [
        "101",
        "200",
    ]


def test_subset_table_empty_result_is_valid(tmp_path: Path) -> None:
    src = tmp_path / "t.tsv.gz"
    with gzip.open(src, "wt") as handle:
        handle.write("Chr\tPos\n")
        handle.write("chr1\t5000\n")
    dst = tmp_path / "out.tsv.gz"
    assert subset_table(src, dst, RegionSet([("chr1", 0, 10)])) == 0
    assert _read_gz_lines(dst) == ["Chr\tPos"]
    assert Path(str(dst) + ".tbi").exists()


def test_subset_tabix_vcf_and_header(tmp_path: Path) -> None:
    text = (
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    text += "".join(
        f"{c}\t{p}\tid{p}\tA\tG\t.\t.\tX=1\n"
        for c, p in [("1", 90), ("1", 150), ("1", 500), ("2", 10)]
    )
    src = tmp_path / "hgmd.vcf.gz"
    _bgzip_index(text, src, preset="vcf")
    dst = tmp_path / "out.vcf.gz"
    kept = subset_tabix(src, dst, RegionSet([("chr1", 100, 200)]), preset="vcf")
    assert kept == 1
    assert _read_gz_lines(dst)[:2] == text.splitlines()[:2]
    with pysam.VariantFile(str(dst)) as vcf:
        assert [(r.contig, r.pos) for r in vcf] == [("1", 150)]


def test_subset_tabix_without_index_streams_same_rows(tmp_path: Path) -> None:
    text = (
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    text += "".join(
        f"{c}\t{p}\tid{p}\tA\tG\t.\t.\tX=1\n"
        for c, p in [("1", 90), ("1", 150), ("1", 200), ("1", 201), ("2", 10)]
    )
    regions = RegionSet([("chr1", 100, 200)])
    indexed = tmp_path / "indexed.vcf.gz"
    _bgzip_index(text, indexed, preset="vcf")
    plain = tmp_path / "noindex.vcf.gz"
    _bgzip_index(text, plain, preset="vcf")
    Path(str(plain) + ".tbi").unlink()
    out_a, out_b = tmp_path / "a.vcf.gz", tmp_path / "b.vcf.gz"
    assert subset_tabix(indexed, out_a, regions, preset="vcf") == 2
    assert subset_tabix(plain, out_b, regions, preset="vcf") == 2
    assert _read_gz_lines(out_a) == _read_gz_lines(out_b)
    assert _read_gz_lines(out_b)[:2] == text.splitlines()[:2]
    with pysam.VariantFile(str(out_b)) as vcf:
        assert [r.pos for r in vcf.fetch("1", 0, 1000)] == [150, 200]


def test_leftover_partial_files_do_not_break_rewrite(tmp_path: Path) -> None:
    src = tmp_path / "t.tsv.gz"
    with gzip.open(src, "wt") as handle:
        handle.write("Chr\tPos\nchr1\t150\n")
    dst = tmp_path / "out.tsv.gz"
    for leftover in (
        dst.name + ".plain.partial",
        dst.name + ".partial.gz",
        dst.name + ".partial.gz.tbi",
    ):
        (tmp_path / leftover).write_bytes(b"junk")
    dst.write_bytes(b"old")
    assert subset_table(src, dst, RegionSet([("chr1", 0, 200)])) == 1
    assert _read_gz_lines(dst) == ["Chr\tPos", "chr1\t150"]
    assert Path(str(dst) + ".tbi").exists()
    assert not [p for p in tmp_path.iterdir() if dataprep.PARTIAL in p.name]


def test_subset_tabix_hash_header_columns(tmp_path: Path) -> None:
    text = "# license\n#CHROM\tPOS\tREF\tALT\n"
    text += "chr1\t150\tA\tG\nchr1\t160\tA\tC\nchr2\t5\tA\tG\n"
    src = tmp_path / "am.tsv.gz"
    _bgzip_index(text, src, seq_col=0, start_col=1, end_col=1)
    dst = tmp_path / "out.tsv.gz"
    assert subset_tabix(src, dst, RegionSet([("chr1", 100, 155)])) == 1
    with pysam.TabixFile(str(dst)) as tbx:
        assert len(tbx.header) == 2
        assert list(tbx.fetch("chr1", 0, 1000)) == ["chr1\t150\tA\tG"]


# ----------------------------------------------------------------- ANNOVAR


def _check_idx(table: Path, allow_interleaved: bool = False) -> None:
    """Every idx range must contain the lines of its bin.

    With ``allow_interleaved`` a range may also hold rows of other bins (an
    unsorted table), but must start and end on a row of its own bin and
    cover all of them.
    """
    lines = open(Path(str(table) + ".idx")).read().splitlines()
    binsize, size = lines[0].split("\t")[1:]
    assert int(size) == table.stat().st_size
    raw = table.read_bytes()
    total = 0
    data = [x for x in raw.splitlines(keepends=True) if x[:1] != b"#"]
    for entry in lines[1:]:
        chrom, first, low, high = entry.split("\t")
        own = 0
        span = raw[int(low) : int(high)].splitlines()
        for line in span:
            fields = line.decode().split("\t")
            start = int(fields[1])
            mine = fields[0] == chrom and start - start % int(binsize) == int(
                first
            )
            assert mine or allow_interleaved
            own += mine
        assert own
        mine_all = [
            x
            for x in data
            if x.split(b"\t")[0].decode() == chrom
            and int(x.split(b"\t")[1]) // int(binsize) * int(binsize)
            == int(first)
        ]
        assert own == len(mine_all)
        total += own
    assert total == len(data)


def test_annovar_subset_with_index_matches_brute_force(tmp_path: Path) -> None:
    src = tmp_path / "hg38_exac03.txt"
    rows = _annovar_rows()
    src.write_text("#Chr\tStart\tEnd\tRef\tAlt\tAF\n" + "".join(rows))
    write_annovar_idx(src, 100)
    regions = RegionSet([("chr1", 500, 900), ("chr1", 1500, 1600)])
    dst = tmp_path / "out" / "hg38_exac03.txt"
    dst.parent.mkdir()
    kept = subset_annovar_table(src, dst, regions)
    expected = [
        r
        for r in rows
        if r.startswith("1\t")
        and (
            500 < int(r.split("\t")[1]) <= 900
            or 1500 < int(r.split("\t")[1]) <= 1600
        )
    ]
    assert kept == len(expected) > 0
    assert dst.read_text().splitlines() == (
        ["#Chr\tStart\tEnd\tRef\tAlt\tAF"] + [r.rstrip() for r in expected]
    )
    _check_idx(dst)
    # The new index is usable: ranges selected from it return the same rows.
    again = tmp_path / "out" / "again.txt"
    assert subset_annovar_table(dst, again, regions) == kept


def test_annovar_headerless_and_chr_prefixed_table(tmp_path: Path) -> None:
    src = tmp_path / "hg38_nci60.txt"
    rows = [
        r.replace("1\t", "chr1\t", 1) if r[0] == "1" else r
        for r in _annovar_rows()
    ]
    src.write_text("".join(rows))
    write_annovar_idx(src, 1000)
    dst = tmp_path / "o.txt"
    kept = subset_annovar_table(src, dst, RegionSet([("chr2", 0, 400)]))
    got = dst.read_text().splitlines()
    assert kept == len(got) > 0
    assert all(
        x.startswith("2\t") and int(x.split("\t")[1]) <= 400 for x in got
    )
    _check_idx(dst)


def test_annovar_stale_or_missing_index_falls_back_to_scan(
    tmp_path: Path,
) -> None:
    src = tmp_path / "t.txt"
    src.write_text("".join(_annovar_rows()))
    regions = RegionSet([("chr1", 500, 900)])
    plain = tmp_path / "plain.txt"
    n_scan = subset_annovar_table(src, plain, regions)
    assert not Path(str(plain) + ".idx").exists()
    write_annovar_idx(src, 100)
    indexed = tmp_path / "indexed.txt"
    assert subset_annovar_table(src, indexed, regions) == n_scan
    assert indexed.read_text() == plain.read_text()
    # Appending a row makes the index stale (size differs): scan, no new idx.
    with open(src, "a") as handle:
        handle.write("1\t600\t600\tA\tG\t0.1\n")
    stale = tmp_path / "stale.txt"
    assert subset_annovar_table(src, stale, regions) == n_scan + 1
    assert not Path(str(stale) + ".idx").exists()


def test_annovar_region_table_uses_overlap_of_zero_based_span(
    tmp_path: Path,
) -> None:
    src = tmp_path / "hg38_rmsk.txt"
    rows = [
        "585\t1\t1\t1\t1\tchr1\t0\t100\t-1\t+\tA\n",  # ends at 1-based 100
        "585\t1\t1\t1\t1\tchr1\t90\t150\t-1\t+\tB\n",  # overlaps
        "585\t1\t1\t1\t1\tchr1\t200\t250\t-1\t+\tC\n",  # starts at 201
        "585\t1\t1\t1\t1\tchr2\t100\t150\t-1\t+\tD\n",
    ]
    src.write_text("".join(rows))
    dst = tmp_path / "o.txt"
    kept = subset_annovar_table(
        src,
        dst,
        RegionSet([("chr1", 99, 200)]),  # 1-based 100..200
        chrom_col=5,
        start_col=6,
        end_col=7,
        zero_based_start=True,
    )
    assert kept == 2
    assert [x.split("\t")[10] for x in dst.read_text().splitlines()] == [
        "A",
        "B",
    ]


def test_annovar_unsorted_table_keeps_rows_and_warns(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    # Contigs interleaved and starts descending, as in hg38_regsnpintron.
    rows = [
        "1\t100\t100\tA\tG\t0\n",
        "1\t5000\t5000\tA\tG\t1\n",
        "2\t120\t120\tA\tG\t2\n",
        "1\t150\t150\tA\tG\t3\n",
        "1\t180\t180\tA\tG\t4\n",
        "2\t900\t900\tA\tG\t5\n",
        "1\t5100\t5100\tA\tG\t6\n",
    ]
    src = tmp_path / "hg38_unsorted.txt"
    src.write_text("#hdr\n" + "".join(rows))
    with caplog.at_level("WARNING"):
        write_annovar_idx(src, 1000)
    assert "hg38_unsorted.txt is not position-sorted" in caplog.text
    # ANNOVAR scans [low, high) of each bin and tests positions itself.
    _check_idx(src, allow_interleaved=True)
    regions = RegionSet([("chr1", 0, 200), ("chr1", 4900, 5200)])
    dst = tmp_path / "out.txt"
    caplog.clear()
    with caplog.at_level("WARNING"):
        kept = subset_annovar_table(src, dst, regions)
    assert "hg38_unsorted.txt is not position-sorted" in caplog.text
    # Same rows, same order as a plain scan (no index) of the source.
    plain_src = tmp_path / "plain" / "hg38_unsorted.txt"
    plain_src.parent.mkdir()
    plain_src.write_text(src.read_text())
    plain = tmp_path / "plain_out.txt"
    assert subset_annovar_table(plain_src, plain, regions) == kept == 5
    assert set(dst.read_text().splitlines()) == set(
        plain.read_text().splitlines()
    )
    assert dst.read_text() == plain.read_text()
    assert [x.split("\t")[5] for x in dst.read_text().splitlines()[1:]] == [
        "0",
        "1",
        "3",
        "4",
        "6",
    ]
    _check_idx(dst, allow_interleaved=True)


def test_select_vep_contigs_uses_synonyms(tmp_path: Path) -> None:
    cache = tmp_path / "116_GRCh38"
    for name in ("1", "KV575253.1"):
        (cache / name).mkdir(parents=True)
    (cache / "chr_synonyms.txt").write_text(
        "KV575253.1\tchr19_KV575253v1_alt\n"
    )
    found, missing = dataprep.select_vep_contigs(
        cache, ["chr1", "chr19_KV575253v1_alt", "chrUn_x"]
    )
    assert found == ["1", "KV575253.1"]
    assert missing == ["chrUn_x"]


def test_subset_table_keep_values_filters_rows(tmp_path: Path) -> None:
    src = tmp_path / "clinvar.tsv.gz"
    with gzip.open(src, "wt") as handle:
        handle.write("Chromosome\tPositionVCF\tSig\tAssembly\n")
        handle.write("1\t150\tBenign\tGRCh37\n")
        handle.write("1\t150\tPathogenic\tGRCh38\n")
        handle.write("1\t160\tBenign\tgrch38\n")
        handle.write("1\t170\tBenign\n")
    dst = tmp_path / "out.tsv.gz"
    kept = subset_table(
        src,
        dst,
        RegionSet([("chr1", 100, 300)]),
        keep_values={"Assembly": ["GRCh38"]},
    )
    assert kept == 2
    assert [x.split("\t")[2] for x in _read_gz_lines(dst)[1:]] == [
        "Pathogenic",
        "Benign",
    ]
    with pytest.raises(ValueError, match="no column"):
        subset_table(
            src,
            dst,
            RegionSet([("chr1", 100, 300)]),
            keep_values={"Build": ["GRCh38"]},
        )


def test_build_rejects_other_build_table_names(sources: DataSources) -> None:
    for bad in (
        {"hgmd_name": "HGMD_Pro_hg19.vcf.gz"},
        {"taiwan_biobank_name": "TB_population.tsv.gz"},
        {"clinvar_name": "ClinVar_GRCh37.tsv.gz"},
    ):
        with pytest.raises(ValueError, match="GRCh38"):
            dataprep._check_source_names(dataclasses.replace(sources, **bad))


# ---------------------------------------------------------------- end to end


@pytest.fixture
def sources(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> DataSources:
    root = tmp_path / "src"
    ref = root / "ref" / "g.fasta"
    ref.parent.mkdir(parents=True)
    ref.write_text(
        ">chr1\n" + "ACGT" * 1250 + "\n>chr2\n" + "ACGT" * 750 + "\n"
    )
    pysam.faidx(str(ref))
    (root / "ref" / "g.dict").write_text(
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:5000\n@SQ\tSN:chr2\tLN:3000\n"
    )
    vep = root / "vep"
    cache = vep / "homo_sapiens_merged" / "116_GRCh38"
    for contig in ("1", "2", "3"):
        (cache / contig).mkdir(parents=True)
        (cache / contig / "all_vars.gz").write_text(f"vars{contig}")
    (cache / "info.txt").write_text(
        "assembly\tGRCh38\nsource_ClinVar\t202509\nother\tx\n"
    )
    (cache / "chr_synonyms.txt").write_text("")
    (vep / "plugins").mkdir()
    (vep / "plugins" / "MaxEntScan.pm").write_text("package X;")
    (vep / "plugins" / "maxentscan").mkdir()
    (vep / "plugins" / "maxentscan" / "me2x5").write_text("model")
    am = vep / "plugins" / "alphamissense" / "AlphaMissense_hg38.tsv.gz"
    am.parent.mkdir()
    _bgzip_index(
        "#CHROM\tPOS\tREF\tALT\nchr1\t150\tA\tG\nchr1\t4000\tA\tG\n",
        am,
        seq_col=0,
        start_col=1,
        end_col=1,
    )
    humandb = root / "humandb"
    humandb.mkdir()
    monkeypatch.setattr(dataprep, "ANNOVAR_FILTER_DATABASES", ("exac03",))
    monkeypatch.setattr(
        dataprep, "ANNOVAR_REGION_DATABASES", (("rmsk", 5, 6, 7),)
    )
    monkeypatch.setattr(dataprep, "ANNOVAR_GENE_FILES", ("hg38_refGene.txt",))
    exac = humandb / "hg38_exac03.txt"
    exac.write_text("#Chr\tStart\n" + "".join(_annovar_rows()))
    write_annovar_idx(exac, 100)
    (humandb / "hg38_rmsk.txt").write_text(
        "585\t1\t1\t1\t1\tchr1\t120\t180\t-1\t+\tA\n"
        "585\t1\t1\t1\t1\tchr1\t4000\t4100\t-1\t+\tB\n"
    )
    (humandb / "hg38_refGene.txt").write_text("gene\n")
    intervar = root / "intervardb"
    intervar.mkdir()
    (intervar / "mim2gene.txt").write_text("x")
    (intervar / "PS4.variants.hg38").write_text("x")
    (intervar / "PS4.variants.hg19").write_text("x")
    profile = root / "profile"
    profile.mkdir()
    for name in dataprep._PROFILE_COPY:
        (profile / name).write_text(name)
    with gzip.open(profile / "ClinVar.tsv.gz", "wt") as handle:
        handle.write(
            "Chromosome\tPositionVCF\tAssembly\n"
            "1\t150\tGRCh38\n1\t150\tGRCh37\n1\t4000\tGRCh38\n"
        )
    # Other-build files sit in the sources but must not be copied.
    (profile / "HGMD_Pro_hg19.vcf.gz").write_text("hg19")
    (profile / "TB_population.tsv.gz").write_text("hg19")
    with gzip.open(profile / "TB_GRCh38.tsv.gz", "wt") as handle:
        handle.write("Chr\tPos\nchr1\t150\nchr1\t4000\n")
    _bgzip_index(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t150\tCM1\tA\tG\t.\t.\t.\n1\t4000\tCM2\tA\tG\t.\t.\t.\n",
        profile / "HGMD.vcf.gz",
        preset="vcf",
    )
    return DataSources(
        reference=ref,
        vep_cache=vep,
        annovar_humandb=humandb,
        intervar_db=intervar,
        profile=profile,
        taiwan_biobank_name="TB_GRCh38.tsv.gz",
        hgmd_name="HGMD.vcf.gz",
    )


def _build(tmp_path: Path, sources: DataSources) -> tuple[Path, dict]:
    bed = tmp_path / "regions.in.bed"
    bed.write_text("chr1\t100\t200\n")
    out = tmp_path / "data"
    manifest = build_data(
        out, sources, bed, padding=100, threads=2, expected_lengths=LENGTHS
    )
    return out, manifest


def test_build_and_verify_roundtrip(
    tmp_path: Path, sources: DataSources
) -> None:
    out, manifest = _build(tmp_path, sources)
    layout = DataLayout(out)
    assert manifest["genome_build"] == "GRCh38"
    assert manifest["regions"] == [["chr1", 0, 300]]
    assert layout.regions_bed.read_text() == "chr1\t0\t300\n"
    # Subsets keep only the in-region row (150) and drop 4000.
    assert _read_gz_lines(layout.clinvar)[1:] == ["1\t150\tGRCh38"]
    assert _read_gz_lines(layout.taiwan_biobank)[1:] == ["chr1\t150"]
    assert [row.split("\t")[1] for row in _read_gz_lines(layout.hgmd)[1:]] == [
        "150"
    ]
    assert [
        row.split("\t")[1] for row in _read_gz_lines(layout.alphamissense)[1:]
    ] == ["150"]
    exac = layout.annovar_humandb / "hg38_exac03.txt"
    _check_idx(exac)
    assert [x.split("\t")[1] for x in exac.read_text().splitlines()[1:]] == [
        str(100 + i * 37) for i in range(6) if 0 < 100 + i * 37 <= 300
    ]
    assert (layout.annovar_humandb / "hg38_rmsk.txt").read_text().count(
        "\n"
    ) == 1
    assert (layout.annovar_humandb / "hg38_refGene.txt").read_text() == "gene\n"
    # VEP: only contig 1 is copied; plugin code and MaxEntScan are present.
    merged = layout.vep_cache / "homo_sapiens_merged" / "116_GRCh38"
    assert (merged / "1" / "all_vars.gz").exists()
    assert not (merged / "2").exists()
    assert (layout.vep_plugins / "MaxEntScan.pm").exists()
    assert (layout.maxentscan / "me2x5").read_text() == "model"
    assert (layout.intervardb / "mim2gene.txt").exists()
    assert (layout.intervardb / "PS4.variants.hg38").exists()
    # No hg19/GRCh37 file of any source reaches the data directory.
    assert not [
        p.name
        for p in out.rglob("*")
        if "hg19" in p.name.lower() or p.name == "TB_population.tsv.gz"
    ]
    assert not (layout.profile / "HGMD_Pro_hg19.vcf.gz").exists()
    assert (layout.reference_fasta.with_suffix(".dict")).exists()
    assert Path(str(layout.reference_fasta) + ".fai").exists()
    assert manifest["vep_cache"] == {
        "assembly": "GRCh38",
        "source_ClinVar": "202509",
    }
    assert "manifest.json" not in manifest["files"]
    assert set(manifest["files"]) == {
        p.relative_to(out).as_posix()
        for p in out.rglob("*")
        if p.is_file() and p.name != "manifest.json"
    }
    assert verify_data(out) == []
    # A calling BED inside the regions passes.
    inside = tmp_path / "in.bed"
    inside.write_text("chr1\t10\t290\n")
    assert verify_data(out, inside) == []


def test_verify_reports_tampering_and_outside_calls(
    tmp_path: Path, sources: DataSources
) -> None:
    out, _ = _build(tmp_path, sources)
    layout = DataLayout(out)
    bed = tmp_path / "call.bed"
    bed.write_text("chr1\t250\t400\nchr2\t0\t5\n")
    problems = verify_data(out, bed)
    assert problems == [
        "outside regions: chr1:300-400",
        "outside regions: chr2:0-5",
    ]
    (layout.profile / "Comments.tsv").write_text("changed!")
    (layout.profile / "extra.txt").write_text("x")
    layout.maxentscan.joinpath("me2x5").unlink()
    problems = verify_data(out)
    assert "checksum differs: profile/Comments.tsv" in problems or (
        "size differs: profile/Comments.tsv" in problems
    )
    assert "missing file: maxentscan/me2x5" in problems
    assert "unlisted file: profile/extra.txt" in problems
    # Same-size edit is caught only by the checksum.
    (layout.profile / "Comments.tsv").write_text("Comments.tsv")
    (layout.profile / "Methylated.tsv").write_text("Methylated.xxx")
    quick = verify_data(out, checksums=False)
    assert not any("Methylated" in p for p in quick)
    assert "checksum differs: profile/Methylated.tsv" in verify_data(out)


def test_verify_rejects_missing_manifest_and_wrong_build(
    tmp_path: Path, sources: DataSources
) -> None:
    empty = tmp_path / "empty"
    empty.mkdir()
    assert verify_data(empty)[0].startswith("manifest missing")
    out, manifest = _build(tmp_path, sources)
    manifest["genome_build"] = "GRCh37"
    DataLayout(out).manifest.write_text(json.dumps(manifest))
    assert "genome build is GRCh37" in verify_data(out)


def test_build_rerun_over_partial_output(
    tmp_path: Path, sources: DataSources
) -> None:
    out, first = _build(tmp_path, sources)
    layout = DataLayout(out)
    layout.manifest.unlink()
    (out / dataprep.BUILD_SENTINEL).write_text("build in progress\n")
    layout.hgmd.write_bytes(b"truncated")
    (layout.hgmd.parent / "HGMD.partial.gz").write_bytes(b"junk")
    (layout.clinvar.parent / "ClinVar.tsv.gz.plain.partial").write_text("x")
    second = build_data(
        out,
        sources,
        tmp_path / "regions.in.bed",
        100,
        2,
        expected_lengths=LENGTHS,
    )
    assert second["data_version"] == first["data_version"]
    assert verify_data(out) == []
    assert not [p for p in out.rglob("*") if dataprep.PARTIAL in p.name]
    assert not (out / dataprep.BUILD_SENTINEL).exists()


def test_build_refuses_unrelated_non_empty_directory(
    tmp_path: Path, sources: DataSources
) -> None:
    bed = tmp_path / "regions.in.bed"
    bed.write_text("chr1\t100\t200\n")
    out = tmp_path / "unrelated"
    out.mkdir()
    keep = out / "notes.partial"
    keep.write_text("mine")
    with pytest.raises(ValueError, match="no manifest.json"):
        build_data(out, sources, bed, 100, 2, expected_lengths=LENGTHS)
    assert keep.read_text() == "mine"
    assert sorted(p.name for p in out.iterdir()) == ["notes.partial"]


def test_build_rejects_wrong_reference(
    tmp_path: Path, sources: DataSources
) -> None:
    bed = tmp_path / "regions.in.bed"
    bed.write_text("chr1\t100\t200\n")
    with pytest.raises(ValueError, match="not GRCh38"):
        build_data(tmp_path / "other", sources, bed, 100, 2)
