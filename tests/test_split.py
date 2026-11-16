# SPDX-License-Identifier: GPL-3.0-or-later
import json
import random
from pathlib import Path

import pysam
import pytest

from paralog2viewer.split import (
    load_native_regions,
    parse_phase_region,
    split_haplotypes,
)

BLUE = "166,206,227"
CHROMS = {"chrT": 2000, "chrN": 3000}


def make_reference(path: Path) -> Path:
    rng = random.Random(1)
    with path.open("w") as out:
        for name, length in CHROMS.items():
            seq = "".join(rng.choice("ACGT") for _ in range(length))
            out.write(f">{name}\n{seq}\n")
    pysam.faidx(str(path))
    return path


def make_bam(path: Path, records: list[dict], lengths=CHROMS) -> Path:
    """Write a sorted, indexed BAM. Each record: name, hp, rn, yc, ..."""
    header = {"SQ": [{"SN": n, "LN": ln} for n, ln in lengths.items()]}
    unsorted = path.with_suffix(".unsorted.bam")
    with pysam.AlignmentFile(str(unsorted), "wb", header=header) as out:
        for rec in records:
            read = pysam.AlignedSegment(out.header)
            read.query_name = rec["name"]
            read.flag = rec.get("flag", 0)
            read.reference_id = list(lengths).index(rec.get("chrom", "chrT"))
            read.reference_start = rec.get("pos", 100)
            read.mapping_quality = 60
            seq = "ACGT" * (rec.get("length", 50) // 4)
            read.query_sequence = seq
            read.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
            read.cigartuples = [(0, len(seq))]
            read.set_tag("RN", rec["rn"])
            read.set_tag("HP", rec["hp"])
            read.set_tag("YC", rec.get("yc", BLUE))
            out.write(read)
    pysam.sort("-o", str(path), str(unsorted))
    unsorted.unlink()
    pysam.index(str(path))
    return path


def make_json(path: Path, genes: dict) -> Path:
    path.write_text(json.dumps(genes))
    return path


def gene(haps: list[str], region="38:chrT:101-1500", **extra) -> dict:
    return {
        "final_haplotypes": {f"h{i}": h for i, h in enumerate(haps)},
        "phase_region": region,
        "haplotype_details": {h: {"boundary": [1, 2]} for h in haps},
        **extra,
    }


def names(path: Path) -> list[str]:
    with pysam.AlignmentFile(str(path), "rb") as bam:
        return sorted(r.query_name for r in bam.fetch(until_eof=True))


@pytest.fixture
def ref(tmp_path):
    return make_reference(tmp_path / "ref.fa")


def test_selects_confident_primary_reads_by_hap(tmp_path, ref):
    records = [
        {"name": "r1", "rn": "g", "hp": "g_hap1", "length": 60},
        {"name": "r1", "rn": "g", "hp": "g_hap1", "flag": 2048, "pos": 900},
        {"name": "r2", "rn": "g", "hp": "g_hap1", "pos": 200},
        {"name": "r3", "rn": "g", "hp": "g_hap2"},
        {"name": "r4", "rn": "g", "hp": "g_hap1", "yc": "0,0,0"},
        {"name": "r5", "rn": "g", "hp": "g_hap1", "flag": 256},
    ]
    bam = make_bam(tmp_path / "p.bam", records)
    js = make_json(
        tmp_path / "p.json",
        {"g": gene(["g_hap1", "g_hap2"]), "empty": {"final_haplotypes": {}}},
    )
    work = tmp_path / "W"
    jobs = split_haplotypes(bam, js, ref, work, 2, sample="S1")
    by_hp = {j.hp: j for j in jobs}
    assert set(by_hp) == {"g_hap1", "g_hap2"}
    job = by_hp["g_hap1"]
    assert job.region == "g"
    assert (job.chrom, job.core_start, job.core_end) == ("chrT", 100, 1500)
    assert job.max_read_length == 60
    assert job.source_bam == work / "hap" / "g_hap1" / "source_reads.bam"
    # r4 (unconfident) and r5 (secondary) excluded; r1 supplementary kept.
    assert names(job.source_bam) == ["r1", "r1", "r2"]
    assert names(by_hp["g_hap2"].source_bam) == ["r3"]
    with pysam.AlignmentFile(str(job.source_bam), "rb") as out:
        assert out.header.to_dict()["RG"][0]["DS"] == "S1"
        tags = {
            (r.get_tag("HP"), r.get_tag("RN"), r.get_tag("RG")) for r in out
        }
    assert tags == {("g_hap1", "g", "g_hap1")}
    audit = json.loads((work / "split_audit.json").read_text())
    assert [g["status"] for g in audit["genes"]] == [
        "no_resolved_haplotypes",
        "resolved",
    ]


def test_chimeric_duplicate_primary_kept_once(tmp_path, ref):
    records = [
        {"name": "r1", "rn": "g", "hp": "g_hap1", "pos": 100},
        {"name": "r1", "rn": "g", "hp": "g_hap1", "pos": 700},
    ]
    bam = make_bam(tmp_path / "p.bam", records)
    js = make_json(tmp_path / "p.json", {"g": gene(["g_hap1"])})
    (job,) = split_haplotypes(bam, js, ref, tmp_path / "W", 1)
    assert names(job.source_bam) == ["r1"]


def test_haplotype_without_reads_is_not_returned(tmp_path, ref):
    bam = make_bam(
        tmp_path / "p.bam", [{"name": "r1", "rn": "g", "hp": "g_hap1"}]
    )
    js = make_json(tmp_path / "p.json", {"g": gene(["g_hap1", "g_hap2"])})
    jobs = split_haplotypes(bam, js, ref, tmp_path / "W", 1)
    assert [j.hp for j in jobs] == ["g_hap1"]


def test_dictionary_mismatch_raises(tmp_path, ref):
    bam = make_bam(
        tmp_path / "p.bam",
        [{"name": "r1", "rn": "g", "hp": "g_hap1"}],
        lengths={"chrT": 1999, "chrN": 3000},
    )
    js = make_json(tmp_path / "p.json", {"g": gene(["g_hap1"])})
    with pytest.raises(ValueError, match="sequence dictionary"):
        split_haplotypes(bam, js, ref, tmp_path / "W", 1)


def test_hp_reused_across_regions_raises(tmp_path, ref):
    bam = make_bam(
        tmp_path / "p.bam",
        [
            {"name": "r1", "rn": "a", "hp": "same"},
            {"name": "r2", "rn": "b", "hp": "same"},
        ],
    )
    js = make_json(
        tmp_path / "p.json", {"a": gene(["same"]), "b": gene(["same"])}
    )
    with pytest.raises(ValueError, match="must be unique"):
        split_haplotypes(bam, js, ref, tmp_path / "W", 1)


def test_conflicting_hp_for_one_read_raises(tmp_path, ref):
    bam = make_bam(
        tmp_path / "p.bam",
        [
            {"name": "r1", "rn": "g", "hp": "g_hap1"},
            {"name": "r1", "rn": "g", "hp": "g_hap2", "pos": 700},
        ],
    )
    js = make_json(tmp_path / "p.json", {"g": gene(["g_hap1", "g_hap2"])})
    with pytest.raises(ValueError, match="Conflicting"):
        split_haplotypes(bam, js, ref, tmp_path / "W", 1)


def test_unsafe_label_raises(tmp_path, ref):
    bam = make_bam(tmp_path / "p.bam", [{"name": "r", "rn": "g", "hp": "a b"}])
    js = make_json(tmp_path / "p.json", {"g": gene(["a b"])})
    with pytest.raises(ValueError, match="Unsafe"):
        split_haplotypes(bam, js, ref, tmp_path / "W", 1)


def test_gene2_uses_native_only_with_boundary_and_reads(tmp_path, ref):
    native = {"g": ("chrN", 501, 1500)}
    records = [
        {
            "name": "n1",
            "rn": "g",
            "hp": "g_gene2hap1",
            "chrom": "chrN",
            "pos": 600,
        },
        {
            "name": "n2",
            "rn": "g",
            "hp": "g_gene2hap2",
            "chrom": "chrN",
            "pos": 600,
        },
        {"name": "n3", "rn": "g", "hp": "g_gene2hap3", "pos": 300},
        {"name": "n4", "rn": "g", "hp": "g_gene2hap4", "pos": 300},
    ]
    bam = make_bam(tmp_path / "p.bam", records)
    haps = [f"g_gene2hap{i}" for i in range(1, 5)]
    result = gene(haps)
    details = result["haplotype_details"]
    details[haps[0]]["boundary_gene2"] = [10, 20]  # complete, reads native
    details[haps[1]]["boundary_gene2"] = [10, None]  # incomplete
    details[haps[2]]["boundary_gene2"] = [10, 20]  # complete, no native reads
    # haps[3]: no boundary_gene2 at all
    js = make_json(tmp_path / "p.json", {"g": result})
    jobs = split_haplotypes(
        bam, js, ref, tmp_path / "W", 1, native_regions=native
    )
    by_hp = {j.hp: (j.chrom, j.core_start, j.core_end) for j in jobs}
    assert by_hp[haps[0]] == ("chrN", 500, 1500)
    for hp in haps[1:]:
        assert by_hp[hp] == ("chrT", 100, 1500)
    audit = json.loads((tmp_path / "W" / "split_audit.json").read_text())
    status = {
        h["hap"]: h["native_projection_status"] for h in audit["haplotypes"]
    }
    assert status == {
        haps[0]: "used",
        haps[1]: "missing_native_boundary",
        haps[2]: "no_native_primary_reads",
        haps[3]: "missing_native_boundary",
    }


def test_packaged_native_regions_table():
    table = load_native_regions()
    assert table["smn1"] == ("chr5", 70040526, 70088546)
    assert len(table) == 10


def test_parse_phase_region_accepts_comma_separated_loci() -> None:
    assert parse_phase_region("38:chrT:101-500, 38:chrT:900-1500") == [
        ("chrT", 101, 500),
        ("chrT", 900, 1500),
    ]
    for bad in ("37:chrT:1-5", "38:chrT:1-5,", "chrT:1-5"):
        with pytest.raises(ValueError, match="phase_region"):
            parse_phase_region(bad)


def test_comma_separated_phase_region_uses_covering_span(tmp_path, ref):
    bam = make_bam(
        tmp_path / "p.bam", [{"name": "r1", "rn": "g", "hp": "g_hap1"}]
    )
    region = "38:chrT:101-500,38:chrT:900-1500"
    js = make_json(tmp_path / "p.json", {"g": gene(["g_hap1"], region)})
    (job,) = split_haplotypes(bam, js, ref, tmp_path / "W", 1)
    assert (job.chrom, job.core_start, job.core_end) == ("chrT", 100, 1500)


def test_phase_region_on_two_contigs_is_rejected(tmp_path, ref):
    bam = make_bam(
        tmp_path / "p.bam", [{"name": "r1", "rn": "g", "hp": "g_hap1"}]
    )
    region = "38:chrT:101-500,38:chrN:900-1500"
    js = make_json(tmp_path / "p.json", {"g": gene(["g_hap1"], region)})
    with pytest.raises(ValueError, match="several contigs"):
        split_haplotypes(bam, js, ref, tmp_path / "W", 1)
