# SPDX-License-Identifier: GPL-3.0-or-later
"""Runner tests with every heavy step replaced by a logging fake.

Workers are forked, so fakes log to a file instead of changing memory.
"""

import json
from pathlib import Path

import pysam
import pytest

from paralog2viewer import (
    annotate,
    call,
    cli,
    config,
    dataprep,
    merge,
    mviewer_table,
    plot,
    rccx,
    realign,
    runner,
    scheduler,
    split,
    vcf_out,
)

HAPS = {"H1": "geneA", "H2": "geneA", "H3": "geneB"}


def _touch(path: Path, text: str = "x") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


class Env:
    """Input files and the call log of one test."""

    def __init__(self, root: Path, *, rccx_region: bool = True):
        self.root = root
        self.log = root / "calls.log"
        self.bam = root / "in" / "S1.paraphase.bam"
        self.bam.parent.mkdir()
        header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
        with pysam.AlignmentFile(str(self.bam), "wb", header=header):
            pass
        pysam.index(str(self.bam))
        genes = {"geneA": {"sample_sex": "female"}}
        if rccx_region:
            genes["rccx"] = {"sample_sex": "female"}
        self.json = self.bam.with_suffix(".json")
        self.json.write_text(json.dumps(genes))
        self.ref = root / "ref.fa"
        self.ref.write_text(">chr1\n" + "A" * 1000 + "\n")
        pysam.faidx(str(self.ref))
        self.data = root / "data"
        self.data.mkdir()
        self.pseudo = _touch(root / "pseudo.tsv")
        self.out = root / "out"
        self.tmp = root / "scratch"
        self.tmp.mkdir()
        _touch(self.tmp / "other-project" / "keep.txt")

    def steps(self) -> list[str]:
        if not self.log.exists():
            return []
        return [
            line.split(" ", 1)[0]
            for line in self.log.read_text().split("\n")
            if line
        ]

    def lines(self, step: str) -> list[str]:
        return [
            line
            for line in self.log.read_text().split("\n")
            if line.startswith(step + " ")
        ]

    def settings(self, **kwargs) -> config.RunSettings:
        values = {
            "bam": self.bam,
            "reference": self.ref,
            "data": self.data,
            "outdir": self.out,
            "threads": 4,
            "memory_gb": 800.0,
            "pseudogene_tsv": self.pseudo,
        }
        values.update(kwargs)
        return config.RunSettings(**values)


@pytest.fixture
def env(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    return Env(tmp_path)


def install_fakes(monkeypatch, env: Env, *, fail_prepare: str | None = None):
    """Replaces every heavy step by a fake that logs its name."""

    def log(step, *info):
        with open(env.log, "a", encoding="utf-8") as handle:
            handle.write(" ".join([step, *map(str, info)]) + "\n")

    def fake_split(bam, js, ref, work, threads, *, sample=""):
        log("split", sample)
        jobs = []
        for hp, region in HAPS.items():
            src = _touch(work / "hap" / hp / "source_reads.bam")
            jobs.append(
                split.HaplotypeJob(region, hp, "chr1", 0, 100, 500, src, sample)
            )
        return jobs

    def fake_prepare(job, ref, work, threads):
        log("prepare", job.hp)
        if job.hp == fail_prepare:
            raise RuntimeError("minimap2 exploded")
        d = work / "hap" / job.hp
        return realign.PreparedHaplotype(
            job.region,
            job.hp,
            "chr1",
            _touch(d / "reads.bam"),
            _touch(d / "region.bed", "chr1\t0\t100\n"),
            _touch(d / "expected_cn.bed"),
            _touch(d / "alignment_audit.json"),
            0,
            100,
        )

    def fake_call(hap_dir, ref, *, threads=1, model_path=""):
        log("call", hap_dir.name, threads)
        return _touch(hap_dir / "calls.vcf.gz")

    def fake_merge(hap_dirs, ref, work):
        log("merge", len(hap_dirs))
        return _touch(work / "merged.vcf.gz")

    def fake_verify(data, calling_bed=None, *, checksums=True):
        if calling_bed is not None:
            log("verify_regions", len(calling_bed))
        return []

    def fake_annotate(merged, ref, layout, work, cache, threads, tools=None):
        log("annotate", cache.name)
        return _touch(work / "annotation" / "records.jsonl.gz")

    def fake_vcf(merged, records, out, mode="default"):
        log("vcf", out.name, mode)
        return _touch(out)

    def fake_table(records, merged, sample, sex, out, header_path=None):
        log("mviewer", out.name, sex)
        _touch(out)
        return 1

    def fake_plot(merged, hap_dirs, keys, pseudo, ref, out, **kwargs):
        log(
            "figures",
            ",".join(sorted(hap_dirs)),
            json.dumps(kwargs["region_of"], sort_keys=True),
        )
        return [_touch(out / "geneA.jpeg")]

    def fake_classify(js, merged, sites, **kwargs):
        log("rccx")
        return [{"hap": "H1"}]

    def fake_plot_rccx(rows, out_path, **kwargs):
        log("rccx_plot")
        return _touch(out_path)

    monkeypatch.setattr(dataprep, "verify_data", fake_verify)
    monkeypatch.setattr(split, "split_haplotypes", fake_split)
    monkeypatch.setattr(realign, "prepare_haplotype", fake_prepare)
    monkeypatch.setattr(call, "call_haplotype", fake_call)
    monkeypatch.setattr(merge, "merge_haplotype_calls", fake_merge)
    monkeypatch.setattr(annotate, "annotate_alleles", fake_annotate)
    monkeypatch.setattr(annotate, "pathogenic_allele_keys", lambda path: set())
    monkeypatch.setattr(vcf_out, "write_annotated_vcf", fake_vcf, raising=False)
    monkeypatch.setattr(
        mviewer_table, "write_mviewer_table", fake_table, raising=False
    )
    monkeypatch.setattr(plot, "plot_regions", fake_plot)
    monkeypatch.setattr(rccx, "classify_rccx", fake_classify)
    monkeypatch.setattr(
        rccx, "write_rccx_table", lambda rows, path: _touch(path)
    )
    monkeypatch.setattr(rccx, "plot_rccx_module", fake_plot_rccx)


def test_stale_figures_of_an_earlier_run_are_replaced(env, monkeypatch):
    install_fakes(monkeypatch, env)
    stale = env.out / "S1" / "figures" / "old_region.jpeg"
    stale.parent.mkdir(parents=True)
    stale.write_text("stale")
    result = runner.run_pipeline(env.settings(tmp=env.tmp))
    assert result.ok
    figures = env.out / "S1" / "figures"
    assert not stale.exists()
    assert (figures / "geneA.jpeg").exists()
    assert not (env.out / "S1" / "figures.tmp").exists()
    assert all(p.exists() for p in result.outputs)


def test_full_run_order_names_and_tmp_cleanup(env, monkeypatch):
    install_fakes(monkeypatch, env)
    result = runner.run_pipeline(env.settings(tmp=env.tmp))
    assert result.ok and result.work_removed
    steps = env.steps()
    assert steps[0] == "split"
    assert steps.index("prepare") < steps.index("call")
    assert steps.count("prepare") == steps.count("call") == len(HAPS)
    last_call = max(i for i, s in enumerate(steps) if s == "call")
    assert steps.index("merge") > last_call
    assert steps.index("verify_regions") > steps.index("merge")
    assert steps.index("annotate") > steps.index("verify_regions")
    for after in ("vcf", "mviewer", "figures"):
        assert steps.index(after) > steps.index("annotate")
    assert steps.index("rccx_plot") > steps.index("rccx") > steps.index("merge")
    # One thread per haplotype call; the per-run cache sits in the data dir.
    assert all(line.endswith("1") for line in env.lines("call"))
    assert env.lines("annotate") == ["annotate annotation_cache.sqlite"]
    sample_dir = env.out / "S1"
    assert (sample_dir / "S1.paralog.vcf.gz").exists()
    assert (sample_dir / "S1.paralog.tsv.gz").exists()
    assert (sample_dir / "figures" / "geneA.jpeg").exists()
    assert (sample_dir / "figures" / "RCCX_module.jpeg").exists()
    assert not (sample_dir / "timings.tsv").exists()
    assert "female" in env.lines("mviewer")[0]
    # region_of maps each HP to its split region.
    assert json.loads(env.lines("figures")[0].split(" ", 2)[2]) == HAPS
    # Only the run's own scratch subdirectory was removed.
    assert [p.name for p in env.tmp.iterdir()] == ["other-project"]
    assert (env.tmp / "other-project" / "keep.txt").exists()
    assert not any(p.is_symlink() for p in env.root.rglob("*"))


def test_format_vcf_only_still_annotates(env, monkeypatch):
    install_fakes(monkeypatch, env)
    runner.run_pipeline(
        env.settings(formats=("vcf",), figures=False, tmp=env.tmp)
    )
    steps = env.steps()
    assert "annotate" in steps and "vcf" in steps
    assert "mviewer" not in steps and "figures" not in steps
    assert "rccx" not in steps
    assert not (env.out / "S1" / "S1.paralog.tsv.gz").exists()


def test_format_mviewer_only_and_vcf_mode(env, monkeypatch):
    install_fakes(monkeypatch, env)
    runner.run_pipeline(
        env.settings(formats=("mviewer",), figures=False, vcf_annotation="full")
    )
    assert "vcf" not in env.steps()
    install_fakes(monkeypatch, env)
    runner.run_pipeline(
        env.settings(
            formats=("vcf",), figures=False, vcf_annotation="full", resume=False
        )
    )
    assert env.lines("vcf")[0].endswith("full")


def test_haplotype_failure_continues_and_exits_nonzero(
    env, monkeypatch, capsys
):
    install_fakes(monkeypatch, env, fail_prepare="H2")
    result = runner.run_pipeline(env.settings(tmp=env.tmp))
    assert not result.ok and not result.work_removed
    assert "haplotype H2" in result.failures
    assert "minimap2 exploded" in result.failures["haplotype H2"]
    assert result.work_dir.exists()
    assert env.lines("merge") == ["merge 2"]
    assert env.lines("figures")[0].split(" ")[1] == "H1,H3"
    assert (env.out / "S1" / "S1.paralog.vcf.gz").exists()
    rows = (result.work_dir / "timings.tsv").read_text()
    assert "prepare:H2\tH2\tfailed" in rows and "call:H2\tH2\tcancelled" in rows
    code = cli.main(
        [
            "run",
            "--bam",
            str(env.bam),
            "--reference",
            str(env.ref),
            "--data",
            str(env.data),
            "--outdir",
            str(env.out),
            "--pseudogene-tsv",
            str(env.pseudo),
            "--threads",
            "4",
            "--no-resume",
        ]
    )
    err = capsys.readouterr().err
    assert code == 1 and "H2" in err


def test_all_haplotypes_failing_is_an_error(env, monkeypatch):
    install_fakes(monkeypatch, env)
    monkeypatch.setattr(
        realign,
        "prepare_haplotype",
        lambda *a: (_ for _ in ()).throw(RuntimeError("bad")),
    )
    with pytest.raises(runner.PipelineError, match="no haplotype produced"):
        runner.run_pipeline(env.settings())
    assert Path("work/S1/timings.tsv").exists()


def test_global_step_failure_raises(env, monkeypatch):
    install_fakes(monkeypatch, env)
    monkeypatch.setattr(
        split,
        "split_haplotypes",
        lambda *a, **k: (_ for _ in ()).throw(ValueError("dict mismatch")),
    )
    with pytest.raises(runner.PipelineError, match="dict mismatch"):
        runner.run_pipeline(env.settings(tmp=env.tmp))
    # work directory kept for inspection after a failure
    assert len(list(env.tmp.glob("paralog2viewer.*"))) == 1


def test_dry_run_prints_plan_and_creates_nothing(env, monkeypatch, capsys):
    install_fakes(monkeypatch, env)
    result = runner.run_pipeline(env.settings(dry_run=True, tmp=env.tmp))
    text = capsys.readouterr().out
    assert result.dry_run and env.steps() == []
    for token in (
        "[1] split_haplotypes",
        "[2-4]",
        "[5]",
        "[6]",
        "S1.paralog.vcf.gz",
        "S1.paralog.tsv.gz",
        "RCCX_module.jpeg",
        "sex female",
    ):
        assert token in text
    assert not env.out.exists()
    assert [p.name for p in env.tmp.iterdir()] == ["other-project"]


def test_resume_skips_finished_steps(env, monkeypatch, caplog):
    install_fakes(monkeypatch, env)
    runner.run_pipeline(env.settings())
    first = env.steps()
    first_rows = Path("work/S1/timings.tsv").read_text().splitlines()
    first_split = next(r for r in first_rows if r.startswith("split\t"))
    assert "\tdone\t" in first_split
    env.log.unlink()
    with caplog.at_level("INFO"):
        result = runner.run_pipeline(env.settings())
    assert result.ok and env.steps() == []
    assert len(first) > 8
    rows = Path("work/S1/timings.tsv").read_text()
    resumed = next(r for r in rows.splitlines() if r.startswith("split\t"))
    assert "\tcached\t" in resumed and "skipped" not in rows
    # Stored values and finish time survive the resume.
    assert resumed.split("\t")[4:] == first_split.split("\t")[4:]
    assert "# total_wall_s\t" in rows and "# cpu_s_cached\t" in rows
    assert "run summary: wall" in caplog.text
    # A changed input BAM index reruns everything that depends on it.
    with open(f"{env.bam}.bai", "ab") as handle:
        handle.write(b"0")
    runner.run_pipeline(env.settings())
    assert env.steps()[0] == "split"


def test_work_directory_lock(env, monkeypatch):
    install_fakes(monkeypatch, env)
    with scheduler.work_lock(Path("work") / "S1"):
        with pytest.raises(runner.PipelineError, match="in use"):
            runner.run_pipeline(env.settings())
    assert env.steps() == []


def test_calling_bed_outside_data_regions_stops_run(env, monkeypatch):
    install_fakes(monkeypatch, env)
    monkeypatch.setattr(
        dataprep,
        "verify_data",
        lambda data, calling_bed=None, *, checksums=True: (
            ["outside regions: chr1:0-100"] if calling_bed else []
        ),
    )
    with pytest.raises(runner.PipelineError, match="outside"):
        runner.run_pipeline(env.settings())
    assert "annotate" not in env.steps()


@pytest.mark.parametrize(
    "mutate, message",
    [
        ("bai", "BAM index not found"),
        ("noref", "reference not found.*pass --reference"),
        ("sex", "inconsistent sample sex"),
        ("dict", "dictionary does not match"),
        ("name", "unsafe sample name"),
        ("data", "data directory"),
        ("suffix", "Paraphase .paraphase.bam"),
    ],
)
def test_input_checks(env, monkeypatch, mutate, message):
    install_fakes(monkeypatch, env)
    settings = env.settings()
    if mutate == "bai":
        Path(f"{env.bam}.bai").unlink()
    elif mutate == "noref":
        settings = env.settings(reference=env.data / "reference/genome.fa")
    elif mutate == "sex":
        env.json.write_text(
            json.dumps(
                {"a": {"sample_sex": "male"}, "b": {"sample_sex": "female"}}
            )
        )
    elif mutate == "dict":
        env.ref.write_text(">chr1\n" + "A" * 999 + "\n")
        pysam.faidx(str(env.ref))
    elif mutate == "name":
        settings = env.settings(sample="bad name")
    elif mutate == "data":
        monkeypatch.setattr(
            dataprep, "verify_data", lambda *a, **k: ["manifest missing: x"]
        )
    elif mutate == "suffix":
        plain = env.bam.with_name("S1.bam")
        plain.write_bytes(env.bam.read_bytes())
        settings = env.settings(bam=plain)
    with pytest.raises(runner.PipelineError, match=message):
        runner.check_inputs(settings)


def test_data_problems_ignore_cache_file(monkeypatch):
    monkeypatch.setattr(
        dataprep,
        "verify_data",
        lambda *a, **k: [
            "unlisted file: annotation_cache.sqlite-wal",
            "unlisted file: other.txt",
        ],
    )
    assert runner.data_problems(Path("d")) == ["unlisted file: other.txt"]


def test_sample_override_is_used_in_names(env, monkeypatch):
    install_fakes(monkeypatch, env)
    runner.run_pipeline(
        env.settings(sample="X_9", formats=("vcf",), figures=False)
    )
    assert (env.out / "X_9" / "X_9.paralog.vcf.gz").exists()


def test_work_dir_option_is_reused_and_kept(env, tmp_path, monkeypatch):
    install_fakes(monkeypatch, env)
    scratch = tmp_path / "scratch_wd"
    scratch.mkdir()
    work = tmp_path / "kept_work"
    result = runner.run_pipeline(
        env.settings(
            sample="X_9",
            tmp=scratch,
            work_dir=work,
            formats=("vcf",),
            figures=False,
        )
    )
    assert result.work_dir == work and work.exists()
    assert not result.work_removed
    assert list(scratch.iterdir()) == []


def test_remove_work_dir_guard(tmp_path):
    tmp = tmp_path / "scratch"
    tmp.mkdir()
    own = Path(tempfile_mkdtemp(tmp))
    other = tmp / "someone-else"
    other.mkdir()
    assert not runner.remove_work_dir(other, tmp, other)
    assert not runner.remove_work_dir(own, tmp, tmp / "paralog2viewer.zzz")
    assert not runner.remove_work_dir(own, tmp_path, own)
    assert other.exists() and own.exists()
    assert runner.remove_work_dir(own, tmp, own)
    assert not own.exists() and other.exists()


def tempfile_mkdtemp(parent: Path) -> str:
    import tempfile

    return tempfile.mkdtemp(prefix=runner.WORK_PREFIX, dir=str(parent))


def test_cache_falls_back_when_not_writable(tmp_path, monkeypatch, caplog):
    data = tmp_path / "data"
    data.mkdir()
    work = tmp_path / "work"
    assert runner.resolve_cache_path(None, data, work) == (
        data / runner.CACHE_NAME
    )
    custom = tmp_path / "mine.sqlite"
    assert runner.resolve_cache_path(custom, data, work) == custom
    monkeypatch.setattr("os.access", lambda *a, **k: False)
    with caplog.at_level("WARNING"):
        assert runner.resolve_cache_path(None, data, work) == (
            work / runner.CACHE_NAME
        )
    assert "not writable" in caplog.text


def test_publish_figures_updates_existing_dir_in_place(tmp_path):
    stage = tmp_path / "work" / "figures"
    stage.mkdir(parents=True)
    (stage / "geneA.jpeg").write_bytes(b"new")
    final = tmp_path / "out" / "figures"
    final.mkdir(parents=True)
    (final / "geneA.jpeg").write_bytes(b"old")
    (final / "stale.jpeg").write_bytes(b"old")
    leftover = tmp_path / "out" / "figures.tmp"
    leftover.mkdir()
    (leftover / "x.jpeg").write_bytes(b"x")
    assert runner._publish_figures(stage, final)
    assert sorted(p.name for p in final.iterdir()) == ["geneA.jpeg"]
    assert (final / "geneA.jpeg").read_bytes() == b"new"
    assert not leftover.exists()
