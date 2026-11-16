# SPDX-License-Identifier: GPL-3.0-or-later
from pathlib import Path

import pytest

from paralog2viewer import cli, config, dataprep, rccx, runner


def test_version_prints_0_9_0(capsys):
    with pytest.raises(SystemExit) as exc:
        cli.main(["--version"])
    assert exc.value.code == 0
    assert capsys.readouterr().out.strip() == "paralog2viewer 0.9.0"


def test_run_defaults_both_formats_and_figures():
    args = cli.build_parser().parse_args(
        ["run", "--bam", "a.bam", "--data", "d"]
    )
    assert args.formats == ("vcf", "mviewer")
    assert args.figures is True
    assert args.threads is None and args.memory is None
    assert args.resume is True and args.dry_run is False
    assert args.vcf_annotation == "default"


def test_run_format_and_no_figures():
    args = cli.build_parser().parse_args(
        ["run", "--bam", "a.bam", "--format", "vcf", "--no-figures"]
    )
    assert args.formats == ("vcf",)
    assert args.figures is False


def test_invalid_format_errors(capsys):
    with pytest.raises(SystemExit) as exc:
        cli.main(["run", "--bam", "a.bam", "--data", "d", "--format", "bam"])
    assert exc.value.code == 2
    assert "unknown --format" in capsys.readouterr().err


def test_missing_data_errors(monkeypatch, capsys):
    monkeypatch.delenv(config.DATA_ENV_VAR, raising=False)
    with pytest.raises(SystemExit) as exc:
        cli.main(["run", "--bam", "a.bam"])
    assert exc.value.code == 2
    assert "--data is required" in capsys.readouterr().err


def test_run_options_reach_settings(monkeypatch):
    seen = {}

    def fake_run(settings):
        seen["s"] = settings
        return runner.RunResult("S", Path("o"), dry_run=True)

    monkeypatch.setenv(config.DATA_ENV_VAR, "/env")
    monkeypatch.setattr(runner, "run_pipeline", fake_run)
    code = cli.main(
        [
            "run",
            "--bam",
            "a.paraphase.bam",
            "--threads",
            "1",
            "--memory",
            "100",
            "--cache",
            "c.db",
            "--pseudogene-tsv",
            "p.tsv",
            "--vcf-annotation",
            "full",
            "--dry-run",
            "--no-resume",
            "--tmp",
            "t",
        ]
    )
    s = seen["s"]
    assert code == 0 and str(s.data) == "/env"
    assert s.threads == 1
    assert s.memory_gb == 100 and str(s.cache) == "c.db"
    assert s.vcf_annotation == "full" and s.dry_run and not s.resume


def test_threads_above_affinity_are_clamped_with_warning(monkeypatch, caplog):
    seen = {}

    def fake_run(settings):
        seen["s"] = settings
        return runner.RunResult("S", Path("o"), dry_run=True)

    monkeypatch.setenv(config.DATA_ENV_VAR, "/env")
    monkeypatch.setattr(runner, "run_pipeline", fake_run)
    usable = config.default_threads()
    with caplog.at_level("WARNING", logger="paralog2viewer"):
        cli.main(["run", "--bam", "a.bam", "--threads", str(usable + 5)])
    assert seen["s"].threads == usable
    assert "exceeds" in caplog.text


def test_bad_threads_and_memory(capsys):
    for extra in (["--threads", "0"], ["--memory", "0"], ["--cpus", "0-1"]):
        with pytest.raises(SystemExit) as exc:
            cli.main(["run", "--bam", "a.bam", "--data", "d", *extra])
        assert exc.value.code == 2


def test_pipeline_error_returns_1(monkeypatch, capsys):
    def boom(settings):
        raise runner.PipelineError("bad input")

    monkeypatch.setattr(runner, "run_pipeline", boom)
    assert cli.main(["run", "--bam", "a.bam", "--data", "d"]) == 1
    assert "error: bad input" in capsys.readouterr().err


def test_build_data_runs_build_then_rccx_and_updates_manifest(
    tmp_path, monkeypatch, capsys
):
    out = tmp_path / "data"
    calls = []

    def fake_build(out_dir, sources, regions, padding, threads):
        calls.append(("build", sources.reference.name, list(regions), padding))
        calls.append(("sources", sources))
        (out_dir / "ref").mkdir(parents=True)
        (out_dir / "ref" / "g.fa").write_text(">c\nA\n")
        (out_dir / "manifest.json").write_text('{"files": {}}')
        return {}

    def fake_sites(reference, genes, out_tsv, *, threads=1):
        calls.append(("sites", reference.name))
        out_tsv.write_text("family\n")
        return 7

    monkeypatch.setattr(dataprep, "build_data", fake_build)
    monkeypatch.setattr(rccx, "build_rccx_sites", fake_sites)
    code = cli.main(
        [
            "build-data",
            "--out",
            str(out),
            "--regions",
            "a.json",
            "b.json",
            "--reference",
            "/x/ref.fa",
            "--vep-cache",
            "/x/vep",
            "--annovar-db",
            "/x/humandb",
            "--intervar-db",
            "/x/intervardb",
            "--mviewer-profile",
            "/x/profile",
            "--padding",
            "10",
        ]
    )
    assert code == 0 and "7 RCCX sites" in capsys.readouterr().out
    assert calls[0] == ("build", "ref.fa", [Path("a.json"), Path("b.json")], 10)
    sources = calls[1][1]
    assert sources.vep_cache == Path("/x/vep")
    assert sources.annovar_humandb == Path("/x/humandb")
    assert sources.intervar_db == Path("/x/intervardb")
    assert sources.profile == Path("/x/profile")
    assert calls[2][0] == "sites"
    import json

    files = json.loads((out / "manifest.json").read_text())["files"]
    assert files["rccx/discriminating_sites.tsv"]["size"] == 7


def test_build_data_requires_regions():
    with pytest.raises(SystemExit) as exc:
        cli.main(["build-data", "--out", "o"])
    assert exc.value.code == 2


def test_build_data_source_directories_are_required(capsys):
    base = ["build-data", "--out", "o", "--regions", "a.bed"]
    sources = {
        "--reference": "r.fa",
        "--vep-cache": "v",
        "--annovar-db": "a",
        "--intervar-db": "i",
        "--mviewer-profile": "m",
    }
    for omitted in sources:
        argv = list(base)
        for flag, value in sources.items():
            if flag != omitted:
                argv += [flag, value]
        with pytest.raises(SystemExit) as exc:
            cli.main(argv)
        assert exc.value.code == 2
        assert omitted in capsys.readouterr().err


def test_run_reference_defaults_to_data_directory(monkeypatch):
    seen = {}

    def fake_run(settings):
        seen["s"] = settings
        return runner.RunResult("S", Path("o"), dry_run=True)

    monkeypatch.setattr(runner, "run_pipeline", fake_run)
    argv = ["run", "--bam", "a.bam", "--data", "/d", "--dry-run"]
    assert cli.main(argv) == 0
    assert seen["s"].reference == Path("/d/reference/genome.fa")
    assert cli.main([*argv, "--reference", "/r/g.fa"]) == 0
    assert seen["s"].reference == Path("/r/g.fa")


def test_run_pseudogene_table_defaults_to_package():
    settings = config.RunSettings(
        bam=Path("a.bam"),
        reference=Path("r.fa"),
        data=Path("d"),
        outdir=Path("o"),
    )
    path = runner.pseudogene_path(settings)
    assert path.parent.name == "resources" and path.is_file()
    assert path.name == "known_gene_pseudogene_variants.tsv"
    custom = Path("mine.tsv")
    assert (
        runner.pseudogene_path(
            config.RunSettings(
                bam=Path("a.bam"),
                reference=Path("r.fa"),
                data=Path("d"),
                outdir=Path("o"),
                pseudogene_tsv=custom,
            )
        )
        == custom
    )


def test_check_reports_missing_tools_and_data(monkeypatch, capsys, tmp_path):
    present = {"samtools": "/bin/samtools"}
    monkeypatch.setattr(cli.shutil, "which", lambda name: present.get(name))
    monkeypatch.setattr(
        runner, "data_problems", lambda d, checksums: ["missing file: x"]
    )
    monkeypatch.delenv(config.DATA_ENV_VAR, raising=False)
    code = cli.main(["check", "--data", str(tmp_path)])
    out = capsys.readouterr().out
    assert code == 1
    assert "samtools" in out and "/bin/samtools" in out
    assert "minimap2      MISSING" in out
    assert "PROBLEMS" in out and "missing file: x" in out


def test_check_ok_when_everything_present(monkeypatch, capsys):
    monkeypatch.setattr(cli.shutil, "which", lambda name: f"/bin/{name}")
    monkeypatch.setattr(cli.os.path, "isfile", lambda path: True)
    monkeypatch.delenv(config.DATA_ENV_VAR, raising=False)
    assert cli.main(["check"]) == 0
    assert "MISSING" not in capsys.readouterr().out


def test_check_reports_missing_annovar_intervar_python(monkeypatch, capsys):
    monkeypatch.setattr(
        cli.shutil, "which", lambda name: None if name == "python3" else "/b"
    )
    monkeypatch.setattr(cli.os.path, "isfile", lambda path: False)
    monkeypatch.delenv(config.DATA_ENV_VAR, raising=False)
    assert cli.main(["check"]) == 1
    out = capsys.readouterr().out
    assert "python3       MISSING" in out
    assert "/opt/annovar/table_annovar.pl MISSING" in out
    assert "/opt/intervar/Intervar.py MISSING" in out


def test_check_quick_skips_checksums(monkeypatch, tmp_path):
    monkeypatch.setattr(cli.shutil, "which", lambda name: f"/bin/{name}")
    monkeypatch.setattr(cli.os.path, "isfile", lambda path: True)
    seen = {}
    monkeypatch.setattr(
        runner,
        "data_problems",
        lambda d, checksums: seen.update(c=checksums) or [],
    )
    assert cli.main(["check", "--data", str(tmp_path), "--quick"]) == 0
    assert seen["c"] is False
