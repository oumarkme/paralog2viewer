# SPDX-License-Identifier: GPL-3.0-or-later
import os
from pathlib import Path

import pytest

from paralog2viewer import config


def test_parse_formats_both_and_single():
    assert config.parse_formats("vcf,mviewer") == ("vcf", "mviewer")
    assert config.parse_formats("MViewer") == ("mviewer",)
    assert config.parse_formats("vcf,Summary") == ("vcf", "summary")


def test_parse_formats_dedup_keeps_order():
    assert config.parse_formats("vcf, vcf") == ("vcf",)


@pytest.mark.parametrize("text", ["", " , ", "vcf,bam"])
def test_parse_formats_invalid(text):
    with pytest.raises(ValueError):
        config.parse_formats(text)


def test_data_env_fallback_and_cli_wins():
    env = {config.DATA_ENV_VAR: "/env"}
    assert str(config.resolve_data_dir(None, env)) == "/env"
    assert str(config.resolve_data_dir("/cli", env)) == "/cli"
    assert config.resolve_data_dir(None, {}) is None


def _settings(**kwargs):
    return config.RunSettings(
        bam=Path("a.bam"),
        reference=Path("r.fa"),
        data=Path("d"),
        outdir=Path("o"),
        **kwargs,
    )


def test_settings_defaults_use_affinity_and_memory():
    s = _settings()
    assert s.threads == len(os.sched_getaffinity(0)) and s.memory_gb > 0
    assert s.resume and not s.dry_run


def test_threads_are_clamped_to_affinity():
    usable = config.default_threads()
    assert config.clamp_threads(None) == (usable, False)
    assert config.clamp_threads(usable + 1) == (usable, True)
    assert config.clamp_threads(1) == (1, False)


def test_default_memory_takes_smaller_of_cgroup_and_available(tmp_path):
    meminfo = tmp_path / "meminfo"
    meminfo.write_text("MemTotal: 9 kB\nMemAvailable: 100000000 kB\n")
    v2, v1 = tmp_path / "v2", tmp_path / "v1"
    v2.write_text("max\n")
    v1.write_text(str(1 << 62))
    assert config.default_memory_gb(v2, v1, meminfo) == pytest.approx(
        80.0 * 1.024
    )
    v2.write_text("50000000000\n")
    assert config.default_memory_gb(v2, v1, meminfo) == pytest.approx(40.0)
    v1.write_text("20000000000\n")
    assert config.default_memory_gb(v2, v1, meminfo) == pytest.approx(16.0)
    with pytest.raises(ValueError, match="--memory"):
        config.default_memory_gb(tmp_path / "x", tmp_path / "y", tmp_path / "z")


@pytest.mark.parametrize(
    "kwargs",
    [
        {"threads": 0},
        {"memory_gb": 0},
        {"vcf_annotation": "all"},
    ],
)
def test_settings_validation(kwargs):
    with pytest.raises(ValueError):
        _settings(**kwargs)
