# SPDX-License-Identifier: GPL-3.0-or-later
from pathlib import Path

import pytest

from paralog2viewer import config


def test_parse_formats_both_and_single():
    assert config.parse_formats("vcf,mviewer") == ("vcf", "mviewer")
    assert config.parse_formats("MViewer") == ("mviewer",)


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


def test_settings_defaults_follow_gs4_budget():
    s = _settings()
    assert s.threads == 40 and s.max_memory_gb == 800
    assert s.resume and not s.dry_run and s.cpus is None


@pytest.mark.parametrize(
    "kwargs",
    [
        {"threads": 0},
        {"max_memory_gb": 0},
        {"vcf_annotation": "all"},
        {"cpus": (0, 1), "threads": 3},
    ],
)
def test_settings_validation(kwargs):
    with pytest.raises(ValueError):
        _settings(**kwargs)
