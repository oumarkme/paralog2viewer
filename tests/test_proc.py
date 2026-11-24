# SPDX-License-Identifier: GPL-3.0-or-later
"""Tests for the external-tool helpers."""

from pathlib import Path

import pytest

from paralog2viewer.proc import ToolError, pipeline, run


def test_run_writes_stdout(tmp_path: Path) -> None:
    out = tmp_path / "out.txt"
    run(["echo", "hello"], stdout=out)
    assert out.read_text() == "hello\n"


def test_run_raises_on_failure() -> None:
    with pytest.raises(ToolError):
        run(["false"])


def test_pipeline_connects_commands(tmp_path: Path) -> None:
    out = tmp_path / "out.txt"
    pipeline([["printf", "b\\na\\n"], ["sort"]], stdout=out)
    assert out.read_text() == "a\nb\n"


def test_pipeline_detects_failure_in_first_command() -> None:
    with pytest.raises(ToolError):
        pipeline([["false"], ["cat"]])
