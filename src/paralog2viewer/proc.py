# SPDX-License-Identifier: GPL-3.0-or-later
"""Helpers for running external tools.

Every module that calls an external program uses these functions, so that
commands are logged, failures raise one error type, and shell pipelines keep
``pipefail`` semantics.
"""

from __future__ import annotations

import logging
import os
import shlex
import subprocess
from collections.abc import Mapping, Sequence
from pathlib import Path

LOGGER = logging.getLogger("paralog2viewer")


class ToolError(RuntimeError):
    """An external tool exited with a non-zero status."""


def run(
    cmd: Sequence[str | os.PathLike[str]],
    *,
    cwd: Path | None = None,
    env: Mapping[str, str] | None = None,
    stdout: Path | None = None,
    stderr: Path | None = None,
) -> None:
    """Run one command and raise ToolError if it fails.

    Args:
        cmd: Program and arguments.
        cwd: Working directory.
        env: Extra environment variables added to the current environment.
        stdout: File that receives standard output; inherited when None.
        stderr: File that receives standard error; inherited when None.
    """
    argv = [str(part) for part in cmd]
    LOGGER.info("run: %s", shlex.join(argv))
    full_env = {**os.environ, **(env or {})}
    out = open(stdout, "wb") if stdout else None  # noqa: SIM115
    err = open(stderr, "wb") if stderr else None  # noqa: SIM115
    try:
        proc = subprocess.run(
            argv, cwd=cwd, env=full_env, stdout=out, stderr=err, check=False
        )
    finally:
        for handle in (out, err):
            if handle:
                handle.close()
    if proc.returncode != 0:
        raise ToolError(
            f"{argv[0]} exited with {proc.returncode}: {shlex.join(argv)}"
        )


def pipeline(
    commands: Sequence[Sequence[str | os.PathLike[str]]],
    *,
    cwd: Path | None = None,
    env: Mapping[str, str] | None = None,
    stdout: Path | None = None,
) -> None:
    """Run commands connected by pipes and raise ToolError if any fails.

    Args:
        commands: Each element is one command (program and arguments).
        cwd: Working directory.
        env: Extra environment variables added to the current environment.
        stdout: File that receives the last command's output; inherited
            when None.
    """
    argvs = [[str(part) for part in cmd] for cmd in commands]
    LOGGER.info("pipeline: %s", " | ".join(shlex.join(a) for a in argvs))
    full_env = {**os.environ, **(env or {})}
    out = open(stdout, "wb") if stdout else None  # noqa: SIM115
    procs: list[subprocess.Popen[bytes]] = []
    try:
        previous = None
        for index, argv in enumerate(argvs):
            last = index == len(argvs) - 1
            proc = subprocess.Popen(
                argv,
                cwd=cwd,
                env=full_env,
                stdin=previous.stdout if previous else None,
                stdout=out if last else subprocess.PIPE,
            )
            if previous and previous.stdout:
                previous.stdout.close()
            procs.append(proc)
            previous = proc
        codes = [proc.wait() for proc in procs]
    finally:
        if out:
            out.close()
    for argv, code in zip(argvs, codes, strict=True):
        if code != 0:
            raise ToolError(f"{argv[0]} exited with {code} in pipeline")
