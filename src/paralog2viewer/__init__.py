# SPDX-License-Identifier: GPL-3.0-or-later
"""paralog2viewer: per-haplotype calls and figures from Paraphase results."""

from importlib import metadata
from pathlib import Path


def _read_version() -> str:
    """Returns the version from package metadata or the VERSION file."""
    try:
        return metadata.version("paralog2viewer")
    except metadata.PackageNotFoundError:
        # Source checkout without installation: VERSION sits at the repo root.
        return (Path(__file__).parents[2] / "VERSION").read_text().strip()


__version__ = _read_version()
