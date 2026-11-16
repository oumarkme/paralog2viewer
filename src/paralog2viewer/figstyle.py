# SPDX-License-Identifier: GPL-3.0-or-later
"""Matplotlib settings shared by the read-pileup and RCCX figures."""

from __future__ import annotations

import functools

from matplotlib import font_manager

# Text size of all figures; 7 pt is the lab minimum for figure text.
FONT_PT = 7


@functools.cache
def font_family() -> str:
    """Return Liberation Sans (Arial metrics) if installed, else DejaVu Sans."""
    families = {f.name for f in font_manager.fontManager.ttflist}
    return "Liberation Sans" if "Liberation Sans" in families else "DejaVu Sans"


def base_rc() -> dict[str, object]:
    """Return the rc settings both figure types use (font and line width)."""
    return {
        "font.family": font_family(),
        "font.size": FONT_PT,
        "axes.linewidth": 0.5,
    }
