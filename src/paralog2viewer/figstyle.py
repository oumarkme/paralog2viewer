# SPDX-License-Identifier: GPL-3.0-or-later
"""Matplotlib settings shared by the read-pileup and RCCX figures."""

from __future__ import annotations

import functools

from matplotlib import font_manager

# Text size of all figures; 7 pt is the lab minimum for figure text.
FONT_PT = 7
# Figure resolution; the HTML summary embeds the files as written.
FIGURE_DPI = 150


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


# One color per haplotype role in every figure (Okabe-Ito and grays).
ROLE_COLOR = {
    "gene": "#0072B2",
    "paralog_gene": "#E69F00",
    "pseudogene": "#D55E00",
    "fusion": "#CC79A7",
    "other": "#595959",
    "unresolved": "#9A9A9A",
}
# Allele carried by haplotypes of more than one role.
MIXED_COLOR = "#009E73"


def tint(color: str, fraction: float) -> str:
    """Return ``color`` blended with white (0 = white, 1 = ``color``)."""
    rgb = [int(color[i : i + 2], 16) for i in (1, 3, 5)]
    blended = [round(255 - (255 - c) * fraction) for c in rgb]
    return "#" + "".join(f"{c:02X}" for c in blended)
