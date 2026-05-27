"""
amoc_plot_style.py  –  Shared matplotlib style for all AMOC resilience paper figures.

Usage in any submodule script:
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).parent.parent))  # adjust depth as needed
    from amoc_plot_style import (
        COL_ON, COL_OFF, COL_EDGE,
        BASIN_ON_FILL, BASIN_OFF_FILL, TRAJ_COLORS,
        make_paper_figure, add_panel_label, savefig_pdf,
    )
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Global rcParams — call apply_style() once at import time
# ---------------------------------------------------------------------------

def apply_style() -> None:
    """Apply paper-quality rcParams globally."""
    matplotlib.rcParams.update({
        # Font
        "font.family":       "serif",
        "font.serif":        ["DejaVu Serif", "Georgia", "Times New Roman", "serif"],
        "font.size":         9,
        "axes.titlesize":    9,
        "axes.labelsize":    9,
        "xtick.labelsize":   8,
        "ytick.labelsize":   8,
        "legend.fontsize":   8,
        # Figure
        "figure.dpi":        150,
        "savefig.dpi":       300,
        "savefig.format":    "pdf",
        "savefig.bbox":      "tight",
        "savefig.pad_inches": 0.02,
        # Lines / patches
        "lines.linewidth":   1.2,
        "axes.linewidth":    0.7,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "xtick.minor.width": 0.5,
        "ytick.minor.width": 0.5,
        # PDF backend
        "pdf.fonttype":      42,   # TrueType fonts in PDF (required by many journals)
        "ps.fonttype":       42,
    })


apply_style()

# ---------------------------------------------------------------------------
# Color constants
# ---------------------------------------------------------------------------

COL_ON   = "#2166AC"   # AMOC-on attractor / low CO2 (blue)
COL_OFF  = "#D6604D"   # AMOC-off attractor / high CO2 (orange-red)
COL_EDGE = "#1A1A1A"   # edge state (near-black)

BASIN_ON_FILL  = "#D1E5F0"   # pale blue fill for on-basin
BASIN_OFF_FILL = "#FDDBC7"   # pale orange fill for off-basin

TRAJ_COLORS = ["#1B7837", "#762A83", "#E08214"]  # green, purple, orange

# ---------------------------------------------------------------------------
# Hard-coded panel size constants (all in inches)
# ---------------------------------------------------------------------------

PANEL_SIZE   = 3.0    # width = height of the square phase-space panels
TOP_HEIGHT   = 1.2    # height of the time-series top row
LEFT_MARGIN  = 0.65   # room for y-axis label
RIGHT_MARGIN = 0.15
TOP_MARGIN   = 0.35
BOT_MARGIN   = 0.50
H_GAP        = 0.45   # vertical gap between the two rows
W_GAP        = 0.05   # horizontal gap between columns (shared y → minimal)

# ---------------------------------------------------------------------------
# Layout helpers
# ---------------------------------------------------------------------------

def make_paper_figure(
    ncols: int = 2,
    panel_size: float = PANEL_SIZE,
    top_height: float = TOP_HEIGHT,
    left_margin: float = LEFT_MARGIN,
    right_margin: float = RIGHT_MARGIN,
    top_margin: float = TOP_MARGIN,
    bot_margin: float = BOT_MARGIN,
    h_gap: float = H_GAP,
    w_gap: float = W_GAP,
    sharey_bottom: bool = True,
    sharey_top: bool = True,
) -> tuple:
    """
    Create the standard 4-panel (2-row × ncols) paper figure with hard-coded
    absolute panel sizes.

    Layout (all dimensions in inches):
      - Top row:    ncols panels, each panel_size wide and top_height tall
      - Bottom row: ncols perfectly square panels (panel_size × panel_size)
      - Bottom panels share the y-axis when sharey_bottom=True

    Returns
    -------
    fig : matplotlib.figure.Figure
    axes_top    : list of ncols Axes  (time-series row)
    axes_bottom : list of ncols Axes  (square phase-space row, shared y)
    """
    fig_w = left_margin + ncols * panel_size + (ncols - 1) * w_gap + right_margin
    fig_h = top_margin + top_height + h_gap + panel_size + bot_margin
    fig = plt.figure(figsize=(fig_w, fig_h))

    def fx(x): return x / fig_w   # inches → figure x-fraction
    def fy(y): return y / fig_h   # inches → figure y-fraction

    # ── Bottom row (create first so sharey reference axis exists) ────────────
    axes_bottom = []
    for c in range(ncols):
        sharey_ax = axes_bottom[0] if (sharey_bottom and c > 0) else None
        ax = fig.add_axes(
            [fx(left_margin + c * (panel_size + w_gap)),
             fy(bot_margin),
             fx(panel_size),
             fy(panel_size)],
            sharey=sharey_ax,
        )
        ax.set_aspect("equal", adjustable="datalim")
        axes_bottom.append(ax)

    # ── Top row ──────────────────────────────────────────────────────────────
    axes_top = []
    for c in range(ncols):
        sharey_ax = axes_top[0] if (sharey_top and c > 0) else None
        ax = fig.add_axes(
            [fx(left_margin + c * (panel_size + w_gap)),
             fy(bot_margin + panel_size + h_gap),
             fx(panel_size),
             fy(top_height)],
            sharey=sharey_ax,
        )
        axes_top.append(ax)

    return fig, axes_top, axes_bottom


def add_panel_label(
    ax: "matplotlib.axes.Axes",
    label: str,
    x: float = 0.03,
    y: float = 0.97,
) -> None:
    """
    Add a panel label (e.g. '(a)') in axes-fraction coordinates.

    Parameters
    ----------
    ax : Axes
    label : str  e.g. '(a)', '(b)', …
    x, y : float  position in axes-fraction units (default: upper-left corner)
    """
    ax.text(
        x, y, label,
        transform=ax.transAxes,
        fontsize=9,
        fontweight="bold",
        va="top",
        ha="left",
    )


def savefig_pdf(fig: "matplotlib.figure.Figure", path: str | Path) -> None:
    """
    Save figure as a tight-layout PDF.

    Parameters
    ----------
    fig : Figure
    path : str or Path  destination file (created with parents if needed)
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, format="pdf", bbox_inches="tight", dpi=300)
    print(f"Figure saved: {path}")
