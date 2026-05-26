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
import matplotlib.gridspec as gridspec

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

# Figure width for double-column LaTeX (inches)
FIGURE_WIDTH = 7.2

# ---------------------------------------------------------------------------
# Layout helpers
# ---------------------------------------------------------------------------

def make_paper_figure(
    ncols: int = 2,
    top_height: float = 1.4,
    bottom_height: float = 3.0,
) -> tuple:
    """
    Create the standard 4-panel (2-row × ncols) paper figure.

    Returns
    -------
    fig : matplotlib.figure.Figure
    axes_top : list of ncols Axes  (AMOC strength vs time)
    axes_bottom : list of ncols Axes  (2-D phase space, square aspect)
    """
    total_height = top_height + bottom_height + 0.35 * FIGURE_WIDTH
    fig = plt.figure(figsize=(FIGURE_WIDTH, total_height))

    gs = gridspec.GridSpec(
        2, ncols,
        figure=fig,
        height_ratios=[top_height, bottom_height],
        hspace=0.35,
        wspace=0.08,
        left=0.10,
        right=0.97,
        top=0.93,
        bottom=0.08,
    )

    axes_top    = [fig.add_subplot(gs[0, c]) for c in range(ncols)]
    axes_bottom = [fig.add_subplot(gs[1, c]) for c in range(ncols)]

    # Square aspect for phase-space panels
    for ax in axes_bottom:
        ax.set_aspect("equal", adjustable="box")

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
