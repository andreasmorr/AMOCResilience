"""
synthesis_figure.py  –  Cross-model AMOC resilience vs CO2 figure.

Reads:
    AMOCBox/data/paper/resilience_vs_co2_boxmodel.csv
        columns: co2_ppm, t_param, measure, value, attractor

    AMOCPlaSim/data/plasim/resilience_metrics.csv
        columns: co2_ppm, state, mean_conv_time_yr, mean_edge_dist,
                 ellipsoid_volume_1sigma, mean_amoc_strength_Sv

    AMOCBoussinesq: placeholder (see TODO below)
    AMOCClimberX:   placeholder (see TODO below)

Output: synthesis_figure.pdf

Run from the AMOCResilience umbrella directory:
    python synthesis_figure.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

UMBRELLA  = Path(__file__).resolve().parent
sys.path.insert(0, str(UMBRELLA))

from amoc_plot_style import (
    COL_ON, COL_OFF, COL_EDGE,
    apply_style, savefig_pdf,
)

FIGURE_WIDTH = 10.0   # inches

BOX_CSV        = UMBRELLA / "AMOCBox"        / "data" / "paper" / "resilience_vs_co2_boxmodel.csv"
PLASIM_CSV     = UMBRELLA / "AMOCPlaSim"    / "data" / "plasim" / "resilience_metrics.csv"
BOUSSINESQ_CSV = UMBRELLA / "AMOCBoussinesq" / "data" / "paper" / "resilience_vs_gamma_boussinesq.csv"
CLIMBERX_CSV   = UMBRELLA / "AMOCClimberX"  / "data" / "paper" / "resilience_vs_co2_climberx.csv"


# ---------------------------------------------------------------------------
# Panel definitions
# ---------------------------------------------------------------------------

# Top panel (full width): AMOC strength
AMOC_PANEL = (
    "amoc_strength_sv",
    "amoc_strength",
    "amoc_strength",
    "mean_amoc_strength_Sv",
    "AMOC strength (Sv)",
    "AMOC strength",
)

# Each entry: (box_measure, boussinesq_measure, climberx_measure, plasim_column, ylabel, title)
# Use None where a model does not provide that measure.
PANELS = [
    (
        "mean_convergence_time",
        "mean_convergence_time",
        "mean_convergence_time",
        "mean_conv_time_yr",
        "Mean convergence time (yr)",
        "Convergence time",
    ),
    (
        "minimal_critical_shock_magnitude",
        None,
        None,
        "mean_edge_dist",
        "Edge\u2013attractor distance (EOF units)",
        "Critical shock / edge distance",
    ),
    (
        "characteristic_return_time",
        "characteristic_return_time",
        "characteristic_return_time",
        "ellipse_long_axis_1sigma",
        "Char. return time / ellipse long axis",
        "Characteristic return time",
    ),
    (
        "basin_stability",
        "basin_volume",
        None,
        None,
        "Basin volume / stability",
        "Basin volume",
    ),
]

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------

def load_box_model() -> pd.DataFrame | None:
    if not BOX_CSV.exists():
        print(f"[box model] CSV not found: {BOX_CSV}")
        print("  Run amoc3box_co2_continuation.jl first.")
        return None
    df = pd.read_csv(BOX_CSV)
    return df


def load_plasim() -> pd.DataFrame | None:
    if not PLASIM_CSV.exists():
        print(f"[PlaSim] CSV not found: {PLASIM_CSV}")
        return None
    return pd.read_csv(PLASIM_CSV)


def load_boussinesq() -> pd.DataFrame | None:
    if not BOUSSINESQ_CSV.exists():
        print(f"[Boussinesq] CSV not found: {BOUSSINESQ_CSV}")
        return None
    return pd.read_csv(BOUSSINESQ_CSV)


def load_climberx() -> pd.DataFrame | None:
    if not CLIMBERX_CSV.exists():
        print(f"[CLIMBER-X] CSV not found: {CLIMBERX_CSV}")
        return None
    return pd.read_csv(CLIMBERX_CSV)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

COL_BOUS    = "#7B3F00"   # brown for Boussinesq
COL_CLIMBERX = "#4D0099"  # purple for CLIMBER-X
COL_PLASIM   = "#D95F02"  # orange for PlaSim


def _plot_panel(ax, box_measure, bous_measure, cx_measure, plasim_col,
                ylabel, panel_title, df_box, df_plasim, df_boussinesq, df_climberx,
                xlabel=False):
    """Plot a single resilience-measure panel onto *ax*."""

    # ── Box model scatter ─────────────────────────────────────────────────
    if df_box is not None and box_measure is not None:
        sub = df_box[
            (df_box["measure"] == box_measure) &
            (df_box["attractor"] == "on")
        ].sort_values("co2_ppm")

        if not sub.empty:
            ax.scatter(
                sub["co2_ppm"].values,
                sub["value"].values,
                color=COL_ON,
                marker="o",
                s=20,
                label="3-box model",
                zorder=3,
            )
        else:
            # Try without attractor filter (e.g. amoc_strength_sv)
            sub_all = df_box[df_box["measure"] == box_measure].sort_values("co2_ppm")
            if not sub_all.empty:
                ax.scatter(
                    sub_all["co2_ppm"].values,
                    sub_all["value"].values,
                    color=COL_ON,
                    marker="o",
                    s=20,
                    label="3-box model",
                    zorder=3,
                )

    # ── PlaSim points ─────────────────────────────────────────────────────
    if df_plasim is not None and plasim_col is not None and plasim_col in df_plasim.columns:
        sub_p = df_plasim[df_plasim["state"] == "AMOC-on"].dropna(
            subset=["co2_ppm", plasim_col]
        )
        if not sub_p.empty:
            ax.scatter(
                sub_p["co2_ppm"].values,
                sub_p[plasim_col].values,
                color=COL_PLASIM,
                marker="^",
                s=50,
                zorder=5,
                label="PlaSim",
                edgecolors="white",
                linewidths=0.5,
            )

    # ── Boussinesq scatter ────────────────────────────────────────────────
    if df_boussinesq is not None and bous_measure is not None:
        sub_b = df_boussinesq[
            (df_boussinesq["measure"] == bous_measure) &
            (df_boussinesq["attractor"] == "on")
        ].sort_values("co2_ppm")
        if not sub_b.empty:
            ax.scatter(
                sub_b["co2_ppm"].values,
                sub_b["value"].values,
                color=COL_BOUS,
                marker="D",
                s=20,
                label="Boussinesq",
                zorder=3,
            )

    # ── CLIMBER-X scatter ─────────────────────────────────────────────────
    if df_climberx is not None and cx_measure is not None:
        sub_cx = df_climberx[
            (df_climberx["measure"] == cx_measure) &
            (df_climberx["attractor"] == "on")
        ].sort_values("co2_ppm")
        if not sub_cx.empty:
            ax.scatter(
                sub_cx["co2_ppm"].values,
                sub_cx["value"].values,
                color=COL_CLIMBERX,
                marker="s",
                s=50,
                zorder=5,
                label="CLIMBER-X",
                edgecolors="white",
                linewidths=0.5,
            )

    ax.set_title(panel_title, fontsize=8)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.tick_params(labelsize=7)
    if xlabel:
        ax.set_xlabel("CO\u2082 concentration (ppm)", fontsize=8)


def main() -> None:
    apply_style()

    df_box       = load_box_model()
    df_plasim    = load_plasim()
    df_boussinesq = load_boussinesq()
    df_climberx  = load_climberx()

    if df_box is None and df_plasim is None and df_boussinesq is None and df_climberx is None:
        print("No data available for synthesis figure. Generate model data first.")
        sys.exit(1)

    # ── Build figure ─────────────────────────────────────────────────────────
    # Layout: AMOC strength spans full width (top row),
    # then 2×2 grid for the four resilience measure panels.
    from matplotlib.gridspec import GridSpec

    ncols = 2
    n_res_panels = len(PANELS)               # 4
    n_res_rows   = (n_res_panels + 1) // ncols  # 2
    nrows_total  = 1 + n_res_rows            # 3

    fig = plt.figure(figsize=(FIGURE_WIDTH, 3.0 + 2.5 * n_res_rows),
                     constrained_layout=True)
    gs  = GridSpec(nrows_total, ncols, figure=fig,
                   height_ratios=[1.0] + [1.0] * n_res_rows)

    # Top panel – AMOC strength
    ax_amoc = fig.add_subplot(gs[0, :])
    _plot_panel(
        ax_amoc,
        *AMOC_PANEL,
        df_box, df_plasim, df_boussinesq, df_climberx,
        xlabel=False,
    )
    ax_amoc.text(0.03, 0.97, "(a)",
                 transform=ax_amoc.transAxes,
                 fontsize=9, fontweight="bold", va="top", ha="left")

    # Resilience-measure panels – 2×2 grid
    panel_labels = ["(b)", "(c)", "(d)", "(e)"]
    for panel_idx, (box_measure, bous_measure, cx_measure, plasim_col, ylabel, panel_title) in enumerate(PANELS):
        row = 1 + panel_idx // ncols
        col = panel_idx % ncols
        is_last_row = (row == nrows_total - 1)

        ax = fig.add_subplot(gs[row, col])
        _plot_panel(
            ax,
            box_measure, bous_measure, cx_measure, plasim_col,
            ylabel, panel_title,
            df_box, df_plasim, df_boussinesq, df_climberx,
            xlabel=is_last_row,
        )
        ax.text(0.03, 0.97, panel_labels[panel_idx],
                transform=ax.transAxes,
                fontsize=9, fontweight="bold", va="top", ha="left")

    # Shared legend – placed below figure
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color=COL_ON,       lw=0, marker="o", markersize=6, label="3-box model"),
        Line2D([0], [0], color=COL_BOUS,     lw=0, marker="D", markersize=6, label="Boussinesq"),
        Line2D([0], [0], color=COL_CLIMBERX, lw=0, marker="s", markersize=6, label="CLIMBER-X"),
        Line2D([0], [0], color=COL_PLASIM,   lw=0, marker="^", markersize=6, label="PlaSim"),
    ]
    fig.legend(
        handles=legend_elements,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.04),
        ncol=4,
        fontsize=7,
        framealpha=0.8,
    )

    out_path = UMBRELLA / "plots" / "synthesis_figure.png"
    (UMBRELLA / "plots").mkdir(exist_ok=True)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    print(f"Figure saved: {out_path}")
    plt.close(fig)


if __name__ == "__main__":
    main()
