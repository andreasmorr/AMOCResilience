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

# Each entry: (box_measure, boussinesq_measure, climberx_measure, plasim_column, ylabel, title)
# Use None where a model does not provide that measure.
PANELS = [
    (
        "amoc_strength_sv",
        "amoc_strength",
        "amoc_strength",
        "mean_amoc_strength_Sv",
        "AMOC strength (Sv)",
        "AMOC strength vs CO\u2082",
    ),
    (
        "mean_convergence_time",
        "mean_convergence_time",
        "mean_convergence_time",
        "mean_conv_time_yr",
        "Mean convergence time (yr)",
        "Convergence time vs CO\u2082",
    ),
    (
        "minimal_critical_shock_magnitude",
        None,
        None,
        "mean_edge_dist",
        "Edge\u2013attractor distance (EOF units)",
        "Critical shock / edge distance vs CO\u2082",
    ),
    (
        "basin_stability",
        None,
        None,
        "ellipsoid_volume_1sigma",
        "Ellipsoid volume (1\u03c3)",
        "Basin stability / variability vs CO\u2082",
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


def main() -> None:
    apply_style()

    df_box       = load_box_model()
    df_plasim    = load_plasim()
    df_boussinesq = load_boussinesq()
    df_climberx  = load_climberx()

    if df_box is None and df_plasim is None and df_boussinesq is None and df_climberx is None:
        print("No data available for synthesis figure. Generate model data first.")
        sys.exit(1)

    # ── Build 2×2 figure ─────────────────────────────────────────────────────
    fig, axes = plt.subplots(
        2, 2,
        figsize=(FIGURE_WIDTH, 5.5),
        constrained_layout=True,
    )
    axes_flat = axes.flatten()

    for panel_idx, (box_measure, bous_measure, cx_measure, plasim_col, ylabel, panel_title) in enumerate(PANELS):
        ax = axes_flat[panel_idx]

        # ── Box model line ────────────────────────────────────────────────────
        if df_box is not None:
            sub = df_box[
                (df_box["measure"] == box_measure) &
                (df_box["attractor"] == "on")
            ].sort_values("co2_ppm")

            if not sub.empty:
                ax.plot(
                    sub["co2_ppm"].values,
                    sub["value"].values,
                    color=COL_ON,
                    lw=1.8,
                    label="3-box model (on)",
                    zorder=3,
                )
            else:
                # Try without attractor filter (e.g. amoc_strength_sv)
                sub_all = df_box[df_box["measure"] == box_measure].sort_values("co2_ppm")
                if not sub_all.empty:
                    ax.plot(
                        sub_all["co2_ppm"].values,
                        sub_all["value"].values,
                        color=COL_ON,
                        lw=1.8,
                        label="3-box model",
                        zorder=3,
                    )

        # ── PlaSim points ─────────────────────────────────────────────────────
        if df_plasim is not None and plasim_col in df_plasim.columns:
            for state, color, marker in [
                ("AMOC-on",  COL_ON,  "^"),
                ("AMOC-off", COL_OFF, "v"),
            ]:
                sub_p = df_plasim[df_plasim["state"] == state].dropna(
                    subset=["co2_ppm", plasim_col]
                )
                if not sub_p.empty:
                    ax.scatter(
                        sub_p["co2_ppm"].values,
                        sub_p[plasim_col].values,
                        color=color,
                        marker=marker,
                        s=50,
                        zorder=5,
                        label=f"PlaSim ({state})",
                        edgecolors="white",
                        linewidths=0.5,
                    )

        # ── Boussinesq line ───────────────────────────────────────────────────
        if df_boussinesq is not None and bous_measure is not None:
            sub_b = df_boussinesq[
                (df_boussinesq["measure"] == bous_measure) &
                (df_boussinesq["attractor"] == "on")
            ].sort_values("co2_ppm")
            if not sub_b.empty:
                ax.plot(
                    sub_b["co2_ppm"].values,
                    sub_b["value"].values,
                    color=COL_BOUS,
                    lw=1.8,
                    linestyle="--",
                    label="Boussinesq (on)",
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
                    label="CLIMBER-X (on)",
                    edgecolors="white",
                    linewidths=0.5,
                )

        # ── Reference lines ───────────────────────────────────────────────────
        ax.axvline(280, color="gray", lw=0.8, ls=":", alpha=0.7)
        ax.axvline(560, color="gray", lw=0.8, ls=":", alpha=0.7)

        ax.set_title(panel_title, fontsize=8)
        ax.set_xlabel("CO\u2082 concentration (ppm)", fontsize=8)
        ax.set_ylabel(ylabel, fontsize=8)
        ax.tick_params(labelsize=7)

        if panel_idx == 0:
            # Add legend in first panel only
            ax.legend(loc="upper right", fontsize=7, framealpha=0.8)

    # Add 1×CO₂ / 2×CO₂ annotations to first panel
    ax0 = axes_flat[0]
    ymin, ymax = ax0.get_ylim()
    y_text = ymin + 0.02 * (ymax - ymin)
    ax0.text(280 + 5, y_text, "1\u00d7CO\u2082", fontsize=7, color="gray", va="bottom")
    ax0.text(560 + 5, y_text, "2\u00d7CO\u2082", fontsize=7, color="gray", va="bottom")

    # Panel labels
    for idx, label in enumerate(["(a)", "(b)", "(c)", "(d)"]):
        ax = axes_flat[idx]
        ax.text(0.03, 0.97, label,
                transform=ax.transAxes,
                fontsize=9, fontweight="bold",
                va="top", ha="left")

    # Shared legend (model identification) — placed below figure
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    legend_elements = [
        Line2D([0], [0], color=COL_ON,      lw=1.8,             label="3-box model (on)"),
        Line2D([0], [0], color=COL_BOUS,    lw=1.8, ls="--",    label="Boussinesq (on)"),
        Line2D([0], [0], color=COL_CLIMBERX, lw=0,  marker="s", markersize=6,
               label="CLIMBER-X (on)"),
        Line2D([0], [0], color=COL_ON,       lw=0,  marker="^", markersize=6,
               label="PlaSim (AMOC-on)"),
        Line2D([0], [0], color=COL_OFF,      lw=0,  marker="v", markersize=6,
               label="PlaSim (AMOC-off)"),
    ]
    fig.legend(
        handles=legend_elements,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.04),
        ncol=5,
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
