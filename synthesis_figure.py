"""
synthesis_figure.py  –  Cross-model AMOC resilience vs CO2 figure.

Layout: 2 columns × 3 rows.  Top row holds the two model-state panels —
(a) AMOC strength and (b) the Tropical−Northern equilibrium salinity gradient
between the two readout boxes — and the four resilience-measure panels (c)–(f)
fill the lower two rows.  The three panels in each column share a CO2 x-axis.

Salinity-gradient sources: Boussinesq and CLIMBER-X carry it as the
"salinity_gradient" measure in their resilience CSVs; the box model and PlaSim
expose both readout boxes in per-CO2 auxiliary files (attractors_*.csv /
state_means_*.csv), from which the gradient is computed here.

Reads:
    AMOCBox/data/paper/resilience_vs_co2_boxmodel.csv
        columns: co2_ppm, t_param, measure, value, attractor

    AMOCPlaSim/data/results/resilience_metrics.csv
        columns: co2_ppm, state, mean_conv_time_yr, mean_edge_dist,
                 ellipsoid_volume_1sigma, ellipse_long_axis_1sigma,
                 local_resilience, mean_amoc_strength_Sv

    AMOCBoussinesq/data/paper/resilience_vs_gamma_boussinesq.csv
        columns: co2_ppm, gamma, measure, value, attractor, units, ...
        Paper values use the shared plotting-unit conversions; native values
        and assumptions are retained in the export and its provenance JSON.

    AMOCClimberX/data/paper/resilience_vs_co2_climberx.csv
        columns: co2_ppm, measure, value, attractor

Outputs: plots/synthesis_figure.png, .pdf and .json

Run from the AMOCResilience umbrella directory:
    python synthesis_figure.py
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
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

_path_parser = argparse.ArgumentParser(add_help=False)
_path_parser.add_argument("--umbrella", type=Path, default=Path(__file__).resolve().parent)
_path_args, _ = _path_parser.parse_known_args()
UMBRELLA = _path_args.umbrella.resolve()
sys.path.insert(0, str(UMBRELLA))

from amoc_plot_style import (
    COL_ON,
    COL_BOUS, COL_CLIMBERX, COL_PLASIM,
    apply_style, add_panel_label,
)

FIGURE_WIDTH = 6.30   # inches (single-column \textwidth: A4, 2.5 cm margins)

BOX_CSV        = UMBRELLA / "AMOCBox"        / "data" / "paper" / "resilience_vs_co2_boxmodel.csv"
PLASIM_CSV     = UMBRELLA / "AMOCPlaSim"    / "data" / "results" / "resilience_metrics.csv"
BOUSSINESQ_CSV = UMBRELLA / "AMOCBoussinesq" / "data" / "paper" / "resilience_vs_gamma_boussinesq.csv"
CLIMBERX_CSV   = UMBRELLA / "AMOCClimberX"  / "data" / "paper" / "resilience_vs_co2_climberx.csv"

# Maximum parameter values; rows beyond these are excluded (None = no cutoff)
BOX_MAX_CO2      = None
BOX_MAX_GAMMA    = None
BOUS_MAX_GAMMA   = None
BOUS_MAX_CO2     = None
PLASIM_MAX_CO2   = None
CLIMBERX_MAX_CO2 = None

# Minimum parameter values; rows below these are excluded (None = no cutoff)
BOX_MIN_CO2      = None
BOX_MIN_GAMMA    = None

# Box-model CO2 levels (ppm) shown only in the basin-volume and minimal-critical-
# shock panels (e, f). 905.8 and 906.64 ppm (t = 2.235, 2.238) lie just before
# the Hopf bifurcation (~908 ppm) where the on-state basin collapses; they are
# kept out of the other panels. Attractor files are tagged with rounded CO2
# (906, 907 ppm), so matching uses BOX_CO2_MATCH_TOL.
BOX_BASIN_ONLY_CO2 = (905.8, 906.64)
BOX_BASIN_ONLY_MEASURES = {"basin_stability", "minimal_critical_shock_magnitude"}
BOX_CO2_MATCH_TOL = 0.5
BOUS_MIN_GAMMA   = None
BOUS_MIN_CO2     = None
PLASIM_MIN_CO2   = None
CLIMBERX_MIN_CO2 = None


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
# Order: (b) local resilience, (c) convergence time, (d) basin volume, (e) critical shock.
# b+d share x-axis (left column); c+e share x-axis (right column).
PANELS = [
    (
        "local_resilience",
        "local_resilience",
        "local_resilience",
        "local_resilience",
        "Local resilience",
        "Local resilience",
    ),
    (
        "mean_convergence_time",
        "mean_convergence_time",
        "mean_convergence_time",
        "mean_conv_time_yr",
        "Convergence time",
        "Convergence time",
    ),
    (
        "basin_stability",
        "basin_volume",
        "basin_volume",
        None,
        "Basin volume",
        "Basin volume",
    ),
    (
        "minimal_critical_shock_magnitude",
        "min_critical_shock",
        "min_critical_shock",
        "mean_edge_dist",
        "Minimal critical shock",
        "Critical shock",
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
    if BOX_MIN_CO2 is not None:
        df = df[df["co2_ppm"] >= BOX_MIN_CO2]
    if BOX_MAX_CO2 is not None:
        df = df[df["co2_ppm"] <= BOX_MAX_CO2]
    if BOX_MIN_GAMMA is not None and "gamma" in df.columns:
        df = df[df["gamma"] >= BOX_MIN_GAMMA]
    if BOX_MAX_GAMMA is not None and "gamma" in df.columns:
        df = df[df["gamma"] <= BOX_MAX_GAMMA]
    return df


def load_plasim() -> pd.DataFrame | None:
    if not PLASIM_CSV.exists():
        print(f"[PlaSim] CSV not found: {PLASIM_CSV}")
        return None
    df = pd.read_csv(PLASIM_CSV)
    if PLASIM_MIN_CO2 is not None:
        df = df[df["co2_ppm"] >= PLASIM_MIN_CO2]
    if PLASIM_MAX_CO2 is not None:
        df = df[df["co2_ppm"] <= PLASIM_MAX_CO2]
    return df


def load_boussinesq() -> pd.DataFrame | None:
    if not BOUSSINESQ_CSV.exists():
        print(f"[Boussinesq] CSV not found: {BOUSSINESQ_CSV}")
        return None
    df = pd.read_csv(BOUSSINESQ_CSV)
    required = {"gamma", "co2_ppm", "measure", "value", "units", "status", "model_value", "model_units"}
    if not required.issubset(df.columns):
        raise ValueError("Use the canonical paper export with converted and native-value columns")
    if not df["status"].isin(["complete", "on_reference_unavailable"]).all():
        raise ValueError("Boussinesq experiment has unfinished parameter slices")
    if df.loc[df["status"] == "on_reference_unavailable", "value"].notna().any():
        raise ValueError("Unavailable ON-reference rows must have blank values")
    if BOUS_MIN_GAMMA is not None and "gamma" in df.columns:
        df = df[df["gamma"] >= BOUS_MIN_GAMMA]
    if BOUS_MAX_GAMMA is not None and "gamma" in df.columns:
        df = df[df["gamma"] <= BOUS_MAX_GAMMA]
    if BOUS_MIN_CO2 is not None and "co2_ppm" in df.columns:
        df = df[df["co2_ppm"] >= BOUS_MIN_CO2]
    if BOUS_MAX_CO2 is not None and "co2_ppm" in df.columns:
        df = df[df["co2_ppm"] <= BOUS_MAX_CO2]
    return df


def load_climberx() -> pd.DataFrame | None:
    if not CLIMBERX_CSV.exists():
        print(f"[CLIMBER-X] CSV not found: {CLIMBERX_CSV}")
        return None
    df = pd.read_csv(CLIMBERX_CSV)
    if CLIMBERX_MIN_CO2 is not None:
        df = df[df["co2_ppm"] >= CLIMBERX_MIN_CO2]
    if CLIMBERX_MAX_CO2 is not None:
        df = df[df["co2_ppm"] <= CLIMBERX_MAX_CO2]
    return df


# ---------------------------------------------------------------------------
# Tropical−Northern equilibrium salinity gradient (between the two readout
# boxes), per model.  Box and PlaSim expose both boxes in per-CO2 auxiliary
# files (gradient computed here); Boussinesq and CLIMBER-X carry it as the
# "salinity_gradient" measure in their resilience CSVs (loaded above), so those
# two need no separate loader — they are plotted through _plot_panel.
# ---------------------------------------------------------------------------

BOX_PAPER_DIR    = UMBRELLA / "AMOCBox"    / "data" / "paper"
PLASIM_PAPER_DIR = UMBRELLA / "AMOCPlaSim" / "data" / "results" / "paper"


def load_box_salinity_gradient() -> pd.DataFrame | None:
    """On-state Tropical−Northern salinity gradient (psu) vs CO2 from box attractors.

    Returned in the long schema (co2_ppm, measure='salinity_gradient', value,
    attractor) so it can be plotted through _plot_panel as the box-model series.
    Box salinity is stored in model units; psu = 35 + 10 * S (matching
    plotting_paper.py), so the +35 offset cancels in the difference and the
    gradient is 10 * (S_T − S_N)."""
    rows = []
    for f in sorted(BOX_PAPER_DIR.glob("attractors_*ppm.csv")):
        m = re.search(r"attractors_(\d+)ppm", f.name)
        if not m:
            continue
        d = pd.read_csv(f)
        on = d[d["state"] == "on"]
        if on.empty:
            continue
        s_n = float(on["S_N"].iloc[0])
        s_t = float(on["S_T"].iloc[0])
        rows.append(dict(co2_ppm=float(m.group(1)), measure="salinity_gradient",
                         value=10.0 * (s_t - s_n), attractor="on"))
    return pd.DataFrame(rows) if rows else None


def load_plasim_salinity_gradient() -> pd.DataFrame | None:
    """On-state Tropical−Northern salinity gradient (psu) vs CO2 from PlaSim state-means.

    Returned in the wide schema (co2_ppm, state='AMOC-on', salinity_gradient) so it
    can be plotted through _plot_panel as the PlaSim series.  x1 (Northern) and x2
    (Tropical) are already in psu, so the gradient is x2 − x1."""
    rows = []
    for f in sorted(PLASIM_PAPER_DIR.glob("state_means_*ppm.csv")):
        m = re.search(r"state_means_(\d+)ppm", f.name)
        if not m:
            continue
        d = pd.read_csv(f)
        on = d[d["state"] == "on"]
        if on.empty:
            continue
        rows.append(dict(co2_ppm=float(m.group(1)), state="AMOC-on",
                         salinity_gradient=float(on["x2"].iloc[0]) - float(on["x1"].iloc[0])))
    return pd.DataFrame(rows) if rows else None


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def _mark_line_end(ax, sub, color):
    """Hollow circle marking the last plotted on-state point.

    The marker is offset to the right by half its width (its radius) and down by
    half its diameter (its radius) via a points-based transform, so the shift is
    independent of the data scaling."""
    if sub is None or sub.empty:
        return
    from matplotlib.transforms import offset_copy
    marker_size = 45                      # points**2
    dx_pt = (marker_size ** 0.5) / 2.0    # half the marker width (= radius), in points
    trans = offset_copy(ax.transData, fig=ax.figure, x=dx_pt, y=-dx_pt, units="points")
    ax.scatter(
        [sub["co2_ppm"].values[-1]], [sub["value"].values[-1]],
        facecolors="none", edgecolors=color, s=marker_size, linewidths=1.3,
        zorder=6, transform=trans,
    )


def _add_up_arrow(ax, x, y, color, gap_frac=0.02, length_frac=0.13):
    """Upward arrow starting just above the data point (x, y), signalling an
    omitted further increase. Call after all series are plotted: the y-axis is
    extended (and then fixed) so the arrow fits inside the panel. Gap and length
    are fractions of the final y-range."""
    ax.autoscale_view()
    bottom, top = ax.get_ylim()
    need = gap_frac + length_frac + 0.02          # arrow plus a small margin
    new_top = max(top, (y - need * bottom) / (1.0 - need))
    ax.set_ylim(bottom, new_top)
    span = new_top - bottom
    ax.annotate(
        "", xy=(x, y + (gap_frac + length_frac) * span), xytext=(x, y + gap_frac * span),
        arrowprops=dict(arrowstyle="-|>", color=color, lw=1.5, mutation_scale=9,
                        shrinkA=0, shrinkB=0),
        zorder=6,
    )


def _plot_panel(ax, box_measure, bous_measure, cx_measure, plasim_col,
                ylabel, panel_title, df_box, df_plasim, df_boussinesq, df_climberx,
                xlabel=False, mark_line_ends=False, box_drop_n=1, box_arrow=False):
    """Plot a single resilience-measure panel onto *ax*.

    When *mark_line_ends* is True, a hollow circle is drawn at the final point of
    the 3-box and Boussinesq lines. For the 3-box model it marks the endpoint
    approaching the AMOC-on bifurcation; for Boussinesq it marks the last
    verified stable ON state, not a certified bifurcation.

    *box_drop_n* is the number of highest-CO2 box-model points omitted from that
    panel: 1 drops the collapsed/off-branch terminal (e.g. the 910 ppm point in
    AMOC strength), 0 shows every on-state point, and 2 truncates one step
    earlier (used for convergence time, whose near-bifurcation values explode).

    When *box_arrow* is True, an upward arrow in the 3-box colour is drawn from
    the last plotted box-model point to signal the omitted further increase."""

    box_tail = slice(None, -box_drop_n) if box_drop_n > 0 else slice(None)
    arrow_at = None

    # Remove the near-bifurcation box-model point(s) from every panel except (e, f),
    # before *box_drop_n* is applied, so the remaining truncation is unchanged.
    if df_box is not None and box_measure not in BOX_BASIN_ONLY_MEASURES:
        near = np.zeros(len(df_box), dtype=bool)
        for c in BOX_BASIN_ONLY_CO2:
            near |= np.abs(df_box["co2_ppm"].to_numpy(dtype=float) - c) <= BOX_CO2_MATCH_TOL
        df_box = df_box[~near]

    # ── Box model line ────────────────────────────────────────────────────
    if df_box is not None and box_measure is not None:
        sub = df_box[
            (df_box["measure"] == box_measure) &
            (df_box["attractor"] == "on")
        ].sort_values("co2_ppm").iloc[box_tail]

        if not sub.empty:
            ax.plot(
                sub["co2_ppm"].values,
                sub["value"].values,
                color=COL_ON,
                lw=1.5,
                label="3-box model",
                zorder=2,
            )
            if mark_line_ends:
                _mark_line_end(ax, sub, COL_ON)
            if box_arrow:
                arrow_at = (sub["co2_ppm"].values[-1], sub["value"].values[-1])
        else:
            # Try without attractor filter (e.g. amoc_strength_sv)
            sub_all = df_box[df_box["measure"] == box_measure].sort_values("co2_ppm").iloc[box_tail]
            if not sub_all.empty:
                ax.plot(
                    sub_all["co2_ppm"].values,
                    sub_all["value"].values,
                    color=COL_ON,
                    lw=1.5,
                    label="3-box model",
                    zorder=2,
                )
                if mark_line_ends:
                    _mark_line_end(ax, sub_all, COL_ON)
                if box_arrow:
                    arrow_at = (sub_all["co2_ppm"].values[-1], sub_all["value"].values[-1])

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
                lw=1.5,
                label="Boussinesq",
                zorder=3,
            )
            if mark_line_ends:
                _mark_line_end(ax, sub_b.dropna(subset=["co2_ppm", "value"]), COL_BOUS)

    # ── CLIMBER-X: two broken series ─────────────────────────────────────
    # The long-equilibrium AMOC state changes branch with CO2: it is the
    # reference "modern" state up to 330 ppm and the "strong" state from
    # 345 ppm onward (attractor column set by export_resilience_csv.py).
    # The two regimes describe different attractors, so they are drawn as two
    # separate series (circles vs stars) and are NOT connected across the gap.
    if df_climberx is not None and cx_measure is not None:
        sub_cx = df_climberx[df_climberx["measure"] == cx_measure]
        for state, marker, ms, lbl in (
            ("modern", "o", 4, "CLIMBER-X (modern)"),
            ("strong", "*", 7, "CLIMBER-X (strong)"),
        ):
            s = sub_cx[sub_cx["attractor"] == state].sort_values("co2_ppm")
            if s.empty:
                continue
            ax.plot(
                s["co2_ppm"].values, s["value"].values,
                color=COL_CLIMBERX, lw=1.5, marker=marker, markersize=ms,
                zorder=4, label=lbl,
            )

    if arrow_at is not None:
        _add_up_arrow(ax, *arrow_at, COL_ON)

    ax.set_ylabel(ylabel, fontsize=8)
    ax.tick_params(labelsize=7)
    if xlabel:
        ax.set_xlabel("CO\u2082 concentration (ppm)", fontsize=8)


def main() -> None:
    global BOUSSINESQ_CSV
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--umbrella", type=Path, default=UMBRELLA)
    parser.add_argument("--boussinesq-csv", type=Path, default=BOUSSINESQ_CSV)
    parser.add_argument("--output", type=Path, default=UMBRELLA / "plots/synthesis_figure.png")
    args = parser.parse_args()
    BOUSSINESQ_CSV = args.boussinesq_csv.resolve()
    apply_style()

    df_box       = load_box_model()
    df_plasim    = load_plasim()
    df_boussinesq = load_boussinesq()
    df_climberx  = load_climberx()

    if df_box is None and df_plasim is None and df_boussinesq is None and df_climberx is None:
        print("No data available for synthesis figure. Generate model data first.")
        sys.exit(1)

    # Tropical−Northern salinity gradient series for the box and PlaSim models
    # (Boussinesq and CLIMBER-X carry it in their resilience CSVs already).
    box_sal_df    = load_box_salinity_gradient()
    plasim_sal_df = load_plasim_salinity_gradient()

    # ── Build figure ─────────────────────────────────────────────────────────
    # Layout: 2 columns × 3 rows.  Top row holds the two model-state panels —
    # (a) AMOC strength and (b) the Tropical−Northern salinity gradient — and the
    # four resilience-measure panels (c)–(f) fill the two lower rows.  The three
    # panels in each column share an x-axis (CO2 concentration).
    from matplotlib.gridspec import GridSpec

    ncols = 2
    nrows_total = 3

    fig = plt.figure(figsize=(FIGURE_WIDTH, 5.4), constrained_layout=True)
    gs  = GridSpec(nrows_total, ncols, figure=fig,
                   height_ratios=[1.0, 1.0, 1.0])

    # ── Top row: AMOC strength (a) and salinity gradient (b) ─────────────────
    ax_amoc = fig.add_subplot(gs[0, 0])
    _plot_panel(
        ax_amoc,
        *AMOC_PANEL,
        df_box, df_plasim, df_boussinesq, df_climberx,
        xlabel=False,
        mark_line_ends=True,
    )
    ax_amoc.margins(y=0.08)  # Leave room for the offset endpoint circles.
    add_panel_label(ax_amoc, "(a)", x=0.99, ha="right")

    ax_sal = fig.add_subplot(gs[0, 1])
    _plot_panel(
        ax_sal,
        "salinity_gradient", "salinity_gradient", "salinity_gradient", "salinity_gradient",
        "Salinity gradient (psu)", "Salinity gradient",
        box_sal_df, plasim_sal_df, df_boussinesq, df_climberx,
        xlabel=False,
        box_drop_n=6,  # end at 826 ppm; the steep near-bifurcation rise is signalled by an arrow
        box_arrow=True,
    )
    # Keep the label clear of the box-model rise at the narrower right edge.
    add_panel_label(ax_sal, "(b)", x=0.02, ha="left")

    # Hide the top-row tick labels; the shared x-axis is labelled on the bottom row.
    col_top = {0: ax_amoc, 1: ax_sal}
    for ax in col_top.values():
        plt.setp(ax.get_xticklabels(), visible=False)

    # ── Rows 1-2: the four resilience-measure panels, shared x per column ─────
    # How many trailing box-model points to omit per measure: local resilience,
    # basin volume and critical shock show every on-state point (0) -- through
    # 903 ppm for local resilience, and through 906.64 ppm for basin volume and
    # critical shock (see BOX_BASIN_ONLY_CO2); convergence time is truncated one step earlier
    # (2) because its near-bifurcation values explode and compress the panel.
    BOX_DROP_N = {
        "local_resilience": 0,
        "basin_stability": 0,
        "minimal_critical_shock_magnitude": 0,
        "mean_convergence_time": 2,
    }
    # Measures whose truncated box-model curve gets an upward arrow signalling
    # the omitted near-bifurcation increase.
    BOX_ARROW = {"mean_convergence_time"}
    panel_labels = ["(c)", "(d)", "(e)", "(f)"]
    for panel_idx, (box_measure, bous_measure, cx_measure, plasim_col, ylabel, panel_title) in enumerate(PANELS):
        row = 1 + panel_idx // ncols
        col = panel_idx % ncols
        is_bottom_row = (row == nrows_total - 1)

        ax = fig.add_subplot(gs[row, col], sharex=col_top[col])

        _plot_panel(
            ax,
            box_measure, bous_measure, cx_measure, plasim_col,
            ylabel, panel_title,
            df_box, df_plasim, df_boussinesq, df_climberx,
            xlabel=is_bottom_row,
            box_drop_n=BOX_DROP_N.get(box_measure, 1),
            box_arrow=(box_measure in BOX_ARROW),
        )
        if not is_bottom_row:
            plt.setp(ax.get_xticklabels(), visible=False)

        if panel_labels[panel_idx] == "(d)":
            add_panel_label(ax, panel_labels[panel_idx], x=0.02, ha="left")
        else:
            add_panel_label(ax, panel_labels[panel_idx], x=0.99, ha="right")

    # Shared legend – placed below figure
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color=COL_ON,       lw=1.5, marker="",  markersize=6, label="3-box model"),
        Line2D([0], [0], color=COL_BOUS,     lw=1.5, marker="",  markersize=6, label="Boussinesq"),
        Line2D([0], [0], color=COL_CLIMBERX, lw=1.5, marker="o", markersize=5, label="CLIMBER-X (modern)"),
        Line2D([0], [0], color=COL_CLIMBERX, lw=1.5, marker="*", markersize=7, label="CLIMBER-X (strong)"),
        Line2D([0], [0], color=COL_PLASIM,   lw=0,   marker="^", markersize=6, label="PlaSim"),
    ]
    fig.legend(
        handles=legend_elements,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.07),
        ncol=5,
        fontsize=7,
        framealpha=0.8,
    )

    out_path = args.output.resolve().with_suffix(".png")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    sources = [BOX_CSV, PLASIM_CSV, BOUSSINESQ_CSV, CLIMBERX_CSV,
               UMBRELLA / "amoc_plot_style.py"]
    sources += sorted(BOX_PAPER_DIR.glob("attractors_*ppm.csv"))
    sources += sorted(PLASIM_PAPER_DIR.glob("state_means_*ppm.csv"))
    conversion_provenance = BOUSSINESQ_CSV.with_suffix(".json")
    conversion_metadata = {}
    if conversion_provenance.is_file():
        sources.append(conversion_provenance)
        conversion_metadata = json.loads(conversion_provenance.read_text())
    provenance = {
        "source_files_sha256": {str(path): hashlib.sha256(path.read_bytes()).hexdigest()
                                for path in sources if path.is_file()},
        "plot_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "layout": "Original 2-column, 3-row layout, sizes and legend; CO2 axis follows updated data; labels b and d moved left to avoid curves",
        "boussinesq_units": "Converted paper values share the existing axes; native values retained in input CSV",
        "boussinesq_conversion_provenance": str(conversion_provenance),
        "boussinesq_transport_mapping": conversion_metadata.get("transport_mapping", {}),
        "amoc_endpoint_circles": {
            "3-box": "Last plotted ON point approaching the AMOC-on bifurcation",
            "Boussinesq": "Last verified stable ON state; not a certified bifurcation",
        },
        "boussinesq_co2_coordinate": conversion_metadata.get("co2_display_relation", "See input CSV provenance"),
        "boussinesq_co2_is_calibrated": conversion_metadata.get("co2_is_calibrated", False),
        "boussinesq_unavailable_on_gammas": sorted(df_boussinesq.loc[
            df_boussinesq["status"] == "on_reference_unavailable", "gamma"].unique().tolist())
            if df_boussinesq is not None else [],
        "other_model_values": "unchanged",
    }
    out_path.with_suffix(".json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"Figure saved: {out_path}")
    plt.close(fig)


if __name__ == "__main__":
    main()
