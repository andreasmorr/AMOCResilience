"""
plotting_perturbations.py  –  Multi-model perturbation & readout overview figure.

4-column × 2-row layout:

  Col 1  CLIMBER-X non-tapered Wood boxes
         (a) Globe  — box footprints on Atlantic
         (e) Section — deep box depth structure (lat × depth at 26°W)

  Col 2  Boussinesq context
         (b) Globe  — same Wood-box footprints for reference
         (f) 2D Boussinesq model domain with North Atlantic & Tropical box masks

  Col 3  CLIMBER-X tapered perturbation
         (c) Globe  — taper-weight field on 5° ocean grid
         (g) Section — tapered shallow-box cross-section, zoomed 0–300 m

  Col 4  PlaSim / EOF weighting  (data placeholders; to be updated)
         (d) Globe  — EOF-1 and EOF-2 spatial patterns
         (h) Section — meridional EOF cross-section

Requires: matplotlib, numpy, scipy, cartopy
Output:   plots/perturbations_overview.pdf
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from scipy.ndimage import distance_transform_edt

# ---------------------------------------------------------------------------
# Cartopy (optional — fall back to flat map if absent)
# ---------------------------------------------------------------------------
try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    HAS_CARTOPY = True
except ImportError:
    HAS_CARTOPY = False
    print("cartopy not found — globe panels will use a flat Plate Carrée map.")

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent       # AMOCResilience/
PLOTS_DIR  = SCRIPT_DIR / "plots"

sys.path.insert(0, str(SCRIPT_DIR))
from amoc_plot_style import apply_style, add_panel_label, savefig_pdf

# ---------------------------------------------------------------------------
# CLIMBER-X 5° ocean grid  (23 levels, 36 lat × 72 lon cells)
# ---------------------------------------------------------------------------
LAT = np.arange(-87.5, 90.0, 5.0)          # (36,) cell centres, °N
LON = np.arange(-177.5, 180.0, 5.0)        # (72,) cell centres, °E
ZRO = np.array([5, 25, 75, 150, 250, 350, 450, 575, 700, 875, 1100, 1350,
                1650, 2000, 2400, 2800, 3200, 3600, 4000, 4400, 4800, 5200,
                5700], dtype=float)         # (23,) approximate depth-level centres (m)

# Grid edges for pcolormesh
LAT_EDGES = np.concatenate([[LAT[0] - 2.5], 0.5 * (LAT[:-1] + LAT[1:]), [LAT[-1] + 2.5]])
LON_EDGES = np.concatenate([[LON[0] - 2.5], 0.5 * (LON[:-1] + LON[1:]), [LON[-1] + 2.5]])
ZRO_EDGES = np.concatenate([[0.0], ZRO])    # top edge of shallowest level = 0 m

# ---------------------------------------------------------------------------
# Land geometry — loaded once at 10 m resolution for smooth clipping and
# at 50 m for the coarse grid-cell ocean mask used by the tapered panel.
# ---------------------------------------------------------------------------
def _load_land(resolution="10m"):
    import cartopy.io.shapereader as shapereader
    import shapely.ops
    path   = shapereader.natural_earth(
        resolution=resolution, category="physical", name="land")
    geoms  = list(shapereader.Reader(path).geometries())
    return shapely.ops.unary_union(geoms)

print("Loading land polygons (10 m) for smooth box clipping …")
try:
    LAND_UNION = _load_land("10m")   # high-res for shapely clipping
except Exception as _exc:
    print(f"  [warning] 10 m land unavailable ({_exc}); trying 50 m")
    LAND_UNION = _load_land("50m")

# Coarse grid-cell ocean mask (True = ocean) — used only for the tapered panel
def _build_ocean_mask(land_union):
    try:
        from shapely.vectorized import contains as _contains
        lon2d, lat2d = np.meshgrid(LON, LAT)
        is_land = _contains(land_union, lon2d.ravel(), lat2d.ravel())
        return (~is_land).reshape(lon2d.shape)
    except Exception as exc:
        print(f"  [warning] ocean mask unavailable ({exc})")
        return np.ones((len(LAT), len(LON)), dtype=bool)

OCEAN_MASK = _build_ocean_mask(LAND_UNION)   # (36, 72) bool

# Fine-resolution (0.5°) grid for smooth globe panels — land-masked via LAND_UNION
_FINE_RES = 0.5
_FINE_LON = np.arange(-179.75, 180.0, _FINE_RES)
_FINE_LAT = np.arange( -89.75,  90.0, _FINE_RES)
_FINE_LON_EDGES = np.concatenate([
    [_FINE_LON[0] - _FINE_RES/2],
    0.5*(_FINE_LON[:-1] + _FINE_LON[1:]),
    [_FINE_LON[-1] + _FINE_RES/2],
])
_FINE_LAT_EDGES = np.concatenate([
    [_FINE_LAT[0] - _FINE_RES/2],
    0.5*(_FINE_LAT[:-1] + _FINE_LAT[1:]),
    [_FINE_LAT[-1] + _FINE_RES/2],
])

def _build_fine_ocean_mask():
    lon2d, lat2d = np.meshgrid(_FINE_LON, _FINE_LAT)
    lons, lats = lon2d.ravel(), lat2d.ravel()
    try:
        from shapely.vectorized import contains as _vc
        is_land = _vc(LAND_UNION, lons, lats)
    except Exception:
        try:
            from shapely import contains_xy
            is_land = contains_xy(LAND_UNION, lons, lats)
        except Exception as exc:
            print(f"  [warning] fine ocean mask unavailable ({exc})")
            return np.ones((len(_FINE_LAT), len(_FINE_LON)), dtype=bool)
    return (~is_land).reshape(lon2d.shape)

print("Building fine ocean mask (0.5°) for globe panels …")
FINE_OCEAN_MASK = _build_fine_ocean_mask()   # (360, 720) bool

# ---------------------------------------------------------------------------
# Wood-box definitions  (full/deep variants and shallow variants)
# ---------------------------------------------------------------------------
BOX_COLOR_NA    = "#2166ac"
BOX_COLOR_TROP  = "#d6604d"
BOX_COLOR_SOUTH = "#4dac26"

BOXES = {
    "NA":    dict(lat_min=37.5,  lat_max=90.0,  lon_min=-72.5, lon_max=20.0,
                  depth_max=None,  color=BOX_COLOR_NA,    label="WOOD_NA"),
    "Trop":  dict(lat_min=-47.5, lat_max=32.5,  lon_min=-72.5, lon_max=20.0,
                  depth_max=875.0, color=BOX_COLOR_TROP,  label="WOOD_TROP"),
    "South": dict(lat_min=-90.0, lat_max=-52.5, lon_min=-180.0, lon_max=180.0,
                  depth_max=450.0, color=BOX_COLOR_SOUTH, label="WOOD_SOUTH"),
}

# Globe panels a and b: all three boxes with gaps closed.
#   Tropical: lat_max extended 32.5°N → 37.5°N (closes gap with NA)
#   South:    lat_max extended -52.5°S → -47.5°S (closes gap with Tropical)
BOXES_GLOBE = {
    "NA":    BOXES["NA"],
    "Trop":  {**BOXES["Trop"], "lat_max": 37.5},
    "South": {**BOXES["South"], "lat_max": -47.5},
}

# Atlantic-only subset (kept for reference)
BOXES_ATLANTIC = {
    "NA":   BOXES["NA"],
    "Trop": BOXES["Trop"],
}

BOXES_SHALLOW = {
    "NA":    dict(lat_min=37.5,  lat_max=90.0,  lon_min=-72.5, lon_max=20.0,
                  depth_max=150.0, color=BOX_COLOR_NA,    label="WOOD_NA_SHALLOW"),
    "Trop":  dict(lat_min=-47.5, lat_max=32.5,  lon_min=-72.5, lon_max=20.0,
                  depth_max=150.0, color=BOX_COLOR_TROP,  label="WOOD_TROP_SHALLOW"),
    "South": dict(lat_min=-90.0, lat_max=-52.5, lon_min=-180.0, lon_max=180.0,
                  depth_max=150.0, color=BOX_COLOR_SOUTH, label="WOOD_SOUTH_SHALLOW"),
}

TAPER_MARGIN_CELLS  = 2
TAPER_MARGIN_LAYERS = 2

# ---------------------------------------------------------------------------
# Mask / taper helpers  (mirrors AMOCClimberX/src/restart.py)
# ---------------------------------------------------------------------------

def horiz_mask(lat_min, lat_max, lon_min, lon_max):
    lat_ok   = (LAT >= lat_min) & (LAT <= lat_max)
    lo_min   = ((lon_min + 180) % 360) - 180
    lo_max   = ((lon_max + 180) % 360) - 180
    lon_norm = ((LON + 180) % 360) - 180
    if (lon_max - lon_min) >= 360:
        lon_ok = np.ones(len(LON), dtype=bool)
    elif lo_min <= lo_max:
        lon_ok = (lon_norm >= lo_min) & (lon_norm <= lo_max)
    else:
        lon_ok = (lon_norm >= lo_min) | (lon_norm <= lo_max)
    return lat_ok[:, None] & lon_ok[None, :]   # (lat, lon)


def horiz_taper(hmask):
    dist = distance_transform_edt(hmask)
    m = float(TAPER_MARGIN_CELLS)
    w = np.where(dist >= m, 1.0, 0.5 * (1.0 - np.cos(np.pi * dist / m)))
    w[~hmask] = 0.0
    return w


def vert_mask(depth_max):
    if depth_max is None:
        return np.ones(len(ZRO), dtype=bool)
    return ZRO <= depth_max


def vert_taper(depth_max):
    inside     = vert_mask(depth_max)
    n          = int(inside.sum())
    if n == 0:
        return np.zeros(len(ZRO), dtype=float)
    inside_idx = np.where(inside)[0]
    w          = np.zeros(len(ZRO), dtype=float)
    m          = TAPER_MARGIN_LAYERS
    for rank, orig in enumerate(inside_idx):
        rank_from_bottom = (n - 1) - rank
        if rank_from_bottom < m:
            frac   = (rank_from_bottom + 1) / (m + 1)
            w[orig] = 0.5 * (1.0 - np.cos(np.pi * frac))
        else:
            w[orig] = 1.0
    return w


def make_3d_mask(box, taper=False):
    hm   = horiz_mask(box["lat_min"], box["lat_max"], box["lon_min"], box["lon_max"])
    vm   = vert_mask(box["depth_max"])
    full = vm[:, None, None] & hm[None, :, :]
    if not taper:
        return full.astype(float)
    hw  = horiz_taper(hm)
    vw  = vert_taper(box["depth_max"])
    w3  = vw[:, None, None] * hw[None, :, :]
    w3[~full] = 0.0
    return w3


# Cross-section longitude: midpoint of [-72.5, 20°E] → 26.25°W
SECTION_IDX = int(np.argmin(np.abs(LON - (-26.25))))
SECTION_LON = float(LON[SECTION_IDX])   # actual grid longitude used

# ---------------------------------------------------------------------------
# Boussinesq 2D model geometry  (M=20 lat cells, N=40 depth cells)
# Coordinate system: x ∈ [0, 5]  →  lat = x/5 * 180 − 90 (°N)
#                    z ∈ [0, 1]  →  normalised depth (0 = surface, 1 = bottom)
# ---------------------------------------------------------------------------
BOUS_M, BOUS_N   = 20, 40
BOUS_SMOOTH_X    = 0.25                         # tanh half-width (model units)
BOUS_XX          = np.linspace(0.0, 5.0, BOUS_M + 1)   # (21,) latitude coords
BOUS_ZZ          = np.linspace(0.0, 1.0, BOUS_N + 1)   # (41,) depth coords

# Surface-intensified weighting h(z) ~ exp(-5z); same role as Boussinesq h(z)
_hz   = np.exp(-5.0 * BOUS_ZZ)
_hz  /= _hz.max()


def bous_lat_mask(x_lo, x_hi):
    """Smooth tanh meridional mask on BOUS_XX."""
    return 0.5 * (np.tanh((BOUS_XX - x_lo) / BOUS_SMOOTH_X)
                  - np.tanh((BOUS_XX - x_hi) / BOUS_SMOOTH_X))


def bous_box_mask(x_lo, x_hi):
    """2D box mask (lat × depth), shape (M+1, N+1), peak-normalised to 1."""
    lm = bous_lat_mask(x_lo, x_hi)
    mask = np.outer(lm, _hz)
    if mask.max() > 0:
        mask /= mask.max()
    return mask


# Geographic bounds of Boussinesq boxes (matching Boussinesq_box.py __main__)
BOUS_NORTH_LO = 130.0 / 36.0   # ≈ 3.611  →  40°N
BOUS_NORTH_HI = 5.0             #             90°N
BOUS_TROP_LO  = 40.0  / 36.0   # ≈ 1.111  → −50°S
BOUS_TROP_HI  = 130.0 / 36.0   # ≈ 3.611  →  40°N
BOUS_SOUTH_LO = -1.0            #             (below domain → clamped to 90°S)
BOUS_SOUTH_HI = 40.0  / 36.0   # ≈ 1.111  → −50°S

BOUS_NORTH_MASK = bous_box_mask(BOUS_NORTH_LO, BOUS_NORTH_HI)   # (21, 41)
BOUS_TROP_MASK  = bous_box_mask(BOUS_TROP_LO,  BOUS_TROP_HI)    # (21, 41)
BOUS_SOUTH_MASK = bous_box_mask(BOUS_SOUTH_LO, BOUS_SOUTH_HI)   # (21, 41)

# Convert Boussinesq x to degrees latitude for axis labels
BOUS_LAT_DEG = BOUS_XX / 5.0 * 180.0 - 90.0    # (21,)


# ---------------------------------------------------------------------------
# Globe-panel setup helpers
# ---------------------------------------------------------------------------

_GLOBE_PROJ    = ccrs.Orthographic(-30, 20) if HAS_CARTOPY else None
_PC_TRANSFORM  = ccrs.PlateCarree()         if HAS_CARTOPY else None

_OCEAN_COLOR   = "#cce5f5"
_LAND_COLOR    = "#e0e0e0"
_GRID_KW       = dict(linewidth=0.3, color="gray", alpha=0.4, linestyle="--")


def setup_globe(ax):
    """Add ocean / land / coastline features to a GeoAxes or flat fallback."""
    if HAS_CARTOPY:
        ax.add_feature(cfeature.OCEAN.with_scale("50m"),
                       facecolor=_OCEAN_COLOR, zorder=0)
        ax.add_feature(cfeature.LAND.with_scale("50m"),
                       facecolor=_LAND_COLOR, zorder=1)
        ax.add_feature(cfeature.COASTLINE.with_scale("50m"),
                       linewidth=0.4, edgecolor="#666666", zorder=2)
        ax.gridlines(**_GRID_KW)
    else:
        ax.set_facecolor(_OCEAN_COLOR)
        ax.set_xlim(-85, 35)
        ax.set_ylim(-75, 90)
        ax.axhline(0, color="gray", linewidth=0.4, linestyle="--")
        ax.set_xlabel("Longitude (°E)", fontsize=7)
        ax.set_ylabel("Latitude (°N)", fontsize=7)


def _pcolor_globe(ax, w2, color_hex):
    """Overlay a (lat × lon) weight field on the globe as a colour fill.
    Land cells are zeroed out so only ocean cells receive colour."""
    w2 = w2 * OCEAN_MASK          # mask land cells to weight = 0
    rgba = np.array(mcolors.to_rgba(color_hex))
    cmap = LinearSegmentedColormap.from_list(
        "", [(rgba[0], rgba[1], rgba[2], 0.0), (*rgba[:3], 0.85)], N=256
    )
    kw = dict(cmap=cmap, vmin=0, vmax=1, shading="flat", zorder=3)
    if HAS_CARTOPY:
        ax.pcolormesh(LON_EDGES, LAT_EDGES, w2, transform=_PC_TRANSFORM, **kw)
    else:
        ax.pcolormesh(LON_EDGES, LAT_EDGES, w2, **kw)


def draw_boxes_smooth(ax, box_dict):
    """High-resolution fill: subtract land from each box polygon and draw
    the ocean-only geometry.  Requires cartopy + shapely."""
    from shapely.geometry import box as shp_box, MultiPolygon, Polygon
    from shapely.ops import unary_union

    for key, box in box_dict.items():
        rect = shp_box(box["lon_min"], box["lat_min"],
                       box["lon_max"], box["lat_max"])
        raw = rect.difference(LAND_UNION)
        if raw.is_empty:
            continue

        # Collect individual Polygon parts and re-orient them so that the
        # exterior ring is CCW (positive area) and holes are CW.  Without
        # this, cartopy sometimes fills the complement (the whole globe
        # minus the polygon) instead of the polygon itself.
        parts = list(raw.geoms) if raw.geom_type == "MultiPolygon" else [raw]
        fixed = []
        for poly in parts:
            if not isinstance(poly, Polygon) or poly.is_empty:
                continue
            # Re-create with correct ring orientation
            ext = list(poly.exterior.coords)
            holes = [list(r.coords) for r in poly.interiors]
            fixed.append(Polygon(ext, holes))

        if not fixed:
            continue

        rgba = mcolors.to_rgba(box["color"])
        ax.add_geometries(
            fixed,
            crs=_PC_TRANSFORM,
            facecolor=(*rgba[:3], 0.55),
            edgecolor="none",
            zorder=3,
        )


def draw_boxes_fine_grid(ax, box_dict):
    """Plot box regions using the 0.5° FINE_OCEAN_MASK — smooth and land-excluded.
    Avoids cartopy polygon-fill artefacts that plague add_geometries on orthographic."""
    for key, box in box_dict.items():
        lat_ok = (_FINE_LAT >= box["lat_min"]) & (_FINE_LAT <= box["lat_max"])
        lo_min = ((box["lon_min"] + 180) % 360) - 180
        lo_max = ((box["lon_max"] + 180) % 360) - 180
        lon_norm = ((_FINE_LON + 180) % 360) - 180
        if (box["lon_max"] - box["lon_min"]) >= 360:
            lon_ok = np.ones(len(_FINE_LON), dtype=bool)
        elif lo_min <= lo_max:
            lon_ok = (lon_norm >= lo_min) & (lon_norm <= lo_max)
        else:
            lon_ok = (lon_norm >= lo_min) | (lon_norm <= lo_max)
        w2 = (lat_ok[:, None] & lon_ok[None, :] & FINE_OCEAN_MASK).astype(float)
        rgba = np.array(mcolors.to_rgba(box["color"]))
        cmap = LinearSegmentedColormap.from_list(
            "", [(rgba[0], rgba[1], rgba[2], 0.0), (*rgba[:3], 0.65)], N=2
        )
        kw = dict(cmap=cmap, vmin=0, vmax=1, shading="flat", zorder=3)
        if HAS_CARTOPY:
            ax.pcolormesh(_FINE_LON_EDGES, _FINE_LAT_EDGES, w2,
                         transform=_PC_TRANSFORM, **kw)
        else:
            ax.pcolormesh(_FINE_LON_EDGES, _FINE_LAT_EDGES, w2, **kw)


def draw_boxes_on_globe(ax, box_dict, taper=False):
    """Plot depth-averaged box weights on the globe as a pcolormesh.
    Used for the tapered panel where per-cell weights must be shown."""
    for key, box in box_dict.items():
        w3 = make_3d_mask(box, taper=taper)
        vm = vert_mask(box["depth_max"])
        n_lev = int(vm.sum())
        if n_lev == 0:
            continue
        w2 = w3[vm, :, :].mean(axis=0)   # depth-averaged weight, (lat, lon)
        _pcolor_globe(ax, w2, box["color"])


def box_legend(ax, box_dict, fontsize=7, loc="lower left"):
    handles = [
        mpatches.Patch(color=b["color"], alpha=0.7, label=b["label"])
        for b in box_dict.values()
    ]
    ax.legend(handles=handles, fontsize=fontsize, loc=loc,
              framealpha=0.85, edgecolor="#cccccc")


# ---------------------------------------------------------------------------
# Cross-section helper: meridional section (lat × depth) using pcolormesh
# ---------------------------------------------------------------------------

def draw_section(ax, box_dict, taper=False, depth_max_plot=None, label_depths=True):
    """
    Meridional cross-section at SECTION_LON.
    depth_max_plot: clip y-axis to this depth (m); None = full depth.
    """
    from matplotlib.colors import LinearSegmentedColormap as LSC
    for key, box in box_dict.items():
        w3 = make_3d_mask(box, taper=taper)   # (lev, lat, lon)
        w2 = w3[:, :, SECTION_IDX]            # (lev, lat)
        if w2.max() == 0:
            continue
        rgba = np.array(mcolors.to_rgba(box["color"]))
        cmap = LSC.from_list(key, [(1, 1, 1, 0), rgba], N=256)
        ax.pcolormesh(LAT_EDGES, ZRO_EDGES, w2,
                      cmap=cmap, vmin=0, vmax=1, shading="flat", zorder=2)

    # Depth boundary lines
    if label_depths:
        for key, box in box_dict.items():
            dmax = box.get("depth_max")
            if dmax is not None:
                ax.axhline(dmax, color=box["color"], linewidth=1.0,
                           linestyle="--", alpha=0.7, zorder=5)
                ax.text(87.0, dmax + 25, f"{int(dmax)} m",
                        fontsize=6.5, color=box["color"], ha="right", va="top")

    # Latitude boundary lines
    for key, box in box_dict.items():
        for lat_edge in [box["lat_min"], box["lat_max"]]:
            ax.axvline(lat_edge, color=box["color"], linewidth=0.8,
                       linestyle=":", alpha=0.5, zorder=5)

    depth_lim = depth_max_plot if depth_max_plot is not None else float(ZRO[-1] + 300)
    ax.set_xlim(-90, 90)
    ax.set_ylim(depth_lim, 0)
    ax.set_facecolor("#f0f0f0")
    ax.set_xlabel("Latitude", fontsize=8)
    ax.axvline(0, color="gray", linewidth=0.4, linestyle="--", alpha=0.5)

    ticks = np.arange(-90, 91, 30)
    ax.set_xticks(ticks)
    ax.set_xticklabels(
        [f"{abs(t)}°{'N' if t >= 0 else 'S'}" for t in ticks], fontsize=6.5
    )


# ---------------------------------------------------------------------------
# Boussinesq 2D panel
# ---------------------------------------------------------------------------

def draw_boussinesq_panel(ax):
    """
    2D contourf plot of the Boussinesq box masks in (latitude, normalised depth).
    Southern box: Greens.  Tropical box: Oranges.  North Atlantic box: Blues.
    """
    XX, ZZ = np.meshgrid(BOUS_LAT_DEG, BOUS_ZZ)   # both (41, 21)
    mask_n = BOUS_NORTH_MASK.T   # (41, 21) so shape matches meshgrid
    mask_t = BOUS_TROP_MASK.T
    mask_s = BOUS_SOUTH_MASK.T

    # Start levels above 0 so that zero-value cells are left unfilled (transparent).
    # Skip the pale end of each colormap so filled regions are visibly dark.
    levels = np.linspace(0.05, 1.0, 20)

    cmap_s = LinearSegmentedColormap.from_list("south_dark",
        plt.cm.Greens(np.linspace(0.50, 1.0, 256)))
    cmap_t = LinearSegmentedColormap.from_list("trop_dark",
        plt.cm.Oranges(np.linspace(0.38, 1.0, 256)))
    cmap_n = LinearSegmentedColormap.from_list("north_dark",
        plt.cm.Blues(np.linspace(0.30, 1.0, 256)))

    # Draw back-to-front so North Atlantic is on top
    ax.contourf(XX, ZZ, mask_s, levels=levels, cmap=cmap_s, alpha=0.95, zorder=2,
                extend="neither")
    ax.contourf(XX, ZZ, mask_t, levels=levels, cmap=cmap_t, alpha=0.95, zorder=3,
                extend="neither")
    ax.contourf(XX, ZZ, mask_n, levels=levels, cmap=cmap_n, alpha=0.95, zorder=4,
                extend="neither")

    # Box boundary lines (latitude only — vertical lines)
    for x_edge, col in [
        (BOUS_SOUTH_HI / 5.0 * 180 - 90, BOX_COLOR_SOUTH),
        (BOUS_TROP_LO  / 5.0 * 180 - 90, BOX_COLOR_TROP),
        (BOUS_TROP_HI  / 5.0 * 180 - 90, BOX_COLOR_TROP),
        (BOUS_NORTH_LO / 5.0 * 180 - 90, BOX_COLOR_NA),
        (BOUS_NORTH_HI / 5.0 * 180 - 90, BOX_COLOR_NA),
    ]:
        ax.axvline(x_edge, color=col, linewidth=1.0, linestyle=":", alpha=0.7, zorder=6)

    ax.set_facecolor("#f0f0f0")
    ax.set_xlim(-90, 90)
    ax.set_ylim(1.0, 0.0)   # depth increases downward
    ax.set_xlabel("Latitude", fontsize=8)
    ax.set_ylabel("Normalised depth", fontsize=8)

    ticks = np.arange(-90, 91, 30)
    ax.set_xticks(ticks)
    ax.set_xticklabels(
        [f"{abs(t)}°{'N' if t >= 0 else 'S'}" for t in ticks], fontsize=6.5
    )
    ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels(["0", "0.25", "0.5", "0.75", "1.0"], fontsize=6.5)
    ax.axvline(0, color="gray", linewidth=0.4, linestyle="--", alpha=0.5)

    pass  # legend handled by figure-level legend


# ---------------------------------------------------------------------------
# EOF placeholder panel
# ---------------------------------------------------------------------------

def draw_eof_placeholder(ax, title_text="", globe=True):
    ax.set_facecolor("#f5f5f5")
    if not globe:
        ax.set_xlim(-90, 90)
        ax.set_ylim(5000, 0)
    style_kw = dict(ha="center", va="center", fontsize=9, color="#888888")
    ax.text(0.5, 0.5, title_text, transform=ax.transAxes, **style_kw)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_edgecolor("#cccccc")


# ---------------------------------------------------------------------------
# Main figure assembly
# ---------------------------------------------------------------------------

def main():
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    apply_style()

    # ── Build figure layout ────────────────────────────────────────────────
    fig = plt.figure(figsize=(17, 9))
    gs  = gridspec.GridSpec(
        2, 4,
        figure=fig,
        height_ratios=[1.35, 1.0],
        hspace=0.40,
        wspace=0.28,
        left=0.06, right=0.97,
        top=0.93,  bottom=0.09,
    )

    # Top row: globe panels (GeoAxes if cartopy, otherwise regular)
    if HAS_CARTOPY:
        axes_top = [fig.add_subplot(gs[0, j], projection=_GLOBE_PROJ)
                    for j in range(4)]
    else:
        axes_top = [fig.add_subplot(gs[0, j]) for j in range(4)]

    # Bottom row: always regular axes
    axes_bot = [fig.add_subplot(gs[1, j]) for j in range(4)]

    panel_labels_top = ["(a)", "(b)", "(c)", "(d)"]
    panel_labels_bot = ["(e)", "(f)", "(g)", "(h)"]

    # ── Column 1: Box model regions ───────────────────────────────────────
    ax = axes_top[0]
    setup_globe(ax)
    draw_boxes_fine_grid(ax, BOXES_GLOBE)
    ax.set_title("Box model regions (Wood et al.)", fontsize=8, fontweight="bold")
    add_panel_label(ax, panel_labels_top[0])

    ax = axes_bot[0]
    draw_section(ax, BOXES, taper=False, depth_max_plot=None, label_depths=True)
    ax.set_ylabel("Depth (m)", fontsize=8)
    ax.set_title("Meridional section (28°W)", fontsize=8)
    add_panel_label(ax, panel_labels_bot[0])

    # ── Column 2: Boussinesq context ──────────────────────────────────────
    ax = axes_top[1]
    setup_globe(ax)
    draw_boxes_fine_grid(ax, BOXES_GLOBE)   # identical to panel a
    ax.set_title("Regions corresponding to Boussinesq model boxes", fontsize=8,
                 fontweight="bold")
    add_panel_label(ax, panel_labels_top[1])

    ax = axes_bot[1]
    draw_boussinesq_panel(ax)
    ax.set_title("Boussinesq 2D box weighting", fontsize=8)
    add_panel_label(ax, panel_labels_bot[1])

    # ── Column 3: CLIMBER-X ───────────────────────────────────────────────
    ax = axes_top[2]
    setup_globe(ax)
    draw_boxes_on_globe(ax, BOXES_SHALLOW, taper=True)
    ax.set_title("CLIMBER-X model boxes", fontsize=8, fontweight="bold")
    add_panel_label(ax, panel_labels_top[2])

    ax = axes_bot[2]
    draw_section(ax, BOXES_SHALLOW, taper=True, depth_max_plot=320, label_depths=True)
    ax.set_ylabel("Depth (m)", fontsize=8)
    ax.set_title("Meridional section (28°W)", fontsize=8)
    add_panel_label(ax, panel_labels_bot[2])

    # ── Column 4: PlaSim ──────────────────────────────────────────────────
    ax = axes_top[3]
    if HAS_CARTOPY:
        setup_globe(ax)
        draw_eof_placeholder(ax, "EOF weighting\n(data to be provided)", globe=True)
    else:
        draw_eof_placeholder(ax, "EOF weighting\n(data to be provided)", globe=False)
    ax.set_title("PlaSim model boxes", fontsize=8, fontweight="bold")
    add_panel_label(ax, panel_labels_top[3])

    ax = axes_bot[3]
    draw_eof_placeholder(ax, "EOF cross-section\n(data to be provided)", globe=False)
    ax.set_title("Meridional EOF section", fontsize=8)
    ax.set_xlabel("Latitude", fontsize=8)
    ax.set_ylabel("Depth (m)", fontsize=8)
    add_panel_label(ax, panel_labels_bot[3])

    # ── Central legend between the two rows ───────────────────────────────
    legend_handles = [
        mpatches.Patch(color=BOX_COLOR_NA,    alpha=0.75, label="North Atlantic Box"),
        mpatches.Patch(color=BOX_COLOR_TROP,  alpha=0.75, label="Tropical Atlantic Box"),
        mpatches.Patch(color=BOX_COLOR_SOUTH, alpha=0.75, label="Southern Ocean Box"),
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.485),   # sits in the hspace gap between the rows
        ncol=3,
        fontsize=8,
        framealpha=0.9,
        edgecolor="#cccccc",
    )

    # ── Save as PNG ───────────────────────────────────────────────────────
    out_path = PLOTS_DIR / "perturbations_overview.png"
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    print(f"Figure saved: {out_path}")
    plt.close(fig)


if __name__ == "__main__":
    main()
