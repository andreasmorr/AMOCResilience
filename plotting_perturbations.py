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
ZRO = np.array([5.25, 18.39, 39.41, 78.82, 131.36, 183.90, 236.45, 288.99,
                367.81, 472.89, 577.98, 735.61, 945.79, 1182.23, 1444.95,
                1707.67, 1970.39, 2364.47, 2889.90, 3415.34, 3940.78,
                4466.21, 4991.65], dtype=float)   # (23,) CLIMBER-X layer-centre depths (m)

# Grid edges for pcolormesh
LAT_EDGES = np.concatenate([[LAT[0] - 2.5], 0.5 * (LAT[:-1] + LAT[1:]), [LAT[-1] + 2.5]])
LON_EDGES = np.concatenate([[LON[0] - 2.5], 0.5 * (LON[:-1] + LON[1:]), [LON[-1] + 2.5]])
ZRO_EDGES = np.concatenate([[0.0], ZRO])    # top edge of shallowest level = 0 m

# ---------------------------------------------------------------------------
# Atlantic ocean basin mask (basin_mask_5x5.nc), aligned to the 5° grid.
# Used by the Atlantic-mask boxes box_na / box_trop (basin_mask == i_atlantic).
# ---------------------------------------------------------------------------
_I_ATLANTIC = 1   # CLIMBER-X climber_grid.f90: i_atlantic = 1

def _load_atlantic_mask():
    """Boolean (lat, lon) Atlantic-basin selection on the 5° LAT/LON grid."""
    path = Path(__file__).resolve().parent / "basin_mask_5x5.nc"
    if not path.exists():
        print(f"WARNING: {path.name} not found — Atlantic-mask boxes will be empty.")
        return np.zeros((len(LAT), len(LON)), dtype=bool)
    import xarray as xr
    with xr.open_dataset(path) as ds:
        bm  = ds["basin_mask"].values
        blat = ds["lat"].values
        blon = ds["lon"].values
    sel = bm == _I_ATLANTIC
    if sel.shape == (len(LAT), len(LON)) and np.allclose(blat, LAT) and np.allclose(blon, LON):
        return sel
    j = np.array([int(np.argmin(np.abs(blat - v))) for v in LAT])
    i = np.array([int(np.argmin(np.abs(((blon - v + 180) % 360) - 180))) for v in LON])
    return sel[np.ix_(j, i)]

ATL_MASK = _load_atlantic_mask()   # (lat, lon) bool

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
_FINE_RES = 0.1
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


def _build_fine_atlantic_mask():
    """Atlantic-basin selection on the fine grid.

    The basin mask is only defined on the coarse 5° grid, so we upsample it to
    the fine grid by nearest-neighbour lookup and then intersect with the fine
    ocean mask.  Coastal edges (Americas, Africa/Europe) therefore follow the
    smooth high-resolution coastline, while the open-ocean basin boundaries
    (Arctic / Southern Ocean / Indian Ocean) inherit the 5° basin partition."""
    j = np.array([int(np.argmin(np.abs(LAT - v))) for v in _FINE_LAT])
    i = np.array([int(np.argmin(np.abs(((LON - v + 180) % 360) - 180)))
                  for v in _FINE_LON])
    atl_fine = ATL_MASK[np.ix_(j, i)]
    return atl_fine & FINE_OCEAN_MASK

FINE_ATL_MASK = _build_fine_atlantic_mask()   # (360, 720) bool

# Southern edge of the Atlantic basin mask (cell boundary, °N).  South of this
# the basin mask is empty (open South Atlantic merging into the Southern Ocean),
# so basin boxes that extend further south fall back to the longitude rectangle.
_atl_rows = np.where(ATL_MASK.any(axis=1))[0]
ATL_SOUTH_EDGE = float(LAT[_atl_rows.min()] - 2.5) if _atl_rows.size else -90.0
# Longitudinal span of the Atlantic at its southernmost row — used to clamp the
# south-of-basin extension so it stays in the South Atlantic (east of South
# America) and does not leak into the Pacific west of Chile.
if _atl_rows.size:
    _south_lons = LON[ATL_MASK[_atl_rows.min()]]
    ATL_SOUTH_LON_MIN = float(_south_lons.min() - 2.5)
    ATL_SOUTH_LON_MAX = float(_south_lons.max() + 2.5)
else:
    ATL_SOUTH_LON_MIN, ATL_SOUTH_LON_MAX = -180.0, 180.0

# ---------------------------------------------------------------------------
# Wood-box definitions  (full/deep variants and shallow variants)
# ---------------------------------------------------------------------------
BOX_COLOR_NA    = "#2166ac"
BOX_COLOR_TROP  = "#d6604d"
BOX_COLOR_SOUTH = "#4dac26"

BOXES = {
    "NA":    dict(lat_min=37.5,  lat_max=90.0,  lon_min=-70.5, lon_max=20.0,
                  depth_max=None,  color=BOX_COLOR_NA,    label="WOOD_NA"),
    "Trop":  dict(lat_min=-47.5, lat_max=32.5,  lon_min=-70.5, lon_max=20.0,
                  depth_max=875.0, color=BOX_COLOR_TROP,  label="WOOD_TROP"),
    "South": dict(lat_min=-90.0, lat_max=-52.5, lon_min=-180.0, lon_max=180.0,
                  depth_max=450.0, color=BOX_COLOR_SOUTH, label="WOOD_SOUTH"),
}

# Globe panels a and b: all three boxes with gaps closed.
#   Tropical: lat_max extended 32.5°N → 37.5°N (closes gap with NA)
#   South:    lat_max extended -52.5°S → -47.5°S (closes gap with Tropical)
# The NA and Tropical boxes follow the Atlantic basin mask so their E/W extent
# stretches to the coasts of the Americas and Africa/Europe (like CLIMBER-X);
# their longitude bounds are widened so the coastlines, not the rectangle,
# define the extent.  Latitude borders are unchanged.
BOXES_GLOBE = {
    "NA":    {**BOXES["NA"],   "lon_min": -100.0, "lon_max": 25.0,
              "basin": "atlantic"},
    "Trop":  {**BOXES["Trop"], "lat_max": 37.5, "lon_min": -100.0,
              "lon_max": 25.0, "basin": "atlantic"},
    "South": {**BOXES["South"], "lat_max": -47.5},
}

# Atlantic-only subset (kept for reference)
BOXES_ATLANTIC = {
    "NA":   BOXES["NA"],
    "Trop": BOXES["Trop"],
}

# Narrow-band version for panel b: each box spans its full latitude range but
# is restricted to a 5°-wide longitude strip centred on the section longitude
# (~28°W).  This symbolises that the Boussinesq model is a 2-D (lat × depth)
# cross-section through the Atlantic.
_SECTION_LON_CTR = -28.0   # centre of the strip (°E)
_BAND_WIDTH = 5.0           # degrees longitude
BOXES_BOUS_NARROW = {
    "NA": dict(
        lat_min=37.5, lat_max=90.0,
        lon_min=_SECTION_LON_CTR - _BAND_WIDTH / 2,
        lon_max=_SECTION_LON_CTR + _BAND_WIDTH / 2,
        depth_max=None, color=BOX_COLOR_NA, label="WOOD_NA",
    ),
    "Trop": dict(
        lat_min=-47.5, lat_max=37.5,
        lon_min=_SECTION_LON_CTR - _BAND_WIDTH / 2,
        lon_max=_SECTION_LON_CTR + _BAND_WIDTH / 2,
        depth_max=875.0, color=BOX_COLOR_TROP, label="WOOD_TROP",
    ),
    "South": dict(
        lat_min=-90.0, lat_max=-47.5,
        lon_min=_SECTION_LON_CTR - _BAND_WIDTH / 2,
        lon_max=_SECTION_LON_CTR + _BAND_WIDTH / 2,
        depth_max=450.0, color=BOX_COLOR_SOUTH, label="WOOD_SOUTH",
    ),
}

BOXES_SHALLOW = {
    "NA":    dict(lat_min=37.5,  lat_max=90.0,  lon_min=-72.5, lon_max=20.0,
                  depth_max=150.0, color=BOX_COLOR_NA,    label="WOOD_NA_SHALLOW"),
    "Trop":  dict(lat_min=-47.5, lat_max=32.5,  lon_min=-72.5, lon_max=20.0,
                  depth_max=150.0, color=BOX_COLOR_TROP,  label="WOOD_TROP_SHALLOW"),
    "South": dict(lat_min=-90.0, lat_max=-52.5, lon_min=-180.0, lon_max=180.0,
                  depth_max=150.0, color=BOX_COLOR_SOUTH, label="WOOD_SOUTH_SHALLOW"),
}

# Current-standard Atlantic-mask boxes (box_na / box_trop / box_south).
# NA and Trop select the Atlantic basin mask (basin="atlantic") intersected with
# a latitude band; the longitudinal edges are coastlines (no longitude limits).
# All three use the top 4 ocean layers (0–105 m).  Latitude tests use the ±35°
# grid borders so the boxes tile exactly (NA >35°N, Trop 35°S–35°N, South <35°S).
BOXES_CLIMBERX = {
    "NA":    dict(lat_min=35.0,  lat_max=90.0,  lon_min=-180.0, lon_max=180.0,
                  depth_max=105.0, basin="atlantic", color=BOX_COLOR_NA,    label="box_na"),
    "Trop":  dict(lat_min=-35.0, lat_max=35.0,  lon_min=-180.0, lon_max=180.0,
                  depth_max=105.0, basin="atlantic", color=BOX_COLOR_TROP,  label="box_trop"),
    "South": dict(lat_min=-90.0, lat_max=-35.0, lon_min=-180.0, lon_max=180.0,
                  depth_max=105.0, basin=None,       color=BOX_COLOR_SOUTH, label="box_south"),
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


def box_hmask(box):
    """Horizontal (lat, lon) selection for a box, honouring an Atlantic basin mask."""
    hm = horiz_mask(box["lat_min"], box["lat_max"], box["lon_min"], box["lon_max"])
    if box.get("basin") == "atlantic":
        hm = hm & ATL_MASK
    return hm


def lat_only_taper(box, hm):
    """Latitudinal-only cosine taper (for basin boxes): constant across longitude."""
    band = np.where((LAT >= box["lat_min"]) & (LAT <= box["lat_max"]))[0]
    w = np.zeros(hm.shape, dtype=float)
    if band.size == 0:
        return w
    bs = band[np.argsort(LAT[band])]          # south → north
    n = bs.size
    m = float(TAPER_MARGIN_CELLS)
    for p, j in enumerate(bs):
        d = min(p, n - 1 - p) + 1
        w[j, :] = 1.0 if d >= m else 0.5 * (1.0 - np.cos(np.pi * d / m))
    w[~hm] = 0.0
    return w


def make_3d_mask(box, taper=False):
    hm   = box_hmask(box)
    vm   = vert_mask(box["depth_max"])
    full = vm[:, None, None] & hm[None, :, :]
    if not taper:
        return full.astype(float)
    # Basin boxes taper in latitude only (E/W edges are coastlines).
    hw  = lat_only_taper(box, hm) if box.get("basin") else horiz_taper(hm)
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
#                    z_model ∈ [0, 1]: z_model=0 = bottom, z_model=1 = surface
# Actual depth grid is tanh-stretched: z_j(j) = 0.5 + tanh(q*(j/N−0.5))/(2·tanh(q/2))
# ---------------------------------------------------------------------------
BOUS_M, BOUS_N   = 20, 40
BOUS_XX          = np.linspace(0.0, 5.0, BOUS_M + 1)   # (21,) latitude coords

# Actual tanh-stretched depth grid (q=3, same as Boussinesq_2DAMOC.py)
_q = 3
_j = np.arange(BOUS_N + 1)
BOUS_Z_MODEL = 0.5 + np.tanh(_q * (_j / BOUS_N - 0.5)) / (2 * np.tanh(_q / 2))
# z_model=0 → bottom, z_model=1 → surface
# For plotting: surface at top (0 m), so flip to z_plot = 1 - z_model
BOUS_Z_PLOT  = 1.0 - BOUS_Z_MODEL   # 0=surface, 1=bottom (for plotting)


# Geographic bounds of Boussinesq boxes (matching Boussinesq_box.py __main__)
# Latitudinal borders snapped to grid cell faces (midpoints between nodes).
BOUS_NORTH_LO = 3.625           # +40.5°N  (cell face)
BOUS_NORTH_HI = 5.0             #  90°N
BOUS_TROP_LO  = 1.625           # −31.5°S  (cell face)
BOUS_TROP_HI  = 3.625           # +40.5°N  (cell face)
BOUS_SOUTH_LO = -1.0            #  (below domain → clamped to 90°S)
BOUS_SOUTH_HI = 1.625           # −31.5°S  (cell face)

# Vertical extent of the finite surface boxes (z-units), matching Boussinesq_box.py:
# cells with z_model ≥ 1 − BOUS_BOX_DEPTH are inside the box.  0.28449 places the
# box bottom on the grid cell face at ~284 m (BOUS_DEPTH_SCALE = 1000 m),
# so the box reaches the ~267 m model level (top 15 levels).
BOUS_BOX_DEPTH = 0.28449

# Build the actual model so the figure stays in sync with the perturbation logic
# in AMOCBoussinesq/Boussinesq_box.py (finite membership + margin-only cosine
# taper + interior boost so the box-mean of the field is 1).
sys.path.insert(0, str(SCRIPT_DIR / "AMOCBoussinesq"))
from Boussinesq_box import BoussinesqBox

_bous_model = BoussinesqBox(BOUS_M, BOUS_N, 1e-3)
_bous_model.make_salinity_forcing(0.1)
_bous_model.register_box("north",    BOUS_NORTH_LO, BOUS_NORTH_HI, box_depth=BOUS_BOX_DEPTH)
_bous_model.register_box("tropical", BOUS_TROP_LO,  BOUS_TROP_HI,  box_depth=BOUS_BOX_DEPTH)
_bous_model.register_box("southern", BOUS_SOUTH_LO, BOUS_SOUTH_HI, box_depth=BOUS_BOX_DEPTH)

# Tapered, interior-boosted perturbation fields, shape (M+1, N+1)
BOUS_NORTH_PERT = _bous_model._box_pert["north"]
BOUS_TROP_PERT  = _bous_model._box_pert["tropical"]
BOUS_SOUTH_PERT = _bous_model._box_pert["southern"]

# Convert Boussinesq x to degrees latitude for axis labels
BOUS_LAT_DEG = BOUS_XX / 5.0 * 180.0 - 90.0    # (21,)


# ---------------------------------------------------------------------------
# Globe-panel setup helpers
# ---------------------------------------------------------------------------

_GLOBE_PROJ    = ccrs.Orthographic(-30, 5) if HAS_CARTOPY else None
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
        sel = lat_ok[:, None] & lon_ok[None, :] & FINE_OCEAN_MASK
        if box.get("basin") == "atlantic":
            # Stretch the box longitudinally to the basin coastlines instead of
            # the nominal lon_min/lon_max rectangle.  South of the basin mask's
            # extent (open South Atlantic) the basin mask is empty, so there we
            # keep the longitude-rectangle fill so the box still reaches its
            # southern latitude bound (where the Southern box begins).  The
            # extension is clamped to the Atlantic's longitudinal span at its
            # southern edge so it does not leak into the Pacific west of Chile.
            ext_lon_ok = (lon_norm >= ATL_SOUTH_LON_MIN) & (lon_norm <= ATL_SOUTH_LON_MAX)
            south_ext = (_FINE_LAT < ATL_SOUTH_EDGE)[:, None] & ext_lon_ok[None, :]
            sel = sel & (FINE_ATL_MASK | south_ext)
        w2 = sel.astype(float)
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


def draw_section_strip_projected(ax, box_dict, lon=_SECTION_LON_CTR,
                                  strip_half_width_m=300_000):
    """Draw each box as a constant-projected-width strip around the section longitude.

    Selects fine-grid cells whose orthographic x-coordinate is within
    strip_half_width_m metres of the section line.  Because the selection is
    done in projected space (not longitude degrees), the strip has constant
    apparent width on screen regardless of latitude.

    Uses the same pcolormesh path as draw_boxes_fine_grid so globe background
    rendering is fully preserved.
    """
    R = 6.371e6  # Earth radius, metres
    LON0_RAD = np.radians(-30.0)    # orthographic projection centre longitude
    lon_sec_rad = np.radians(lon)

    lon2d_rad = np.radians(np.meshgrid(_FINE_LON, _FINE_LAT)[0])   # (lat, lon)
    lat2d_rad = np.radians(np.meshgrid(_FINE_LON, _FINE_LAT)[1])

    # Orthographic x for every grid cell and for the section line at each latitude
    x_cells   = R * np.cos(lat2d_rad) * np.sin(lon2d_rad  - LON0_RAD)
    x_section = R * np.cos(lat2d_rad) * np.sin(lon_sec_rad - LON0_RAD)

    in_strip = np.abs(x_cells - x_section) < strip_half_width_m   # (lat, lon) bool

    for key, box in box_dict.items():
        lat_ok = (_FINE_LAT >= box["lat_min"]) & (_FINE_LAT <= box["lat_max"])
        w2 = (lat_ok[:, None] & in_strip & FINE_OCEAN_MASK).astype(float)

        rgba = np.array(mcolors.to_rgba(box["color"]))
        cmap = LinearSegmentedColormap.from_list(
            "", [(rgba[0], rgba[1], rgba[2], 0.0), (*rgba[:3], 0.75)], N=2
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

def draw_section(ax, box_dict, taper=False, depth_max_plot=None, label_depths=True,
                 box_limits=False, depth_lines=True, lat_lines=True):
    """
    Meridional cross-section at SECTION_LON.
    depth_max_plot: clip y-axis to this depth (m); None = full depth.
    box_limits: if True, outline each box's finite extent (latitude range ×
        surface→depth_max) as a dashed rectangle instead of loose edge lines.
    depth_lines: if True, draw a horizontal dashed line at each box's depth_max
        (the depth text label is controlled separately by label_depths).
    lat_lines: if True, draw the vertical box-edge lines and the gray centre line.
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

    depth_lim = depth_max_plot if depth_max_plot is not None else float(ZRO[-1] + 300)

    # Depth boundary labels
    if label_depths:
        for key, box in box_dict.items():
            dmax = box.get("depth_max")
            if dmax is not None:
                if depth_lines and not box_limits:
                    ax.axhline(dmax, color=box["color"], linewidth=1.0,
                               linestyle="--", alpha=0.7, zorder=5)
                ax.text(87.0, dmax + 25, f"{int(dmax)} m",
                        fontsize=6.5, color=box["color"], ha="right", va="top")

    if box_limits:
        # Finite box limits as dashed rectangles: lat range × [surface, depth_max]
        for key, box in box_dict.items():
            dmax = box.get("depth_max")
            depth_bottom = dmax if dmax is not None else depth_lim
            lat_lo = max(box["lat_min"], -90.0)
            lat_hi = min(box["lat_max"],  90.0)
            ax.add_patch(mpatches.Rectangle(
                (lat_lo, 0.0), lat_hi - lat_lo, depth_bottom,
                fill=False, edgecolor=box["color"], linewidth=1.3,
                linestyle="--", alpha=0.9, zorder=6))
    elif lat_lines:
        # Latitude boundary lines
        for key, box in box_dict.items():
            for lat_edge in [box["lat_min"], box["lat_max"]]:
                ax.axvline(lat_edge, color=box["color"], linewidth=0.8,
                           linestyle=":", alpha=0.5, zorder=5)

    ax.set_xlim(-90, 90)
    ax.set_ylim(depth_lim, 0)
    ax.set_facecolor("#f0f0f0")
    ax.set_xlabel("Latitude", fontsize=8)
    if lat_lines:
        ax.axvline(0, color="gray", linewidth=0.4, linestyle="--", alpha=0.5)

    ticks = np.arange(-90, 91, 30)
    ax.set_xticks(ticks)
    ax.set_xticklabels(
        [f"{abs(t)}°{'N' if t >= 0 else 'S'}" for t in ticks], fontsize=6.5
    )


# ---------------------------------------------------------------------------
# Boussinesq 2D panel
# ---------------------------------------------------------------------------

BOUS_DEPTH_SCALE = 1000.0   # normalised depth 1 corresponds to this many metres

def _bous_lat_deg(x):
    """Boussinesq model latitude coord x ∈ [0, 5] → degrees north, clamped to domain."""
    return np.clip(x, 0.0, 5.0) / 5.0 * 180.0 - 90.0


def draw_boussinesq_panel(ax):
    """
    pcolormesh plot of the Boussinesq finite tapered perturbation fields using
    the model's actual tanh-stretched depth grid.  Each box field is the
    margin-only cosine taper with the interior boost (box-mean = 1), normalised
    to its peak for display so the taper ramp is visible.  The finite box limits
    (latitude range × surface→box-bottom depth) are drawn as dashed rectangles.
    Y-axis is in physical metres (z_plot × BOUS_DEPTH_SCALE) to match the shared
    axis with the meridional section panels.
    """
    from matplotlib.colors import LinearSegmentedColormap as LSC

    # Convert tanh-stretched z_plot to metres (0=surface, 1=bottom in plot)
    ZZ_m = BOUS_Z_PLOT * BOUS_DEPTH_SCALE   # (41,) in metres

    pert_n = BOUS_NORTH_PERT.T   # (N+1, M+1)
    pert_t = BOUS_TROP_PERT.T
    pert_s = BOUS_SOUTH_PERT.T

    for pert, color, zorder in [
        (pert_s, BOX_COLOR_SOUTH, 2),
        (pert_t, BOX_COLOR_TROP,  3),
        (pert_n, BOX_COLOR_NA,    4),
    ]:
        disp = pert / pert.max() if pert.max() > 0 else pert   # normalise for display
        rgba = np.array(mcolors.to_rgba(color))
        cmap = LSC.from_list("", [(rgba[0], rgba[1], rgba[2], 0.0),
                                   (*rgba[:3], 0.95)], N=256)
        # Node coords are cell *centres*: shading="nearest" puts cell boundaries
        # on the midpoints (the grid cell faces), so the painted cells line up
        # exactly with the face-aligned box-limit rectangles drawn below.
        ax.pcolormesh(BOUS_LAT_DEG, ZZ_m, disp,
                      cmap=cmap, vmin=0, vmax=1,
                      shading="nearest", zorder=zorder)

    # Finite box limits: latitude range × [surface, box bottom] as dashed rectangles
    box_bottom_m = BOUS_BOX_DEPTH * BOUS_DEPTH_SCALE
    for (x_lo, x_hi), col in [
        ((BOUS_NORTH_LO, BOUS_NORTH_HI), BOX_COLOR_NA),
        ((BOUS_TROP_LO,  BOUS_TROP_HI),  BOX_COLOR_TROP),
        ((BOUS_SOUTH_LO, BOUS_SOUTH_HI), BOX_COLOR_SOUTH),
    ]:
        lat_lo, lat_hi = _bous_lat_deg(x_lo), _bous_lat_deg(x_hi)
        ax.add_patch(mpatches.Rectangle(
            (lat_lo, 0.0), lat_hi - lat_lo, box_bottom_m,
            fill=False, edgecolor=col, linewidth=1.3, linestyle="--",
            alpha=0.9, zorder=6))

    ax.set_facecolor("#f0f0f0")
    ax.set_xlim(-90, 90)
    # ylim set externally to match shared axis (1000 m); depth increases downward
    ax.set_xlabel("Latitude", fontsize=8)
    ax.set_ylabel("Depth (m)", fontsize=8)

    lat_ticks = np.arange(-90, 91, 30)
    ax.set_xticks(lat_ticks)
    ax.set_xticklabels(
        [f"{abs(t)}°{'N' if t >= 0 else 'S'}" for t in lat_ticks], fontsize=6.5
    )

    ax.axvline(0, color="gray", linewidth=0.4, linestyle="--", alpha=0.5)


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
    draw_section(ax, BOXES, taper=False, depth_max_plot=2000, label_depths=False,
                 depth_lines=False, lat_lines=False)
    ax.set_ylabel("Depth (m)", fontsize=8)
    ax.set_title("Meridional section (28°W)", fontsize=8)
    add_panel_label(ax, panel_labels_bot[0], y=0.04, va="bottom")

    # ── Column 2: Boussinesq context ──────────────────────────────────────
    ax = axes_top[1]
    setup_globe(ax)
    draw_section_strip_projected(ax, BOXES_BOUS_NARROW)   # constant-width strip → 2-D cross-section
    ax.set_title("Regions corresponding to Boussinesq model boxes", fontsize=8,
                 fontweight="bold")
    add_panel_label(ax, panel_labels_top[1])

    ax = axes_bot[1]
    draw_boussinesq_panel(ax)
    ax.set_title("Boussinesq 2D tapered perturbation", fontsize=8)
    add_panel_label(ax, panel_labels_bot[1], y=0.04, va="bottom")

    # ── Column 3: CLIMBER-X ───────────────────────────────────────────────
    # Atlantic-mask boxes box_na / box_trop / box_south (top 4 layers, 0–105 m);
    # NA & Trop follow the Atlantic basin mask, tapered in latitude only.
    ax = axes_top[2]
    setup_globe(ax)
    draw_boxes_on_globe(ax, BOXES_CLIMBERX, taper=True)
    ax.set_title("CLIMBER-X model boxes", fontsize=8, fontweight="bold")
    add_panel_label(ax, panel_labels_top[2])

    ax = axes_bot[2]
    draw_section(ax, BOXES_CLIMBERX, taper=True, depth_max_plot=300, label_depths=True,
                 box_limits=True)
    ax.set_ylabel("Depth (m)", fontsize=8)
    ax.set_title("Meridional section (28°W)", fontsize=8)
    add_panel_label(ax, panel_labels_bot[2], y=0.04, va="bottom")

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
    add_panel_label(ax, panel_labels_bot[3], y=0.04, va="bottom")

    # ── Shared y-axis for bottom panels e, f, g (depth 0–1000 m) ──────────
    for ax in axes_bot[:3]:
        ax.set_ylim(1000, 0)
    # Only the leftmost panel (e) keeps the y-axis label and tick labels;
    # f, g, h suppress both.
    for ax in axes_bot[1:]:
        ax.set_ylabel("")
        ax.tick_params(labelleft=False)

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
