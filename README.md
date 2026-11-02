# AMOCResilience

This repository collects four subprojects that together investigate the
resilience of the Atlantic Meridional Overturning Circulation (AMOC) across a
hierarchy of climate models.

The central question is whether different notions of resilience yield
consistent signals under climate change, and whether they can serve as robust
early-warning indicators of an AMOC tipping point.

The four submodules each contain a `PAPERINFO.md` file. Treat those files as the
paper-facing source of truth for model-specific methods, metric definitions,
and wording. The root files collect the paper draft and the cross-model figures.

## Subprojects

### [AMOCBox](https://github.com/andreasmorr/AMOCBox)
Resilience analysis of the reduced three-box AMOC salinity model. The current
paper outputs follow a linear parameter continuation from 1xCO2 (280 ppm)
toward and beyond the 2xCO2 parameter vector in 14 ppm steps, with one extra
point at 903 ppm. Complete paper-facing resilience metrics are exported from
280 through 903 ppm; the on-state is lost before the next step (910 ppm). At
each step, 10,000 initial conditions are sampled uniformly in a +/-2 psu box
centred on the on-state. The submodule contributes AMOC strength, local
resilience, convergence time, basin volume, and minimal critical shock metrics
for the conceptual-model member of the hierarchy.

### [AMOCPlaSim](https://github.com/andreasmorr/AMOCPlaSim)
Analysis of existing PlaSim-LSG edge-track and equilibrium output at 285 ppm
and 360 ppm CO2 (Börner et al.; doi:10.5281/zenodo.17053348). The current
reduction uses custom deep-box salinity readouts for the North Atlantic
(35-80N, 0-1000 m) and Tropical Atlantic (35S-35N, 0-500 m), not the removed
EOF coordinate workflow. 4-sigma Gaussian covariance ellipses in this
two-dimensional salinity space are used for convergence times,
edge-to-attractor distances, and the covariance-based local resilience proxy.
AMOC strength (maximum Atlantic overturning between 46N and 66N) comes from a
separate sidecar file. No new PlaSim runs are performed and no basin volume is
estimated.

### [AMOCBoussinesq](https://github.com/andreasmorr/AMOCBoussinesq)
Equilibrium-field salinity perturbation experiment in the 2D latitude-depth
Boussinesq model. At each forcing value, the complete stationary ON fields
are perturbed on a 20 x 20 grid over +/-2 display-salinity units in northern
and tropical boxes extending to 0.28449 of the basin depth. Salt compensation
is distributed outside both boxes, including the deeper ocean. The canonical
experiment requests 20 gamma values from 0 to 0.11 with denser upper-end
spacing. The completed run contains 6,800 trajectories and 17 controls, with
verified ON references through gamma = 0.1085; metrics at 0.109, 0.1095 and
0.11 are unavailable because stationary correction did not verify an ON
reference. This is not evidence of a zero-sized basin or a certified fold.
Full fields distinguish ON, reversed and two-cell outcomes. Local resilience
comes from the full-state tangent spectrum, and convergence times use the
last observed entry into a fixed-normalization full-field ON neighbourhood.
Canonical results retain model units. The paper figures preserve their
committed layouts and shared axes, using fixed reporting conversions in a
separate paper export. One model time unit is reported as 100 model years;
return times multiply by 100 and inverse-time resilience rates divide by 100.
The salinity-gradient conversion uses the same affine PSU display convention
as the perturbation grid, with its offset cancelling in the difference.
Critical shocks already use this coordinate and are not rescaled. Overturning
uses the existing fixed Sv display factor. These are reporting conventions,
not a simultaneous physical calibration of the idealized model.
Its CO2 coordinate `C = 280 * 2^(gamma / 0.05)` is an illustrative display
relation, not a calibrated atmospheric CO2 response.

### [AMOCClimberX](https://github.com/andreasmorr/AMOCClimberX)
Trajectory-based perturbation experiments in the Earth system model of
intermediate complexity CLIMBER-X. Fixed-CO2 equilibria at 15 levels between
285 and 600 ppm are perturbed in North Atlantic (35-80N) and Tropical Atlantic
(35S-35N) shallow salinity boxes (upper 105 m) on a 9 x 9 grid over +/-2 psu,
and each run is integrated for 10,000 years. The added salt is compensated
uniformly over the rest of the upper 105 m of the world ocean (global
compensation). A dedicated 9-run stability ensemble per CO2 level provides
local linear stability estimates via a fitted 2 x 2 discrete linear map. The
reference equilibrium lies on the "modern" branch up to 330 ppm and on a
stronger-AMOC ("strong") branch from 345 ppm onward; pure salinity
perturbations do not reach an AMOC off-state. Raw model output lives on the PIK
HPC system, so `plotting_paper.py` has to be run there.

## Umbrella-level files

Several shared Python files live at the repository root and support the paper
figures:

| File | Purpose |
|------|---------|
| `amoc_plot_style.py` | Shared matplotlib style: color constants, `make_paper_figure`, `add_panel_label`, `savefig_pdf`. Imported by every `plotting_paper.py`. |
| `synthesis_figure.py` | Cross-model resilience figure in the committed six-panel layout, with AMOC strength, Tropical-minus-North-Atlantic salinity gradient, local resilience, convergence time, basin volume, and critical shock. Boussinesq uses fixed paper-unit conversions on the existing shared axes and a heuristic CO2 coordinate. Reads paper CSV exports and PlaSim's summary CSV; writes PNG, PDF and input-hash provenance. |
| `plotting_perturbations.py` | Multi-model perturbation and readout overview figure showing the geographic perturbation/readout geometry and reduced PlaSim salinity coordinates. Requires `cartopy` for map panels and falls back to a flat map if needed. |
| `basin_mask_5x5.nc` | Small Atlantic-basin mask used by `plotting_perturbations.py`. |
| `main.tex` | Paper draft. It is useful for narrative context, but model-specific details should be checked against each submodule's `PAPERINFO.md`. |

Run umbrella figures from the repository root:
```bash
python AMOCBoussinesq/scripts/plotting_paper.py
python AMOCBoussinesq/scripts/export_paper_resilience.py
python synthesis_figure.py
python plotting_perturbations.py
```

The Boussinesq four-panel paper figure compares gamma = 0 and 0.09, with
four shared grid points selected by `AMOCBoussinesq/configs/paper_figure.json`:
`(3,3), (15,16), (18,18), (17,10)`. They give reversed-to-two-cell, ON-to-two-cell,
reversed-to-two-cell and ON-to-ON transitions, respectively. The committed
panel geometry, trajectory colours, absolute PSU axes and reference circles
are retained. The figure writes `AMOCBoussinesq/plots/paper/boussinesq_paper.png`,
a matching PDF and provenance JSON.

For smooth trajectories, run
`python AMOCBoussinesq/scripts/refine_paper_trajectories.py` before plotting.
This reruns only the eight selected cases from their exact saved initial fields,
with unchanged timestep .001 and readouts at every step (0.1 model years).
The plot automatically selects the matching completed dense output. All original
sampled observations and final fields are checked against the canonical run.
The full ensemble and synthesis statistics retain their .25 observation spacing.
All observations from time zero are included. Reference circles are projected
visual neighbourhoods; endpoint classification uses all model fields.

The Boussinesq export reads the completed canonical run selected by
`AMOCBoussinesq/data/paper_experiment/latest.json`, preserves each original
numerical measure and unit in `model_value` and `model_units`, and writes
`AMOCBoussinesq/data/paper/resilience_vs_gamma_boussinesq.csv` with provenance.
The separate `value` and `units` columns hold paper conversions. All 20
requested gamma values are retained; unavailable references have blank values.
The fixed factors are salinity slope 10.7273456518, streamfunction factor
18/1.5171475119, and 100 model years per unit model time. Under the
[`H^2/kappa_T` time scale of Soons et al. (2025)](https://doi.org/10.1017/jfm.2025.248),
the latter corresponds formally to an effective
diffusivity of approximately 3.17e-4 m²/s at a nominal 1000 m depth; it preserves
the committed reporting convention rather than fitting the resilience curves.
Pass `--input` to the export script to select a run explicitly. The synthesis
script accepts `--boussinesq-csv`, `--output`, and `--umbrella` overrides.
Other models retain their existing exported values; their local-resilience
and critical-shock definitions differ and should be interpreted using each
subproject's `PAPERINFO.md`.

Current paper-facing figure outputs are written below `plots/`, which is ignored
by git.

## Structure

This repo uses [git submodules](https://git-scm.com/book/en/v2/Git-Tools-Submodules). Clone with:

```bash
git clone --recurse-submodules https://github.com/andreasmorr/AMOCResilience.git
```

Or, if already cloned:

```bash
git submodule update --init --recursive
```
