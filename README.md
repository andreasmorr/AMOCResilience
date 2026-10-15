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
Geographic box perturbation sweep of the 2D latitude-depth Boussinesq model.
Salinity perturbations are applied on a 20 x 20 grid over +/-2 PSU in North
Atlantic and Tropical Atlantic box salinity (top 15 grid levels), with a
Southern box compensation. The deterministic experiment sweeps
`gamma = 0.00 ... 0.06`, equivalent to 280-643 ppm CO2 under the adopted
mapping `gamma = 0.05 log2(C / 280 ppm)` (2xCO2 is `gamma = 0.05`). A cached
`gamma = 0.07` sweep (~739 ppm) is excluded from all paper outputs because most
of its runs time out. Local resilience comes from a separate eight-direction
return-time ensemble around the symmetric-forcing (`beta = 0`) on-state, and
the reported AMOC strength is that of this on-state.

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
| `synthesis_figure.py` | Cross-model resilience vs CO2 figure with AMOC strength and Tropical-minus-North-Atlantic salinity gradient panels plus local resilience, convergence time, basin volume, and minimal critical shock panels. Reads the long-format paper CSV exports from the submodules plus PlaSim's summary CSV. |
| `plotting_perturbations.py` | Multi-model perturbation and readout overview figure showing the geographic perturbation/readout geometry and reduced PlaSim salinity coordinates. Requires `cartopy` for map panels and falls back to a flat map if needed. |
| `basin_mask_5x5.nc` | Small Atlantic-basin mask used by `plotting_perturbations.py`. |
| `main.tex` | Paper draft. It is useful for narrative context, but model-specific details should be checked against each submodule's `PAPERINFO.md`. |

Run umbrella figures from the repository root:
```bash
python synthesis_figure.py
python plotting_perturbations.py
```

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
