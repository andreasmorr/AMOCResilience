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
toward and beyond the 2xCO2 parameter vector, with complete paper-facing
resilience metrics currently exported through 896 ppm. The submodule contributes
AMOC strength, local resilience, convergence time, basin volume, and minimal
critical shock metrics for the conceptual-model member of the hierarchy.

### [AMOCPlaSim](https://github.com/andreasmorr/AMOCPlaSim)
Analysis of existing PlaSim-LSG edge-track and equilibrium output at 285 ppm
and 360 ppm CO2. The current reduction uses custom deep-box salinity readouts
for the North Atlantic and Tropical Atlantic, not the removed EOF coordinate
workflow. Gaussian covariance ellipses in this two-dimensional salinity space
are used for convergence times, edge-to-attractor distances, and the
covariance-based local resilience proxy.

### [AMOCBoussinesq](https://github.com/andreasmorr/AMOCBoussinesq)
Geographic box perturbation sweep of the 2D latitude-depth Boussinesq model.
Salinity perturbations are applied on a 20 x 20 grid over +/-2 PSU in North
Atlantic and Tropical Atlantic box salinity, with a Southern box compensation.
The current deterministic experiment sweeps `gamma = 0.00 ... 0.07`, equivalent
to roughly 280-455 ppm CO2 under the adopted mapping.

### [AMOCClimberX](https://github.com/andreasmorr/AMOCClimberX)
Trajectory-based perturbation experiments in the intermediate-complexity model
CLIMBER-X. Fixed-CO2 equilibria are perturbed in North Atlantic and Tropical
Atlantic shallow salinity boxes on a 9 x 9 grid, with salt compensation in the
Southern Ocean. A dedicated 9-run stability ensemble provides local linear
stability estimates via a fitted 2 x 2 discrete linear map.

## Umbrella-level files

Several shared Python files live at the repository root and support the paper
figures:

| File | Purpose |
|------|---------|
| `amoc_plot_style.py` | Shared matplotlib style: color constants, `make_paper_figure`, `add_panel_label`, `savefig_pdf`. Imported by every `plotting_paper.py`. |
| `synthesis_figure.py` | Cross-model resilience vs CO2 figure with AMOC strength plus local resilience, convergence time, basin volume, and minimal critical shock panels. Reads the long-format paper CSV exports from the submodules plus PlaSim's summary CSV. |
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
