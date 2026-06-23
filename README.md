# AMOCResilience

This repository collects four subprojects that together investigate the resilience of the Atlantic Meridional Overturning Circulation (AMOC) across a hierarchy of climate models.

The central question is whether different notions of resilience yield consistent signals under climate change — and whether they can serve as robust early-warning indicators of an AMOC tipping point.

## Subprojects

### [AMOCBox](https://github.com/andreasmorr/AMOCBox)
Resilience analysis of the Alkhayuon et al. 3-box AMOC model. The model is parametrized for 1×CO₂ (280 ppm) and 2×CO₂ (560 ppm) settings; the parameter curve is interpolated between them and extrapolated beyond 2×CO₂ up to t=2.4 (≈952 ppm). Multiple stability measures from [Attractors.jl](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/attractors/stable/) are computed along this continuum to identify coherent resilience trends.

### [AMOCPlaSim](https://github.com/andreasmorr/AMOCPlaSim)
Edge-tracking and equilibrium analysis of the AMOC in the general circulation model PlaSim. Bisection trajectories from the saddle (edge) state to either the AMOC-on or AMOC-off attractor are available at pre-industrial (285 ppm) and current (360 ppm) CO₂ levels. Gaussian covariance ellipsoids fitted to equilibrium runs in EOF space are used to compute convergence times, edge-to-attractor distances, and local stability metrics.

### [AMOCBoussinesq](https://github.com/andreasmorr/AMOCBoussinesq)
Geographic box perturbation sweep of the 2D latitude-depth Boussinesq model. Salinity perturbations are applied to coherent geographic boxes (North Atlantic, tropical, southern ocean) across a range of Arctic amplification values γ, mapping out basins of attraction, mean convergence times, and local linear stability in the Boussinesq AMOC.

### [AMOCClimberX](https://github.com/andreasmorr/AMOCClimberX)
Trajectory-based basin-of-attraction mapping in the intermediate-complexity model CLIMBER-X. Many runs are started from perturbed salinity states to characterize basin geometry and transient dynamics. A dedicated 9-run stability ensemble (8 isotropic directions + 1 reference) additionally provides local linear stability estimates via a fitted 2×2 discrete linear map.

## Umbrella-level files

Several shared Python files live at the repository root and are used by all submodules:

| File | Purpose |
|------|---------|
| `amoc_plot_style.py` | Shared matplotlib style: color constants, `make_paper_figure`, `add_panel_label`, `savefig_pdf`. Imported by every `plotting_paper.py`. |
| `synthesis_figure.py` | Cross-model resilience vs CO₂ figure (5 panels): (a) AMOC strength (full-width top); then a 2×2 grid: (b) local resilience / characteristic return time, (c) convergence time, (d) basin volume, (e) minimal critical shock. Reads CSV exports from each submodule: `resilience_vs_co2_boxmodel.csv`, `resilience_vs_gamma_boussinesq.csv`, `resilience_vs_co2_climberx.csv`, `resilience_metrics.csv` (PlaSim). The Boussinesq AMOC strength is scaled to Sv anchored at 18 Sv pre-industrial; PlaSim uses `ellipse_long_axis_1sigma` as a local-resilience stand-in. |
| `plotting_perturbations.py` | Multi-model perturbation and readout overview figure (4 columns × 2 rows): Wood-box footprints on globe, meridional cross-sections, Boussinesq 2D domain, and PlaSim EOF panels. Requires `cartopy` (falls back to flat map). |

Run umbrella figures from the repository root:
```bash
python synthesis_figure.py
python plotting_perturbations.py
```

## Structure

This repo uses [git submodules](https://git-scm.com/book/en/v2/Git-Tools-Submodules). Clone with:

```bash
git clone --recurse-submodules https://github.com/andreasmorr/AMOCResilience.git
```

Or, if already cloned:

```bash
git submodule update --init --recursive
```
