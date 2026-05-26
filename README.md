# AMOCResilience

This repository collects four subprojects that together investigate the resilience of the Atlantic Meridional Overturning Circulation (AMOC) across a hierarchy of climate models.

The central question is whether different notions of resilience yield consistent signals under climate change — and whether they can serve as robust early-warning indicators of an AMOC tipping point.

## Subprojects

### [AMOCBox](https://github.com/andreasmorr/AMOCBox)
Resilience analysis of the Alkhayuon et al. 3-box AMOC model. The model is parametrized for 1×CO₂ and 2×CO₂ settings; interpolating between them gives a continuum of AMOC states under increasing CO₂. Multiple stability measures from [Attractors.jl](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/attractors/stable/) are computed along this continuum to identify coherent resilience trends.

### [AMOCPlaSim](https://github.com/andreasmorr/AMOCPlaSim)
Edge-tracking and equilibrium analysis of the AMOC in the general circulation model PlaSim. Bisection trajectories from the saddle (edge) state to either the AMOC-on or AMOC-off attractor are available at pre-industrial (285 ppm) and current (360 ppm) CO₂ levels. Gaussian covariance ellipsoids fitted to equilibrium runs in EOF space are used to compute convergence times, edge-to-attractor distances, and local stability metrics.

### [AMOCBoussinesq](https://github.com/andreasmorr/AMOCBoussinesq)
Geographic box perturbation sweep of the 2D latitude-depth Boussinesq model. Salinity perturbations are applied to coherent geographic boxes (Arctic, tropical, southern ocean) across a range of Arctic amplification values γ, mapping out basins of attraction and mean convergence times in the Boussinesq AMOC.

### [AMOCClimberX](https://github.com/andreasmorr/AMOCClimberX)
Trajectory-based basin-of-attraction mapping in the intermediate-complexity model CLIMBER-X. Many runs are started from perturbed salinity states to characterize basin geometry and transient dynamics. Equilibrium runs additionally provide local linear stability estimates comparable to the box model results.

## Structure

This repo uses [git submodules](https://git-scm.com/book/en/v2/Git-Tools-Submodules). Clone with:

```bash
git clone --recurse-submodules https://github.com/andreasmorr/AMOCResilience.git
```

Or, if already cloned:

```bash
git submodule update --init --recursive
```
