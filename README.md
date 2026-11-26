# AMOCResilience

Experiments and manuscript for AMOC resilience across four climate models.
The four model repositories retain their experiment code, numerical results,
and `PAPERINFO.md` descriptions. The separate private
[AMOCResiliencePaper](https://github.com/andreasmorr/AMOCResiliencePaper)
repository at `paper/` contains the manuscript and its figures for Overleaf.

## Repository layout

| Path | Purpose |
| --- | --- |
| `AMOCBox/` | Three-box salinity model, continuation, and resilience analysis |
| `AMOCBoussinesq/` | Latitude-depth model, certified equilibria, perturbation and local-response experiments |
| `AMOCClimberX/` | CLIMBER-X experiment setup, HPC workflows, and analysis |
| `AMOCPlaSim/` | Analysis of existing PlaSim-LSG edge-track and equilibrium output |
| `paper/` | Independent manuscript repository: `main.tex`, `SI.tex`, bibliography, eight figures |
| `amoc_plot_style.py` | Shared plotting style |
| `paper_figures.py` | Explicit paper-figure registry, publication helper, synchronization and checks |
| `synthesis_figure.py` | Cross-model resilience figure and input-hash provenance |
| `plotting_perturbations.py` | Cross-model perturbation/readout geometry figure |
| `basin_mask_5x5.nc` | Basin mask required by the geometry figure |
| `plots/` | Local outputs of the two cross-model plotters |

Use each model's `PAPERINFO.md` as the source of truth for methods and metric
interpretation. Reproduction commands and data requirements are in its README.
Numerical results, source snapshots, run manifests, and scientific validation
remain with the models. Diagnostic-only plotting utilities and old build
products have been removed.

## Clone and manuscript editing

```sh
git clone --recurse-submodules https://github.com/andreasmorr/AMOCResilience.git
# For an existing checkout:
git submodule update --init --recursive
```

The `paper` submodule requires access to the private paper repository. Import
that repository into Overleaf and select `main.tex`; select `SI.tex` to build
the supplement. Both compile independently of this experiment workspace.

Edit the TeX files only in `paper/` or Overleaf. Before updating figures, sync
Overleaf edits back to GitHub and pull the paper repository:

```sh
git -C paper switch main
git -C paper pull --ff-only
```

Plot generation and synchronization change only paper figure files. They never
replace manuscript text or the bibliography. Commit and push changes inside
`paper/`, then record the updated paper commit in this parent repository. Model
changes likewise belong in their respective repositories before updating the
parent's submodule references.

## Paper figure workflow

Run these from the workspace root, after producing the model exports described
in each model README:

```sh
python AMOCBox/scripts/plotting_paper.py
python AMOCPlaSim/scripts/plotting_paper.py
python AMOCBoussinesq/scripts/plotting_paper.py
python AMOCBoussinesq/scripts/plot_equilibria.py
python synthesis_figure.py
python plotting_perturbations.py
```

Each canonical plot saves to its original plotting directory and also copies
the exact asset used by the manuscript into `paper/figures/`. Local alternate
formats and provenance remain with the experiment. Custom Boussinesq/synthesis
output paths are kept local so exploratory figures do not replace paper assets.

CLIMBER-X paper regeneration needs raw output on the HPC system. Its current
defaults describe newer deep-box experiments, while this manuscript still
uses historical shallow experiments. Follow its README's explicit historical
command to update `climberx_paper.png`. The supplementary
`box_traj_315_10ky_gc_narrow.png` is retained as a historical snapshot in
`AMOCClimberX/plots/`; its generating revision is documented in that README.

After retrieving remote figure outputs, refresh all eight paper copies:

```sh
python paper_figures.py sync
python paper_figures.py check
```

To publish into another local paper checkout, set `AMOC_PAPER_DIR` to its
absolute path. No Git operation runs automatically when a figure is generated.
For HPC use, provide the shared `amoc_plot_style.py` and `paper_figures.py`
modules in the parent directory of the model checkout, as in this workspace.

Build from the standalone paper repository:

```sh
cd paper
latexmk -pdf main.tex
latexmk -pdf SI.tex
```

Before regenerating the synthesis figure, export Boussinesq paper units using
`python AMOCBoussinesq/scripts/export_paper_resilience.py`. Canonical model
values and provenance remain unchanged; see that model's `PAPERINFO.md` for
reporting conversions and the illustrative CO2 coordinate. Keep the CLIMBER-X
historical paper CSV selected until the manuscript is revised for deep-box runs.
