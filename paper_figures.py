"""Publish only manuscript figures, without copying code or data into the paper.

Run ``python paper_figures.py sync`` after pulling results from another machine,
or ``python paper_figures.py check`` to check dependencies and figure freshness.
Set AMOC_PAPER_DIR to use a separate checkout of AMOCResiliencePaper.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import re
import shutil
import tempfile

ROOT = Path(__file__).resolve().parent
FIGURES = {
    "perturbations_overview.png": "plots/perturbations_overview.png",
    "amocbox_paper.png": "AMOCBox/plots/amocbox_paper.png",
    "boussinesq_paper.png": "AMOCBoussinesq/plots/paper/boussinesq_paper.png",
    "climberx_paper.png": "AMOCClimberX/plots/climberx_paper.png",
    "plasim_paper.png": "AMOCPlaSim/plots/plasim_paper.png",
    "synthesis_figure.png": "plots/synthesis_figure.png",
    "boussinesq_equilibria.pdf": "AMOCBoussinesq/plots/paper/boussinesq_equilibria.pdf",
    "box_traj_315_10ky_gc_narrow.png": "AMOCClimberX/plots/box_traj_315_10ky_gc_narrow.png",
}


def paper_directory() -> Path:
    return Path(os.environ.get("AMOC_PAPER_DIR", ROOT / "paper")).expanduser().resolve()


def publish_paper_figure(source: str | Path, *, name: str | None = None) -> Path:
    """Atomically mirror an approved figure; leave the original plot in place.

    Plotters decide whether their selected experiment belongs to the manuscript.
    The explicit allowlist prevents diagnostic outputs entering the paper repo.
    """
    source = Path(source).resolve()
    name = source.name if name is None else name
    if name not in FIGURES:
        raise ValueError(f"Not a manuscript figure: {name}")
    if not source.is_file() or source.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty figure: {source}")
    destination = paper_directory() / "figures" / name
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and source.read_bytes() == destination.read_bytes():
        return destination
    # Replace atomically so interrupted copying cannot leave a partial image.
    with tempfile.NamedTemporaryFile(dir=destination.parent, delete=False) as handle:
        temporary = Path(handle.name)
    try:
        shutil.copyfile(source, temporary)
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)
    print(f"Paper figure saved: {destination}")
    return destination


def sync_figures() -> None:
    # Check all sources before changing any paper copies.
    missing = [str(ROOT / path) for path in FIGURES.values()
               if not (ROOT / path).is_file() or (ROOT / path).stat().st_size == 0]
    if missing:
        raise FileNotFoundError("Missing paper figure sources:\n" + "\n".join(missing))
    for name, path in FIGURES.items():
        publish_paper_figure(ROOT / path, name=name)


def check_figures() -> None:
    directory = paper_directory()
    errors = []
    referenced = set()
    for tex_name in ("main.tex", "SI.tex"):
        tex = directory / tex_name
        if not tex.is_file():
            errors.append(f"Missing manuscript: {tex}")
            continue
        content = re.sub(r"(?<!\\)%[^\n]*", "", tex.read_text())
        for name in re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", content):
            referenced.add(name)
            figure = directory / "figures" / name
            if not figure.is_file() or figure.stat().st_size == 0:
                errors.append(f"{tex_name}: missing or empty figure {name}")
        for bibliography in re.findall(r"\\bibliography\{([^}]+)\}", content):
            for name in bibliography.split(","):
                if not (directory / f"{name.strip()}.bib").is_file():
                    errors.append(f"{tex_name}: missing bibliography {name}")
    for name in referenced:
        if name not in FIGURES:
            errors.append(f"Figure has no registered experiment source: {name}")
            continue
        source, destination = ROOT / FIGURES[name], directory / "figures" / name
        if not source.is_file():
            errors.append(f"Missing original figure: {source}")
        elif destination.is_file() and source.read_bytes() != destination.read_bytes():
            errors.append(f"Paper copy differs from plot: {name}; run sync")
    if errors:
        raise ValueError("\n".join(errors))
    print(f"Verified {len(referenced)} manuscript figures and both bibliographies.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("sync", "check"))
    args = parser.parse_args()
    if args.command == "sync":
        sync_figures()
    check_figures()


if __name__ == "__main__":
    main()
