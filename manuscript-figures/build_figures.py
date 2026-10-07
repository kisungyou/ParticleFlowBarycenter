#!/usr/bin/env python3
"""Build all ten current manuscript figures from saved repository inputs."""
from pathlib import Path
import argparse
import hashlib
import importlib.metadata
import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
PUBLISHED = HERE / "published"
GROUPS = [
    ("gaussian_figures.py", ["fig-sim-gauss-2", "fig-sim-gauss-4"]),
    ("runtime_figure.py", ["fig-sim-gauss-3"]),
    ("posterior_figures.py", ["fig-sim-normal-1", "fig-sim-normal-2", "fig-sim-normal-3"]),
    ("application_figures.py", ["fig-real-digits-1", "fig-real-digits-2", "fig-real-clustering-1", "fig-real-clustering-2"]),
]
FIGURES = {stem for _, stems in GROUPS for stem in stems}


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def raster_names(path):
    return re.findall(r"\\includegraphics\[[^\]]*\]\{([^}]+)\}", path.read_text())


def check_published():
    manifest = json.loads((PUBLISHED / "manifest.json").read_text())["files"]
    for name, item in manifest.items():
        path = PUBLISHED / name
        if not path.is_file() or sha256(path) != item["sha256"]:
            raise RuntimeError("Published artifact is missing or changed: " + name)
    print(f"Verified {len(manifest)} published assets against their manuscript SHA-256 hashes.", flush=True)
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE / "build",
                        help="output directory (default: manuscript-figures/build)")
    parser.add_argument("--check-published", action="store_true",
                        help="verify copied publication assets and exit; no plotting dependencies needed")
    args = parser.parse_args()
    manifest = check_published()
    if args.check_published:
        return
    for module in ["matplotlib", "numpy", "pandas", "PIL"]:
        if importlib.util.find_spec(module) is None:
            raise RuntimeError(f"Missing Python dependency: {module}. See README.md.")
    if shutil.which("pdflatex") is None:
        raise RuntimeError("pdflatex was not found on PATH. Install/configure a TeX distribution as described in README.md.")
    for name in ["refresh/gaussian_resolution.csv", "gaussian/paired_results.csv",
                 "refresh/wasp_metrics.csv", "refresh/wasp_full_reference.csv"]:
        if not (REPO / "revision-work" / name).is_file():
            raise RuntimeError("Required saved result is missing: revision-work/" + name)
    output = args.output.resolve()
    if output == PUBLISHED.resolve():
        raise RuntimeError("Choose a build directory instead of overwriting the preserved published assets.")
    output.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["OTTO_FIGURE_OUTPUT"] = str(output)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    started = time.perf_counter()
    timings = {}
    for script, stems in GROUPS:
        step = time.perf_counter()
        subprocess.run([sys.executable, str(HERE / script)], cwd=HERE, env=env, check=True)
        timings[script] = round(time.perf_counter() - step, 3)
        for stem in stems:
            for suffix in [".pdf", ".pgf"]:
                path = output / (stem + suffix)
                if not path.is_file() or path.stat().st_size == 0:
                    raise RuntimeError("Figure output is missing: " + path.name)
            for name in raster_names(output / (stem + ".pgf")):
                if not (output / name).is_file():
                    raise RuntimeError("PGF raster dependency is missing: " + name)
    comparison = {}
    for stem in sorted(FIGURES):
        comparison[stem] = {
            "pgf_matches_published": sha256(output / (stem + ".pgf")) == manifest[stem + ".pgf"]["sha256"],
            "pdf_sha256": sha256(output / (stem + ".pdf")),
            "raster_dependencies_match_published": all(
                name in manifest and sha256(output / name) == manifest[name]["sha256"]
                for name in raster_names(output / (stem + ".pgf")))
        }
    report = {
        "figures_built": len(FIGURES),
        "seconds": round(time.perf_counter() - started, 3),
        "python_version": sys.version.split()[0],
        "packages": {name: importlib.metadata.version(name)
                     for name in ["matplotlib", "numpy", "pandas", "Pillow"]},
        "script_seconds": timings,
        "comparison": comparison,
        "pdf_hash_note": "PDF creation metadata can change hashes between otherwise identical builds. PGF and raster hashes compare the plotted contents.",
    }
    (output / "build-report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"Built {len(FIGURES)} figures in {report['seconds']:.1f} seconds.")
    print("Outputs and build-report.json: " + str(output))


if __name__ == "__main__":
    main()
