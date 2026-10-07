# Current manuscript figures

This directory contains the ten figures used by the current manuscript, their plotting scripts, saved plotting inputs, and the raster layers required by their PGF files. The figure build reads saved results. It does not rerun any statistical experiment or require access to files outside this repository.

## Build all ten figures

From the repository root:

```sh
python manuscript-figures/build_figures.py
```

The command writes PDF and PGF files to `manuscript-figures/build/`, together with their required raster layers and `build-report.json`. It leaves the preserved manuscript figures in `published/` unchanged. The output directory can be selected with `--output PATH`.

The saved-data build requires Python with `matplotlib`, `numpy`, `pandas`, and `Pillow`, plus `pdflatex` on `PATH`. The TeX installation needs PGF, `amsmath`, `amssymb`, `underscore`, and `scrextend`. No R execution or network access is needed for this build. Install the repository's Python dependencies using the root setup instructions.

The build reads the following experiment summaries from the repository:

- `revision-work/refresh/gaussian_resolution.csv`
- `revision-work/gaussian/paired_results.csv`
- `revision-work/refresh/wasp_metrics.csv`
- `revision-work/refresh/wasp_full_reference.csv`

The remaining plotting inputs are supplied in `data/` and `legacy-panels/`. See the source limitations below before interpreting a successful figure build as a rerun of an experiment.

## Preserved publication assets

`published/` contains the exact ten PDF and PGF files used in the manuscript, with every raster dependency referenced by the PGF sources. `published/manifest.json` records their SHA-256 hashes. Verify this preserved copy without plotting dependencies or TeX:

```sh
python manuscript-figures/build_figures.py --check-published
```

The PDF files are ready for inclusion with `\includegraphics`. Keep each PGF file beside its `*-img*.png` dependencies. When including a PGF from another directory, use a LaTeX import mechanism that resolves those dependencies relative to the PGF file.

| Manuscript figure | File stem | Generator | Plot source |
|---|---|---|---|
| 1 | `fig-sim-gauss-2` | `gaussian_figures.py` | 36 saved Gaussian fits |
| 2 | `fig-sim-gauss-3` | `runtime_figure.py` | Preserved runtime plot interiors |
| 3 | `fig-sim-gauss-4` | `gaussian_figures.py` | 150 saved paired solver fits |
| 4 | `fig-sim-normal-1` | `posterior_figures.py` | Saved posterior density grids |
| 5 | `fig-sim-normal-2` | `posterior_figures.py` | 54 saved posterior aggregation fits |
| 6 | `fig-sim-normal-3` | `posterior_figures.py` | Same 54 fits and saved full-data reference |
| 7 | `fig-real-digits-1` | `application_figures.py` | Original digit image and Otsu threshold |
| 8 | `fig-real-digits-2` | `application_figures.py` | Preserved particle and density plot interiors |
| 9 | `fig-real-clustering-1` | `application_figures.py` | Saved PBMC/FASHION PCA coordinates and preserved NEWS plot interior |
| 10 | `fig-real-clustering-2` | `application_figures.py` | Preserved clustering curves and error bars |

## Source limitations and graphical fidelity

Figures 1, 3–7 and panels 9A/9C are plotted from saved numerical inputs. Their means, error bars, density grids, coordinates, and transformations are unchanged from the manuscript figures. Figures 1 and 3 use means and sample standard deviations. Figure 1 contains both step sizes, with nested open and filled markers exposing their near overlap. Figures 5 and 6 use small horizontal display offsets to separate coincident summaries. Those offsets do not change the means or standard deviations.

Complete historical inputs were not available for Figure 2, Figure 8, Figure 10, or panel 9B after checking both original and revised experiment folders, the submitted code archive, and the saved experiment objects. Their exact original plot interiors are retained in `legacy-panels/`. The scripts replace axes, headings and legends with TeX typography. Figure 8 density shading is converted to grayscale. No missing coordinates, numerical results, replicates, or error bars are estimated from these images.

In particular, the available runtime smoke rows cannot reproduce the complete runtime sweep. The available digit smoke fits and earlier 50-image fits cannot replace the displayed 100-image prototypes. The saved NEWS UMAP coordinates cannot replace its missing processed PCA features. The incomplete saved clustering runs cannot reproduce all of Figure 10. The successful build therefore reproduces these figure graphics, not the unavailable historical experiments.

All ten figures use shared `rcParams` in `plot_style.py`, the manuscript's Computer Modern font, and their final display widths. Labels use 12-point type, most tick and legend text uses 10-point type, and dense application axes use 9-point ticks. The multiple-series plots retain line patterns, marker shapes, or open and filled symbols so their interpretation does not depend on color. Figure 9 displays dataset geometry rather than identifying reference classes in a monochrome printout. Its dense regenerated scatter layers use 600 dpi rasterization to avoid pdfTeX memory limits. Axes and typography remain vector elements.

## Refresh optional plotting exports

The cached exports in `data/` make the standard build independent of R. The scripts below can regenerate them from the repository's saved R objects and installed R packages. Run them from any working directory:

```sh
Rscript manuscript-figures/export_digit_image.R
Rscript manuscript-figures/export_application_pca.R
Rscript manuscript-figures/export_posterior.R
```

`export_digit_image.R` needs `T4transport` and reproduces the original 256-bin Otsu rule for the first digit-8 image. `export_application_pca.R` uses the stored PBMC/FASHION processed matrices, the original sampling seeds, and R's `prcomp`. `export_posterior.R` needs `MASS` and `mvtnorm`, the saved posterior arrays, and `revision-work/refresh/wasp_1_5_100_1.rds`. These commands accept an optional output directory to compare a regenerated export with the committed cache before replacing it.

`audit_gaussian_supports.R` optionally verifies the differences between full-step and half-step Gaussian supports using the saved RDS fits and `T4transport`. Its committed result is `data/gaussian_support_audit.csv`.

`data/manuscript-tables/` preserves all eleven CSV summaries accompanying the current manuscript. Gaussian and WASP summaries reflect the refreshed revision results. Digit and clustering summaries retain the historical application results used in the manuscript. These application summaries are not substitutes for the missing individual runs or full plotting inputs. The figure scripts use the run-level CSVs listed above whenever those are available.

## Reproduction checks

The build checks that all ten PDF/PGF pairs and their raster dependencies exist. `build-report.json` records the runtime, package versions, and whether regenerated PGF files and raster dependencies match the preserved publication assets. PDF creation metadata can change their binary hashes between otherwise identical builds. Differences caused by plotting-library or TeX versions should be reviewed in the rendered figures before replacing publication assets.

The repository port was verified with Python 3.11.4, Matplotlib 3.11.1, NumPy 2.4.6, pandas 3.0.5, Pillow 12.3.0, and pdfTeX from TeX Live 2023. All ten figures built in about 73 seconds on the verification machine. Every regenerated PGF and required raster asset matched the corresponding publication asset byte for byte. Rendered PDFs also matched at 120 dpi, and all fourteen optional R plotting exports matched their committed caches. `reproduction-check.json` records these checks. `input-manifest.json` records hashes for the supplied plotting data, table summaries and preserved plot layers.
