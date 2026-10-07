# Candidate-restricted MNIST point-cloud medoid

This folder contains the same-representation prototype check reported in the second-revision manuscript and supplement. It preserves the three completed split results, numerical checkpoints, original timing metadata, and input serialization. It does not contain manuscript or reviewer-response drafts.

## Design and reported result

The 5,000 `T4transport::digits` images are converted with the repository's existing Otsu helper: 256 bins, foreground coordinates in `[-1,1]^2`, and uniform atom masses. Split seeds 800001–800003 give 100 training and 250 test images per digit. The first ten class members in each shuffled training order form a prespecified candidate set. The selected representative minimizes mean squared exact Wasserstein distance to all 100 class-training images. Classification uses the closest of the ten selected class representatives under exact W2, with class-order tie breaking.

Mean accuracy is **0.5968**, with sample SD **0.0281623863** across the three splits. These describe split variability conditional on this fixed corpus. The representatives have their original foreground atom counts; this is a same-representation, one-prototype-per-class comparison, **not a matched atom-budget comparison or a full-search medoid**.

## Verify the saved results

Run these commands from the repository root. Scripts resolve repository files relative to their own locations, so they also work from other directories.

```sh
python3 revision-work/mnist/03_report.py --check
Rscript revision-work/mnist/04_verify.R
python3 revision-work/mnist/02_medoid.py --validate-only
```

The first command uses only Python's standard library and independently reconstructs classification metrics from the 7,500 saved predictions, then checks split metrics, means, and SDs. The R check uses `jsonlite` and the existing scoring helper to verify split balance, train/test separation, candidate membership, selected candidates, and classification metrics. The final command requires NumPy, SciPy, and POT. It checks input hashes and solves only 20 exact-OT problems against the saved R reference; it does not reclassify the test sets or replace historical runtime metadata.

## Reproduce the experiment

Dependencies are R packages `T4transport` and `jsonlite`, and Python packages `numpy`, `scipy`, and `POT`. Recorded versions are in `r_metadata.json` and `python_metadata.json`; the latter records POT 0.9.7.post1, NumPy 2.4.6, and SciPy 1.17.1.

The bundled `point_clouds.csv.gz` is a losslessly compressed version of the original CSV. Its decompressed bytes match the historical SHA-256 in `python_metadata.json`. `images.csv` and `splits.json` preserve image identifiers, labels, and R-generated split indices, avoiding cross-language RNG assumptions.

To replay checkpointed processing and rebuild the numerical summaries:

```sh
python3 revision-work/mnist/02_medoid.py
python3 revision-work/mnist/03_report.py
```

Existing candidate costs and classification distances are reused. `execution_metadata.json` records the current invocation without replacing the original `python_metadata.json` timing provenance. Use `--rep 1`, `--rep 2`, or `--rep 3` to select a split.

For a fresh, separately timed computation, use a new directory. The environment variable must be set for each command or exported as below:

```sh
export MNIST_WORK_DIR="$(mktemp -d)"
Rscript revision-work/mnist/01_prepare.R
python3 revision-work/mnist/02_medoid.py
python3 revision-work/mnist/03_report.py
Rscript revision-work/mnist/04_verify.R
```

Preparation uses [`code-real-digits/common/rev_utils.R`](../../code-real-digits/common/rev_utils.R), recording its checksum. By default it preserves existing prepared inputs; `--force` explicitly regenerates them. Use a fresh work directory whenever changing inputs so that old checkpoints cannot be mistaken for new results.

POT uses float64 network simplex, one solver thread, and at most 1,000,000 iterations. Construction and classification times exclude preprocessing, file access, and independent validation. Retained timings describe the original Python execution and are not a controlled speed comparison with historical R barycenter runs.

## Files and limits

- `rep1/`–`rep3/`: candidate objectives, selected image IDs, candidate-training costs, prototype-test distances, predictions, confusion matrices, and metrics.
- `summary.json`, `split_metrics.csv`: reported means, sample SDs, and individual split results.
- `distance_validation.json`, `r_reference_distances.csv`, `verification.json`: recorded solver and split/metric checks.
- `provenance.json`: hashes of the retained source evidence before portable packaging. The point-cloud CSV is compressed, and scripts are adapted; no numerical study was rerun during packaging.

Historical barycenter core classification and centroid outputs are not included in the available legacy digit directory (which contains smoke examples). This folder therefore does not claim to reconstruct paired medoid-minus-barycenter differences or certify convergence of the old 50-iteration barycenter runs. The old experiment scripts remain under `code-real-digits/`.
