# Experiment inputs

The repository includes these original inputs:

- `normal-simplified-reduced.RData`: stored subset posterior arrays and the full-data MCMC reference for WASP. The reported runs use three partition/sampling realizations of one simulated dataset, not three independently generated datasets.
- `processed_pbmc.RData`: 2,638 quality-controlled cells represented in 20 dimensions, with nine reference labels.
- `processed_fmnist.RData`: 10,000 Fashion-MNIST test images represented in 20 dimensions, with ten reference labels.

`processed_news.RData` is absent. Complete recomputation of the original three-dataset clustering study requires this missing matrix. Saved figure layers and reported summaries are preserved without claiming that they replace the missing observations.

Gaussian inputs are synthetic. The current paired Gaussian comparison bundles its exact serialized inputs under `../revision-work/gaussian/inputs`. The MNIST medoid comparison bundles the thresholded point clouds and split records under `../revision-work/mnist`, with their source and validation described there. The original MNIST notebooks load the `digits` data from `T4transport`.

The current posterior draw-count sensitivity comparison has a compact input record under `../revision-work/refresh`, so it does not require an external full posterior archive. See that directory's README for its exact filename and contents.
