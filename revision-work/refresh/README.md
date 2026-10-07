# Gaussian resolution and WASP revision records

This directory contains the exact input arrays, initial supports, saved fits,
and numerical summaries used by the current Gaussian-resolution figure,
WASP figures, and supplementary retained-draw sensitivity check. It is
self-contained: no parent project directory, original submission, or external
33 MB posterior archive is required. Plotting is handled by
`../../manuscript-figures/`.

From the repository root, verify the saved calculations with:

```sh
Rscript revision-work/refresh/refresh.R --check
Rscript revision-work/refresh/sensitivity.R --check
```

These use public T4transport transport calculations. The first checks all 36
Gaussian and 54 WASP per-run records, representative accuracy calculations,
and the full-data reference. The second recomputes all 12 sensitivity summaries
from saved supports. They do not optimize new barycenters. All scripts locate
their files relative to themselves and work from any current directory.

To rebuild summary files without fitting:

```sh
Rscript revision-work/refresh/refresh.R --assemble
Rscript revision-work/refresh/sensitivity.R --assemble
```

The resulting publication inputs retain their names:
`gaussian_resolution.csv`, `wasp_metrics.csv`, `wasp_full_reference.csv`, and
`draw_sensitivity.csv`. The full-data reference is bundled and independently
checked by `--check`; it does not need refitting.

## Full reruns and the recorded package build

Full optimization used a T4transport 0.1.9 build exposing the compiled routine
`cpp_free_bary_gradient_damped_init`. Its installed metadata did not retain a
source commit, so the version number alone is not a reliable feature check.
`helpers.R` explicitly checks for the routine and explains its requirement if
it is missing. The saved-fit checks, summary assembly, and publication plots
do not require this private routine. Stored initial supports remove any
dependency on reproducing the package's private initializer.

With the compiled routine available:

```sh
Rscript revision-work/refresh/refresh.R --check-native
Rscript revision-work/refresh/refresh.R --run
Rscript revision-work/refresh/sensitivity.R --run
```

`--check-native` additionally repeats one small published Gaussian fit and
compares its support and iteration count. `--run` computes only missing fits.
Use `--recompute` instead to replace all fits and runtimes in that study; this
is substantially more expensive. The recorded run used R 4.5.1 and
T4transport 0.1.9. New timing values will depend on the machine and build.

The stopping code is preserved: blocks of at most 50 native updates use
internal tolerance 1e-12. At each checkpoint, exact plans give a projection
M and trial update Z+ = (1-alpha)Z + alpha M. Acceptance of convergence requires
both `abs(F(Z+)-F(Z))/max(s²,F(Z)) <= 1e-8` and
`sqrt(sum_i v_i ||M_i-Z_i||² / s²) <= 1e-4` at the returned support, or an
iteration-limit status at 500 updates. Here s² is average within-input atomic
variance with a floor of 1e-12. Timing includes these diagnostic evaluations.
Atomic covariance divides by the total number of equally weighted support
atoms; it does not use the sample-covariance correction.

## Inputs and provenance

`gaussian_inputs.rds` contains three independent Gaussian problems (seeds
100001–100003), four 300-observation inputs per problem, population Gaussian
barycenter references, and 500 common oracle draws per problem (seed 300000+r).
The resolution grid is m=10,25,50,100,200,500 and alpha=1,0.5. Initial supports
use seeds 200000+1000r+m and are shared by both step sizes.

`wasp_inputs.rds` contains exactly the nine input configurations used in the
main WASP figures: the first three partition/sampling realizations for each
K=5,10,20, selected from the repository's reduced posterior archive at index
`(K-2)*50+r`. Each subset contains 100 retained draws. It also stores the 200
full-data reference draws, analytic posterior parameters, 1,000 common oracle
draws, and 2,000 test covariates. Seed 700000 generated the oracle sample and
then the test covariates. The grid uses m=50,100,200 and alpha=1,0.5, with shared
initialization seeds 600000+10000r+100K+m.

`sensitivity_inputs.rds` contains only the first three K=5 configurations from
the unreduced 200-draw archive, their fixed permutations and initial supports,
and the common oracle sample. Seed 950000+r initialized m=100; subset j was
permuted using seed 960000+100r+j. The resulting nested samples have
25,50,100,200 draws. The nested 100-draw input is not necessarily the reduced
100-draw input used in the main figures. All four fits within a realization
share the initial support from the complete 200-draw inputs.

`initial_supports.rds` stores the 18 Gaussian and 27 WASP starting supports.
They were captured from the recorded Gaussian initializer, which averages
input means and sample covariance diagonals, adds 1e-8 to the diagonal, and
draws the requested number of atoms. The sensitivity initial supports are
stored directly in its input bundle. `input-provenance.json` and
`input-checksums.sha256` document the extraction and bundled arrays.

The posterior realizations all come from one generated regression dataset
(n=1000, seed 496, standardized covariates, beta=(1,1), unit observation
variance, N(0,16I) prior). This port preserves the retained posterior measures;
it does not regenerate Stan chains or claim verified MCMC mixing. The interval
summaries describe conditional means and coefficient inclusion within that
fixed dataset. They exclude future observation noise and are not estimates
of predictive coverage or repeated-data calibration.

Saved files comprise 36 `gaussian_*.rds`, 54 `wasp_*.rds`, and 12
`sensitivity_*.rds` fit records. All published fits met their stopping checks.
The uncertainty bars are sample SD across three Gaussian problems or three
partition/sampling realizations, respectively. Sensitivity references use the
fit obtained from all 200 draws for that realization, without a global-optimality
claim. No old figure renderer, manuscript drafts, build logs, or review files
are included here.
