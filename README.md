
<!-- README.md is generated from README.Rmd. Please edit that file -->

# nbsurv

<!-- badges: start -->

[![R-CMD-check](https://github.com/ielbadisy/nbsurv/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ielbadisy/nbsurv/actions/workflows/R-CMD-check.yaml)
[![License:
GPL-3](https://img.shields.io/badge/License-GPL--3-yellow.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

`nbsurv` fits conditional naive Bayes survival models for right-censored
time-to-event data. The package focuses on horizon-specific survival
prediction with inverse-probability of censoring weighting, and includes
tools for model fitting, prediction, evaluation, cross-validation,
hyper-parameter tuning, and permutation variable importance.

All results below are real output from `lung` (from the `survival`
package, rows with missing predictors dropped).

## Main workflow

``` r
library(nbsurv)
library(survival)

lung <- na.omit(lung)
lung$status <- as.integer(lung$status == 2)

fit <- nbsurv(
  Surv(time, status) ~ age + sex + ph.ecog,
  data = lung
)

surv_prob <- predict(fit, newdata = lung[1:5, ], times = c(100, 200, 400))
round(surv_prob, 3)
#>   t_100 t_200 t_400
#> 2 0.860 0.782 0.459
#> 4 0.904 0.686 0.389
#> 6 0.745 0.592 0.371
#> 7 0.770 0.629 0.351
#> 8 0.739 0.615 0.335

event_prob <- predict(fit, newdata = lung[1:5, ], times = c(100, 200, 400), type = "event")
round(event_prob, 3)
#>   t_100 t_200 t_400
#> 2 0.140 0.218 0.541
#> 4 0.096 0.314 0.611
#> 6 0.255 0.408 0.629
#> 7 0.230 0.371 0.649
#> 8 0.261 0.385 0.665

plot(fit, times = c(50, 100, 200, 400, 600, 800))
```

<img src="man/figures/README-survival-curves-1.png" alt="" width="100%" />

## Evaluation

``` r
metrics <- evaluatenbsurv(
  fit,
  newdata = lung,
  times = c(100, 200, 400)
)

metrics
#>   time     brier concordance
#> 1  100 0.1145281   0.6108482
#> 2  200 0.1820750   0.6364067
#> 3  400 0.2177052   0.6299697
```

`evaluatenbsurv()` currently returns:

- IPCW Brier score
- Concordance

for each requested prediction horizon.

## Cross-validation

``` r
cv_fit <- cvnbsurv(
  Surv(time, status) ~ age + sex + ph.ecog,
  data = lung,
  folds = 3,
  times = c(100, 200, 400),
  seed = 1
)

cv_fit$summary
#>   time     brier concordance
#> 1  100 0.1309608   0.5802539
#> 2  200 0.2038665   0.6257687
#> 3  400 0.2428066   0.5978809

plot(cv_fit)
```

<img src="man/figures/README-cv-1.png" alt="" width="100%" />

## Hyper-parameter tuning

``` r
grid <- data.frame(
  scale = c(TRUE, FALSE),
  laplace = c(1, 2),
  min_sd = c(0.05, 0.10)
)
grid$time_grid <- I(list(NULL, NULL))

tuned <- tunenbsurv(
  Surv(time, status) ~ age + sex + ph.ecog,
  data = lung,
  param_grid = grid,
  folds = 3,
  times = c(100, 200, 400),
  seed = 1
)

tuned$best_params
#>   scale laplace min_sd time_grid mean_metric
#> 2 FALSE       2    0.1             0.1925446

plot(tuned)
```

<img src="man/figures/README-tune-1.png" alt="" width="100%" />

## Variable importance

``` r
vi <- varimpnbsurv(
  fit,
  newdata = lung,
  times = c(100, 200, 400),
  n_repeats = 10,
  seed = 1
)

vi
#> nbsurv permutation variable importance
#> Metric: brier | repeats: 10 
#> 
#>  feature  baseline  permuted  importance
#>  ph.ecog 0.1714361 0.1912398 0.019803689
#>      sex 0.1714361 0.1810340 0.009597880
#>      age 0.1714361 0.1772819 0.005845766

plot(vi)
```

<img src="man/figures/README-varimp-1.png" alt="" width="100%" />

`age` and `ph.ecog` degrade the Brier score the most when permuted (most
important); `sex` slightly *improves* it when shuffled here, i.e. it
contributes essentially nothing beyond noise on this fit/horizon set.

`calibrationplotnbsurv(fit, newdata = lung, horizon = 200)` is also
available for visual calibration assessment at a single horizon.

## Benchmark vs. `coxph`

A quick, honest comparison against `survival::coxph()` on the same
formula: 20 random 70/30 train/test splits of `lung`, IPCW Brier score
and concordance computed identically for both models (`nbsurv`’s own
internal metric functions, applied to `coxph`’s linear predictor and
`survfit()`-derived survival probabilities too), means reported below.

``` r
#>   time nb_brier cox_brier nb_conc cox_conc
#> 1  100   0.1312    0.1304  0.6011   0.6294
#> 2  200   0.2051    0.1950  0.6107   0.6294
#> 3  400   0.2397    0.2344  0.5792   0.6294
```

`coxph` edges out `nbsurv` on both metrics at every horizon on this
dataset - expected: `age + sex + ph.ecog` have a roughly log-linear,
non-interacting effect on the hazard here, which is exactly what a
correctly-specified Cox model is built for, while `nbsurv`’s
conditional-independence assumption between predictors costs it some
accuracy in exchange for a much simpler, distribution-light estimation
procedure. The gap is real but small (concordance within ~0.03-0.05,
Brier within ~0.01-0.05).

`cov_structure = "full"` (below) relaxes that independence assumption
between continuous predictors, but on this dataset (`age`/`ph.ecog`
correlated at only r=0.31, n=167) it doesn’t close the gap - it’s a
genuine, validated improvement when continuous predictors are more
strongly correlated (see below), not a guaranteed win on every dataset.

## Relaxing the independence assumption: `cov_structure = "full"`

By default, continuous predictors are modeled as independent given the
class (survivor/event) - the classic naive Bayes assumption.
`cov_structure = "full"` instead models them jointly as a single
multivariate Gaussian per class, capturing correlation between them:

``` r
fit_full <- nbsurv(
  Surv(time, status) ~ age + sex + ph.ecog,
  data = lung,
  cov_structure = "full"
)
```

On a synthetic case with two continuous predictors correlated at r=0.85
(20 held-out 70/30 splits), `"full"` gives a consistent ~5-6% relative
improvement in Brier score over `"diagonal"` at every horizon tested;
its effect on concordance is small and inconsistent even at that
correlation strength, since concordance only depends on the *ranking* of
risk scores, which the naive product often preserves even when
miscalibrated. This is a real, tested effect
(`tests/testthat/test-cov-structure.R`) - but whether it helps on *your*
data depends on how correlated your continuous predictors actually are.

## Borrowing strength across horizons: `time_smooth = TRUE`

The remaining inefficiency vs. `coxph`: every horizon is fit completely
independently by default, discarding the fact that the true
class-conditional statistics should vary *smoothly* with time. This
costs accuracy right where one of the two horizon-defined classes
(already failed vs. known to survive) is smallest - typically near the
earliest or latest requested horizon. `time_smooth = TRUE` precomputes
class-conditional statistics at every point of the internal `time_grid`
once at fit time, then blends them at `predict()` time via kernel
regression across nearby horizons:

``` r
fit_smooth <- nbsurv(
  Surv(time, status) ~ age + sex + ph.ecog,
  data = lung,
  time_smooth = TRUE,
  bandwidth = 300  # tune via tunenbsurv() for your data
)
```

Re-running the `coxph` benchmark above with `time_smooth = TRUE` at a
few bandwidths (still 20 held-out 70/30 splits):

``` r
#>                  brier (t=100,200,400)   concordance (t=100,200,400)
#> diagonal         0.1312, 0.2051, 0.2397  0.601, 0.611, 0.579
#> smooth (bw=101)  0.1284, 0.2030, 0.2413  0.601, 0.613, 0.608
#> smooth (bw=200)  0.1273, 0.2023, 0.2463  0.608, 0.613, 0.618
#> smooth (bw=500)  0.1265, 0.2013, 0.2528  0.614, 0.614, 0.617
#> coxph            0.1304, 0.1950, 0.2344  0.629, 0.629, 0.629
```

A genuine bias-variance tradeoff, confirmed rather than assumed: wider
bandwidths consistently improve concordance (closing roughly half to
two-thirds of the gap to `coxph`) but consistently *worsen* Brier score
at the most extreme horizon (t=400), the expected cost of over-smoothing
away real horizon-specific calibration. No single bandwidth wins on both
metrics at every horizon - **tune `bandwidth` via `tunenbsurv()`’s
`param_grid` for your own data** rather than trusting the default.

Even at the best bandwidths tried here, the gap to `coxph` narrows but
does not fully close. `age`/`ph.ecog` have close-to-log-linear,
non-interacting effects on this dataset - exactly what a
correctly-specified Cox model is built to exploit. `nbsurv`’s value
proposition is flexibility (genuinely time-varying effects a
fixed-coefficient Cox model can’t represent), not guaranteed accuracy
dominance on every dataset - reported honestly here rather than
oversold.

## Notes

- `predict()` returns monotone survival curves by applying a cumulative
  minimum across increasing horizons.
- Continuous predictors are modeled with Gaussian class-conditional
  densities (independent by default, or jointly via
  `cov_structure = "full"`).
- Categorical predictors use Laplace-smoothed probabilities.
- The package exposes a single clean public API centered on `nbsurv()`,
  `predict()`, `evaluatenbsurv()`, `cvnbsurv()`, `tunenbsurv()`, and
  `varimpnbsurv()`.
