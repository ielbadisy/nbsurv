# nbsurv

<!-- badges: start -->
[![R-CMD-check](https://github.com/ielbadisy/nbsurv/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ielbadisy/nbsurv/actions/workflows/R-CMD-check.yaml)
[![License: GPL-3](https://img.shields.io/badge/License-GPL--3-yellow.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

`nbsurv` fits conditional naive Bayes survival models for right-censored
time-to-event data. The package focuses on horizon-specific survival
prediction with inverse-probability of censoring weighting, and includes
tools for model fitting, prediction, evaluation, cross-validation,
hyper-parameter tuning, and permutation variable importance.

## Main workflow

```r
library(nbsurv)
library(survival)

lung$status <- as.integer(lung$status == 2)

fit <- nbsurv(
  Surv(time, status) ~ age + sex + ph.ecog,
  data = lung
)

surv_prob <- predict(fit, newdata = lung[1:5, ], times = c(100, 200, 400))
event_prob <- predict(fit, newdata = lung[1:5, ], times = c(100, 200, 400), type = "event")
```

## Evaluation

```r
metrics <- evaluate_nbsurv(
  fit,
  newdata = lung,
  times = c(100, 200, 400)
)

metrics
```

`evaluate_nbsurv()` currently returns:

- IPCW Brier score
- Concordance

for each requested prediction horizon.

## Cross-validation

```r
cv_fit <- cv_nbsurv(
  Surv(time, status) ~ age + sex + ph.ecog,
  data = lung,
  folds = 3,
  times = c(100, 200, 400),
  seed = 1
)

cv_fit$summary
```

## Hyper-parameter tuning

```r
grid <- data.frame(
  scale = c(TRUE, FALSE),
  laplace = c(1, 2),
  min_sd = c(0.05, 0.10)
)
grid$time_grid <- I(list(NULL, NULL))

tuned <- tune_nbsurv(
  Surv(time, status) ~ age + sex + ph.ecog,
  data = lung,
  param_grid = grid,
  folds = 3,
  times = c(100, 200, 400),
  seed = 1
)

tuned$best_params
```

## Notes

- `predict()` returns monotone survival curves by applying a cumulative minimum
  across increasing horizons.
- Continuous predictors are modeled with Gaussian class-conditional densities.
- Categorical predictors use Laplace-smoothed probabilities.
- The package exposes a single clean public API centered on `nbsurv()`,
  `predict()`, `evaluate_nbsurv()`, `cv_nbsurv()`, and `tune_nbsurv()`.
