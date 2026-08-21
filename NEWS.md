# nbsurv 0.3.1

## Bug fixes

* **`plot.nbsurv()` silently double-scaled continuous predictors before
  plotting, distorting every survival curve it drew.** `nbsurv()` only ever
  stored the *scaled* training data on the fit object; `plot.nbsurv()`
  passed that already-scaled data into `predict()`, whose `prepare_newdata()`
  step scales its input again using the same stored means/sds. The result
  was visually "chunky" curves with implausible flat segments and near-zero
  collapses for some rows, especially past t=400 in the README example -
  not real model behavior, an artifact of predicting from
  twice-standardized covariates. Fixed by storing the raw (unscaled)
  training data separately (`fit$training_data_raw`) and having
  `plot.nbsurv()` use that instead. `predict()`/`evaluate_nbsurv()`/
  `cv_nbsurv()`/`tune_nbsurv()`/`varimp_nbsurv()` were never affected - they
  operate on user-supplied `newdata`, which is only ever scaled once. New
  regression test locks `plot.nbsurv()`'s curves to match `predict()`
  called directly on the same raw rows.

# nbsurv 0.3.0

## New features

* New `varimp_nbsurv()`: permutation-based variable importance. For each
  predictor, its values are independently shuffled (`n_repeats` times),
  `evaluate_nbsurv()` is recomputed on the permuted data, and importance is
  the resulting degradation in the chosen metric (`"brier"` or
  `"concordance"`), averaged over the requested horizons and repeats.
  Returns a `data.frame` (class `varimp_nbsurv`) sorted by decreasing
  importance, with its own `print()` method.
* New base-graphics `plot()` methods: `plot.varimp_nbsurv()` (horizontal
  bar chart of importance), `plot.cv_nbsurv()` (mean cross-validation
  metric vs. horizon), and `plot.tune_nbsurv()` (mean tuning metric per
  candidate, best candidate highlighted).

# nbsurv 0.2.0

## New features

* `evaluate_nbsurv()` now computes the **integrated Brier score** (IBS) via the
  trapezoidal rule in addition to horizon-specific Brier scores and concordance.
* New `calibration_plot_nbsurv()` function for visual calibration assessment at
  a single prediction horizon.

## Bug fixes and improvements

* Continuous predictor standard deviations are now bounded away from zero by
  `min_sd` in both the survivor and event classes, preventing numerical collapse
  in degenerate training subsets.
* `predict.nbsurv()` now enforces monotone survival curves via cumulative minimum
  across increasing horizons.
* Categorical likelihoods now use Laplace smoothing controlled by the `laplace`
  argument, with IPCW-weighted event-class counts.

# nbsurv 0.1.0

* Initial release (prototype; class `CNB`, single-horizon prediction).
