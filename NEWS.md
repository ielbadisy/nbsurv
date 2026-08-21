# nbsurv 0.4.0

## New features

* New `cov_structure = c("diagonal", "full")` argument to `nbsurv()`/
  `cv_nbsurv()` (default `"diagonal"`, fully backward-compatible). With
  `"full"`, continuous predictors are modeled jointly as a single
  multivariate Gaussian per class (survivor/event), relaxing naive Bayes's
  conditional-independence assumption between them, instead of the product
  of independent univariate Gaussians. The covariance matrix is shrunk
  toward its diagonal by the new `shrinkage` argument (default `0.2`) for
  numerical stability at small sample sizes. Categorical predictors are
  unaffected and remain conditionally independent either way. With fewer
  than 2 continuous predictors, `"full"` is numerically identical to
  `"diagonal"` (verified by regression test).
* **Empirical validation** (see `tests/testthat/test-cov-structure.R`): on
  a synthetic case with two continuous predictors correlated at r=0.85,
  `cov_structure = "full"` gives a consistent ~5-6% relative improvement in
  held-out IPCW Brier score over `"diagonal"` at every horizon tested. Its
  effect on **concordance is small and inconsistent** even under strong
  correlation - concordance depends only on the *ranking* of risk scores,
  which a monotonic transform of the naive product often preserves even
  when the underlying probabilities are miscalibrated; Brier score, which
  scores the actual probability values, is the metric this correction
  reliably improves. On the real `lung` dataset (`age`/`ph.ecog` correlated
  at r=0.31, n=167), the improvement was smaller and mixed (Brier improved
  at 2 of 3 horizons tested, concordance improved at 2 of 3) - `"full"` is
  a genuine, validated improvement when continuous predictors are
  correlated, not a guaranteed win on every dataset, and did not close the
  gap to `coxph` on `lung` specifically (see README benchmark).

# nbsurv 0.3.2

## Bug fixes (scientific validity)

* **`evaluate_nbsurv()`/`cv_nbsurv()`/`tune_nbsurv()` reported the
  *complement* of concordance (`1 - true concordance`), not concordance
  itself.** `concordance_at_horizon()` calls
  `survival::concordance(Surv(time, status) ~ risk)`, whose formula
  interface defaults to the *opposite* direction convention from a risk
  score (it expects higher `risk` to mean *longer* survival unless told
  otherwise). Since every internal caller passes a genuine risk-direction
  score (higher = more likely to fail), the returned value was
  systematically inverted - a model with a true C-index of 0.63 was
  reported as 0.37, making a genuinely discriminating model look
  *worse than random*. Fixed by adding `reverse = TRUE` to the
  `survival::concordance()` call. Verified against `survival::concordance()`
  on both a `coxph` linear predictor and an `nbsurv` fit (matches to full
  precision either way), and a new regression test on a clearly-separable
  synthetic case that would fail loudly if the direction ever flips back.
  Brier scores, predictions, and every other metric were unaffected -
  this only touched the concordance calculation.

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
