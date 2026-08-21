# nbsurv 0.5.0

## New features: kernel-smoothed horizons (`time_smooth`)

* New `time_smooth = TRUE` / `bandwidth` arguments to `nbsurv()`/
  `cv_nbsurv()` (default `FALSE`, fully backward-compatible). The real
  remaining inefficiency identified after `cov_structure = "full"` (0.4.0):
  every horizon was fit *completely independently*, discarding the fact
  that class-conditional statistics should vary smoothly with time -
  costly right where one of the two horizon-defined classes (already
  failed vs. known to survive) is smallest, typically near the earliest or
  latest requested horizon. `time_smooth = TRUE` precomputes class-
  conditional statistics once at every point of `time_grid` at fit time,
  then combines them at `predict()` time via Nadaraya-Watson kernel
  regression across horizons (Gaussian kernel, width `bandwidth`) -
  borrowing strength from neighboring horizons while still allowing
  genuinely time-varying effects, unlike a proportional-hazards model
  that shares one coefficient across all of time. Because the combination
  is a convex (weight-normalized) average of already-valid per-grid-point
  statistics, smoothed covariance matrices stay positive semi-definite and
  smoothed categorical probabilities stay valid automatically - no PSD
  repair or renormalization needed. `bandwidth` is now a tunable column in
  `tune_nbsurv()`'s `param_grid`.
* **Bug fix found and fixed during implementation**: the first working
  version returned `NA` for every smoothed prediction. `0 * NA` is still
  `NA` in R, so zero-weighting an invalid (all-`NA`) grid point before
  summing did not neutralize it; `colSums()` propagated the `NA`. Fixed by
  zeroing the invalid rows' *values*, not just their weights, before
  summing. Regression test locks this in
  (`tests/testthat/test-time-smooth.R`).
* **Empirical validation** against `coxph` on `lung` (20 held-out 70/30
  splits, `age + sex + ph.ecog`): a genuine, principled bias-variance
  tradeoff, confirmed rather than assumed. Concordance improves as
  `bandwidth` widens (0.60 at the default heuristic up to ~0.615-0.619 as
  bandwidth grows toward pooling most of the time range), closing roughly
  half to two-thirds of the gap to `coxph`'s 0.629 - but IPCW Brier score
  at the most extreme horizon tested (t=400) *worsens* monotonically as
  bandwidth widens (0.240 to 0.254), the expected cost of over-smoothing
  away genuine horizon-specific calibration. No single bandwidth wins on
  both metrics at every horizon; the default heuristic (one quarter of the
  fitted time range) is a reasonable balance, not a universally optimal
  choice - **use `tune_nbsurv()`'s `bandwidth` grid column to tune it for
  a specific dataset**, exactly as any other hyperparameter.
* Even at the best bandwidths tried, the `lung` benchmark gap to `coxph`
  narrows but does not fully close - reported honestly in the README
  rather than oversold. `age`/`ph.ecog`'s covariate effects are close to
  log-linear and non-interacting on this dataset, which is exactly what a
  correctly-specified Cox model is built to exploit; `nbsurv`'s value
  proposition is flexibility (genuinely time-varying, non-proportional
  effects a fixed-coefficient Cox model cannot represent), not raw
  accuracy dominance on every dataset.
* **No C++ was added.** Benchmarked at n=5000 (3270 distinct event times):
  `time_smooth = TRUE` fit time is ~2.5s vs ~0.03s without it - a real but
  tolerable one-time cost at realistic survival-analysis sample sizes,
  with predict() itself remaining fast (~17ms for 500 predictions) since
  it only does cheap kernel-weighted averaging over the precomputed grid.
  The existing `time_grid` argument already lets users pass a coarser
  custom grid to cut this cost further if needed. Rcpp would be
  unjustified complexity for a bottleneck that doesn't clearly exist at
  the scale this package targets.

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
