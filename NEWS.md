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
