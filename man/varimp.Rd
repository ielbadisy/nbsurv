\name{varimp_nbsurv}
\alias{varimp_nbsurv}
\alias{print.varimp_nbsurv}
\title{Permutation Variable Importance for nbsurv Models}
\usage{
varimp_nbsurv(
  object,
  newdata,
  times,
  metric = c("brier", "concordance"),
  n_repeats = 10,
  seed = NULL
)
}
\arguments{
  \item{object}{A fitted \code{"nbsurv"} model.}
  \item{newdata}{Evaluation data containing the same outcome and predictors as the training formula.}
  \item{times}{Prediction horizons averaged over when computing the metric.}
  \item{metric}{The metric to degrade under permutation. One of \code{"brier"} (higher = worse) or
    \code{"concordance"} (higher = better).}
  \item{n_repeats}{Number of independent permutations per feature; importance is the mean degradation
    across repeats.}
  \item{seed}{Optional random seed for reproducible permutations.}
}
\value{
A data frame (class \code{"varimp_nbsurv"}) with one row per predictor, sorted by decreasing
importance, with columns \code{feature}, \code{baseline} (the unpermuted metric, averaged over
\code{times}), \code{permuted} (the metric after permuting that feature, averaged over
\code{times} and \code{n_repeats}), and \code{importance} (the degradation: \code{permuted -
baseline} for \code{"brier"}, \code{baseline - permuted} for \code{"concordance"}, so higher is
always more important).
}
\description{
Computes permutation-based variable importance: for each predictor, its values are independently
shuffled \code{n_repeats} times, \code{\link{evaluate_nbsurv}} is recomputed on the permuted data,
and importance is the resulting degradation in the chosen metric, averaged over \code{times} and
repeats. Predictors the model relies on more heavily degrade performance more when permuted.
}
\examples{
library(survival)

lung <- stats::na.omit(lung)
lung$status <- as.integer(lung$status == 2)
fit <- nbsurv(Surv(time, status) ~ age + sex + ph.ecog, data = lung)
varimp_nbsurv(fit, newdata = lung, times = c(100, 200, 400), n_repeats = 5, seed = 1)
}
