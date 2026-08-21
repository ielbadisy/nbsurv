\name{cv_nbsurv}
\alias{cv_nbsurv}
\alias{evaluate_nbsurv}
\alias{tune_nbsurv}
\title{Evaluate, Cross-Validate, and Tune nbsurv Models}
\usage{
cv_nbsurv(
  formula,
  data,
  folds = 5,
  times = NULL,
  seed = NULL,
  scale = TRUE,
  laplace = 1,
  min_sd = 0.05,
  time_grid = NULL,
  eps = 1e-06,
  cov_structure = c("diagonal", "full"),
  shrinkage = 0.2,
  time_smooth = FALSE,
  bandwidth = NULL
)

evaluate_nbsurv(object, newdata, times, metrics = c("brier", "concordance"), ibs = FALSE)

tune_nbsurv(
  formula,
  data,
  param_grid,
  folds = 5,
  times = NULL,
  metric = c("brier", "concordance"),
  maximize = NULL,
  seed = NULL,
  eps = 1e-06
)
}
\arguments{
  \item{formula}{A model formula with a \code{survival::Surv()} response.}
  \item{data}{A data frame containing outcome and predictors.}
  \item{folds}{Number of cross-validation folds.}
  \item{times}{Prediction horizons to evaluate. If \code{NULL}, event-time quantiles are used.}
  \item{seed}{Optional random seed for reproducible fold assignment.}
  \item{scale}{Logical; if \code{TRUE}, continuous predictors are standardized.}
  \item{laplace}{Laplace smoothing constant for categorical likelihoods.}
  \item{min_sd}{Lower bound used for Gaussian standard deviations.}
  \item{time_grid}{Optional prediction grid stored with each fitted model.}
  \item{eps}{Small numerical constant used in probability clipping and IPCW calculations.}
  \item{cov_structure}{Passed to each fold's \code{\link{nbsurv}} fit; see there.}
  \item{shrinkage}{Passed to each fold's \code{\link{nbsurv}} fit; see there.}
  \item{time_smooth}{Passed to each fold's \code{\link{nbsurv}} fit; see there.}
  \item{bandwidth}{Passed to each fold's \code{\link{nbsurv}} fit; see there.}
  \item{object}{A fitted \code{"nbsurv"} model.}
  \item{newdata}{Evaluation data containing the same outcome and predictors as the training formula.}
  \item{metrics}{Performance metrics to compute. One or both of \code{"brier"} and \code{"concordance"}.}
  \item{ibs}{Logical; if \code{TRUE} and \code{"brier"} is in \code{metrics}, the integrated Brier score
    (trapezoidal rule over \code{times}) is attached as \code{attr(result, "ibs")}.}
  \item{param_grid}{A data frame with columns \code{scale}, \code{laplace}, \code{min_sd}, and \code{time_grid}.}
  \item{metric}{The metric used to choose the best hyper-parameter setting.}
  \item{maximize}{Logical; if \code{TRUE}, larger metric values are preferred.}
}
\value{
\code{evaluate_nbsurv()} returns a data frame of horizon-specific metrics. When
\code{ibs = TRUE} the integrated Brier score is attached via \code{attr(result, "ibs")}.
\code{cv_nbsurv()} returns a list with fold-level and mean cross-validation summaries.
\code{tune_nbsurv()} returns a list with all candidate results and the best row.
}
\description{
Provides resampling-based utilities for evaluating, cross-validating, and
tuning \code{nbsurv} survival models using horizon-specific Brier score,
integrated Brier score (IBS), and concordance.

The integrated Brier score summarises calibration over a range of horizons:
\deqn{
  \mathrm{IBS} = \frac{1}{\tau - t_1} \int_{t_1}^{\tau} \mathrm{BS}(t)\, dt
}
approximated by the trapezoidal rule over the supplied \code{times}.
}
\examples{
library(survival)

lung <- stats::na.omit(lung)
lung$status <- as.integer(lung$status == 2)
fit <- nbsurv(Surv(time, status) ~ age + sex, data = lung)
res <- evaluate_nbsurv(fit, newdata = lung[1:50, ], times = c(100, 200, 400), ibs = TRUE)
res
attr(res, "ibs")

cv_nbsurv(Surv(time, status) ~ age + sex, data = lung, folds = 3, times = c(100, 200))
}
