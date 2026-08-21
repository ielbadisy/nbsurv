\name{nbsurv}
\alias{nbsurv}
\title{Fit a Conditional Naive Bayes Survival Model}
\usage{
nbsurv(
  formula,
  data,
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
}
\arguments{
  \item{formula}{A model formula with a \code{survival::Surv()} response.}
  \item{data}{A data frame containing the outcome and predictors.}
  \item{scale}{Logical; if \code{TRUE}, continuous predictors are standardized.}
  \item{laplace}{Laplace smoothing constant for categorical likelihoods.}
  \item{min_sd}{Lower bound used for Gaussian standard deviations.}
  \item{time_grid}{Optional numeric vector of prediction horizons stored with the fit.}
  \item{eps}{Small numerical constant used to clip probabilities away from 0 and 1.}
  \item{cov_structure}{\code{"diagonal"} (default): continuous predictors are modeled as
    conditionally independent (the classic naive Bayes assumption), each with its own
    class-conditional Gaussian. \code{"full"}: continuous predictors are modeled jointly as a
    single multivariate Gaussian per class (with covariance shrunk toward diagonal by
    \code{shrinkage}), relaxing the independence assumption between them. Categorical
    predictors are always treated as conditionally independent regardless of this setting.
    Requires 2+ continuous predictors to have any effect; with only one, \code{"full"} is
    numerically identical to \code{"diagonal"}. Empirically improves the IPCW Brier score
    (calibration) when continuous predictors are correlated, with a smaller and less
    consistent effect on concordance (ranking) - see \code{vignette("nbsurv-workflow")}.}
  \item{shrinkage}{Shrinkage weight toward the diagonal covariance, in `[0, 1]`, used only
    when \code{cov_structure = "full"}. `0` uses the raw sample covariance (more variance,
    especially at small `n`); `1` is equivalent to `"diagonal"`. Default `0.2`.}
  \item{time_smooth}{Logical, default \code{FALSE}. If \code{TRUE}, class-conditional
    statistics (means, covariances, categorical probabilities) are estimated once at every
    point of \code{time_grid} at fit time, then combined across horizons via Nadaraya-Watson
    kernel regression (Gaussian kernel, width \code{bandwidth}) when \code{predict()} is
    called at a given horizon, instead of being re-estimated from scratch, independently, at
    each requested horizon. This borrows statistical strength from neighboring horizons,
    reducing estimation variance where one of the two horizon-defined classes (subjects who
    already failed vs. subjects known to survive past it) is sparse - typically near the
    earliest or latest horizons requested - while still allowing genuinely time-varying
    effects (unlike a proportional-hazards model, which shares a single coefficient across
    all of time). A convex (weight-normalized) combination of the per-grid-point statistics
    is used, so smoothed covariance matrices stay positive semi-definite and smoothed
    categorical probabilities stay valid automatically, with no extra correction needed.
    Increases fit time (every grid point's statistics are computed once, up front) but not
    predict time. Empirically (see `NEWS.md`) trades ranking robustness (concordance) for
    calibration at extreme horizons as `bandwidth` widens, and vice versa as it narrows -
    tune `bandwidth` via \code{\link{tune_nbsurv}} rather than relying on the default for a
    specific dataset.}
  \item{bandwidth}{Kernel bandwidth (in the same units as \code{time}), used only when
    \code{time_smooth = TRUE}. Defaults to one quarter of the fitted \code{time_grid}'s range
    when \code{NULL}.}
}
\value{
An object of class \code{"nbsurv"}.
}
\description{
Fits a conditional naive Bayes model for right-censored survival data.
The model combines the marginal Kaplan-Meier survival estimate with
horizon-specific likelihood terms estimated among subjects who fail
before the horizon and subjects known to survive beyond it.
}
\examples{
library(survival)

lung <- stats::na.omit(lung)
lung$status <- as.integer(lung$status == 2)
fit <- nbsurv(
  Surv(time, status) ~ age + sex,
  data = lung
)

predict(fit, newdata = lung[1:3, ], times = c(100, 200, 300))
}
