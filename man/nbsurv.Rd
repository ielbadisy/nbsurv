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
  shrinkage = 0.2
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
