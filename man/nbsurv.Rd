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
  eps = 1e-06
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
