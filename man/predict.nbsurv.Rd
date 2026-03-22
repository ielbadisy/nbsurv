\name{predict.nbsurv}
\alias{predict.nbsurv}
\alias{predictSurvProb}
\alias{predictSurvProb.nbsurv}
\title{Predict Survival Probabilities from a nbsurv Model}
\usage{
\method{predict}{nbsurv}(object, newdata, times, type = c("survival", "event"), ...)

predictSurvProb(object, ...)

\method{predictSurvProb}{nbsurv}(object, newdata, times, ...)
}
\arguments{
  \item{object}{A fitted \code{"nbsurv"} model.}
  \item{newdata}{A data frame of predictor values.}
  \item{times}{Numeric vector of prediction horizons.}
  \item{type}{Either \code{"survival"} for survival probabilities or \code{"event"} for event probabilities.}
  \item{...}{Unused additional arguments.}
}
\value{
A matrix with one row per observation in \code{newdata} and one column per
prediction horizon.
}
\description{
Computes horizon-specific survival or event probabilities from a fitted
conditional naive Bayes survival model. The returned survival curves are
forced to be non-increasing across increasing time horizons.
}
\examples{
library(survival)

lung <- stats::na.omit(lung)
lung$status <- as.integer(lung$status == 2)
fit <- nbsurv(Surv(time, status) ~ age + sex, data = lung)
predict(fit, newdata = lung[1:2, ], times = c(100, 200))
predict(fit, newdata = lung[1:2, ], times = c(100, 200), type = "event")
}
