\name{nbsurv-class}
\alias{print.nbsurv}
\alias{plot.nbsurv}
\title{Conditional Naive Bayes Survival Model Object}
\usage{
\method{print}{nbsurv}(x, ...)

\method{plot}{nbsurv}(x, times = NULL, n_curves = 5, ...)
}
\arguments{
  \item{x}{A fitted \code{"nbsurv"} model object.}
  \item{times}{Optional vector of horizons used for plotting.}
  \item{n_curves}{Number of training-set curves to display in the plot method.}
  \item{...}{Additional arguments passed to the underlying plotting functions.}
}
\description{
\code{nbsurv} objects are created by \code{\link{nbsurv}} and store the training
data, outcome, feature types, scaling information, and the Kaplan-Meier
estimates required for prediction.
}
