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
\value{
\code{print.nbsurv} is called for its side effect of printing a concise summary
of the fitted model (formula, number of training rows, predictors, prediction
grid size, covariance structure, and time-smoothing settings) to the console;
it returns its argument \code{x} invisibly.

\code{plot.nbsurv} is called for its side effect of drawing predicted
survival curves for the first \code{n_curves} training observations over
\code{times}; it returns, invisibly, the matrix of predicted survival
probabilities used for the plot (one row per displayed observation, one column
per element of \code{times}).
}
\description{
\code{nbsurv} objects are created by \code{\link{nbsurv}} and store the training
data, outcome, feature types, scaling information, and the Kaplan-Meier
estimates required for prediction.
}
