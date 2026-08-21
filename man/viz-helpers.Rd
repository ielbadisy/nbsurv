\name{plot.varimp_nbsurv}
\alias{plot.varimp_nbsurv}
\alias{plot.cv_nbsurv}
\alias{plot.tune_nbsurv}
\title{Visualization Helpers for nbsurv Results}
\usage{
\method{plot}{varimp_nbsurv}(x, ...)

\method{plot}{cv_nbsurv}(x, metric = c("brier", "concordance"), ...)

\method{plot}{tune_nbsurv}(x, ...)
}
\arguments{
  \item{x}{A \code{"varimp_nbsurv"}, \code{"cv_nbsurv"}, or \code{"tune_nbsurv"} object.}
  \item{metric}{For \code{plot.cv_nbsurv()}: which mean cross-validation metric to plot
    against \code{time}. One of \code{"brier"} or \code{"concordance"}.}
  \item{...}{Additional arguments passed to the underlying base graphics function
    (\code{\link[graphics]{barplot}}, \code{\link[graphics]{matplot}}, or
    \code{\link[graphics]{plot}}).}
}
\value{
The input object, invisibly. Called for its plotting side effect.
}
\description{
Base-graphics plot methods for the result objects returned by
\code{\link{varimp_nbsurv}}, \code{\link{cv_nbsurv}}, and \code{\link{tune_nbsurv}}.

\code{plot.varimp_nbsurv()} draws a horizontal bar chart of permutation importance,
sorted by increasing importance (largest bar at top).

\code{plot.cv_nbsurv()} draws the mean cross-validation metric against the evaluation
horizons.

\code{plot.tune_nbsurv()} draws the mean tuning metric for every candidate row of
\code{param_grid}, highlighting the selected best candidate in red.
}
\examples{
library(survival)

lung <- stats::na.omit(lung)
lung$status <- as.integer(lung$status == 2)
fit <- nbsurv(Surv(time, status) ~ age + sex + ph.ecog, data = lung)

vi <- varimp_nbsurv(fit, newdata = lung, times = c(100, 200, 400), n_repeats = 5, seed = 1)
plot(vi)

cvfit <- cv_nbsurv(Surv(time, status) ~ age + sex, data = lung, folds = 3, times = c(100, 200, 300))
plot(cvfit)

grid <- expand.grid(scale = TRUE, laplace = c(0.5, 1, 2), min_sd = 0.05)
tf <- tune_nbsurv(
  Surv(time, status) ~ age + sex, data = lung,
  param_grid = grid, folds = 3, times = c(100, 200)
)
plot(tf)
}
