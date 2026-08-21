\name{calibration_plot_nbsurv}
\alias{calibration_plot_nbsurv}
\title{Calibration Plot for a nbsurv Model}
\usage{
calibration_plot_nbsurv(object, newdata, horizon, n_groups = 10, ...)
}
\arguments{
  \item{object}{A fitted \code{"nbsurv"} model.}
  \item{newdata}{A data frame containing the outcome and predictors.}
  \item{horizon}{A single numeric prediction horizon.}
  \item{n_groups}{Number of quantile-based calibration groups.}
  \item{...}{Additional graphical arguments passed to \code{plot()}.}
}
\value{
Invisibly returns a data frame with columns \code{mean_pred} (mean predicted
event probability per group), \code{observed} (observed event probability from
Kaplan-Meier), and \code{n} (group size). The calibration plot is drawn as a
side effect.
}
\description{
Assesses calibration of a \code{nbsurv} model at a single prediction horizon
by grouping observations into quantile bins of predicted event probability and
comparing the group mean prediction against the Kaplan-Meier observed event
probability within each bin.
}
\examples{
library(survival)

lung <- stats::na.omit(lung)
lung$status <- as.integer(lung$status == 2)
fit <- nbsurv(Surv(time, status) ~ age + sex, data = lung)
calibration_plot_nbsurv(fit, newdata = lung, horizon = 300)
}
