#' Deprecated functions
#'
#' These underscore names are kept for one release and will be removed in the next.
#' Use the glued lowercase names instead.
#'
#' @param ... Passed to the new function.
#' @return The value of the new function.
#' @name nbsurv-deprecated
#' @keywords internal
NULL

#' @rdname nbsurv-deprecated
#' @export
calibration_plot_nbsurv <- function(...) {
  .Deprecated("calibrationplotnbsurv", package = "nbsurv")
  calibrationplotnbsurv(...)
}

#' @rdname nbsurv-deprecated
#' @export
cv_nbsurv <- function(...) {
  .Deprecated("cvnbsurv", package = "nbsurv")
  cvnbsurv(...)
}

#' @rdname nbsurv-deprecated
#' @export
evaluate_nbsurv <- function(...) {
  .Deprecated("evaluatenbsurv", package = "nbsurv")
  evaluatenbsurv(...)
}

#' @rdname nbsurv-deprecated
#' @export
tune_nbsurv <- function(...) {
  .Deprecated("tunenbsurv", package = "nbsurv")
  tunenbsurv(...)
}

#' @rdname nbsurv-deprecated
#' @export
varimp_nbsurv <- function(...) {
  .Deprecated("varimpnbsurv", package = "nbsurv")
  varimpnbsurv(...)
}
