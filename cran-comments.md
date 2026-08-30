# CRAN submission comments — nbsurv 0.5.1

## Resubmission

This is a resubmission addressing the review of 0.5.0.

* Added references describing the methods to the `Description` field of
  `DESCRIPTION`, in the `authors (year) <doi:...>` form:
  Wolfson et al. (2015) <doi:10.1002/sim.6526> for the censored naive
  Bayes model implemented by `nbsurv()`, and Gerds and Schumacher (2006)
  <doi:10.1002/bimj.200610301> for the inverse-probability-of-censoring
  weighted Brier score used by `evaluate_nbsurv()` / `cv_nbsurv()`.
* Added a `\value` section to `man/nbsurv-class.Rd`, documenting the
  return values (structure/class and meaning) of the exported `print()`
  and `plot()` methods for `"nbsurv"` objects. All other exported
  functions already had `\value`.

## Test environments

* Local: Ubuntu 24.04, R 4.4.x
* win-builder: R-devel
* rhub: windows-x86_64-devel, ubuntu-gcc-devel

## R CMD check results

0 ERRORs | 0 WARNINGs | 1 NOTE

The NOTE is "New submission". The `Description` field may also flag the
author surnames (Wolfson, Gerds, Schumacher) and "Kaplan-Meier" as
possibly misspelled words; these are proper names and are spelled
correctly.

## Reverse dependencies

None (new package).
