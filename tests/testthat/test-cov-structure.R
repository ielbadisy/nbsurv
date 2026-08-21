rmvnorm2 <- function(n, rho) {
  z <- matrix(rnorm(2 * n), ncol = 2)
  L <- chol(matrix(c(1, rho, rho, 1), 2, 2))
  z %*% L
}

test_that("cov_structure = 'full' reduces to 'diagonal' with a single continuous predictor", {
  dat <- simulate_nbsurv_data(n = 200, seed = 21)
  dat$x2 <- NULL
  dat$grp <- NULL
  fit_diag <- nbsurv(survival::Surv(time, status) ~ x1, data = dat)
  fit_full <- nbsurv(survival::Surv(time, status) ~ x1, data = dat, cov_structure = "full")

  p_diag <- predict(fit_diag, newdata = dat, times = c(1, 2, 3))
  p_full <- predict(fit_full, newdata = dat, times = c(1, 2, 3))

  expect_equal(p_diag, p_full)
})

test_that("cov_structure errors on an invalid value", {
  dat <- simulate_nbsurv_data(n = 100, seed = 22)
  expect_error(
    nbsurv(survival::Surv(time, status) ~ x1 + x2, data = dat, cov_structure = "bogus"),
    "arg"
  )
})

test_that("cov_structure = 'full' improves Brier score under correlated continuous predictors", {
  set.seed(42)
  n <- 600
  Z <- rmvnorm2(n, rho = 0.85)
  x1 <- Z[, 1]
  x2 <- Z[, 2]
  lp <- 0.9 * x1 + 0.9 * x2
  event_time <- rexp(n, rate = exp(lp) / 50)
  censor_time <- rexp(n, rate = 1 / 150)
  dat <- data.frame(
    time = pmin(event_time, censor_time),
    status = as.integer(event_time <= censor_time),
    x1 = x1,
    x2 = x2
  )

  train_idx <- 1:400
  train <- dat[train_idx, ]
  test <- dat[-train_idx, ]

  fit_diag <- nbsurv(survival::Surv(time, status) ~ x1 + x2, data = train)
  fit_full <- nbsurv(survival::Surv(time, status) ~ x1 + x2, data = train, cov_structure = "full")

  times <- c(20, 40, 80)
  m_diag <- evaluate_nbsurv(fit_diag, newdata = test, times = times)
  m_full <- evaluate_nbsurv(fit_full, newdata = test, times = times)

  expect_lt(mean(m_full$brier), mean(m_diag$brier))
})

test_that("print.nbsurv reports the covariance structure", {
  dat <- simulate_nbsurv_data(n = 100, seed = 23)
  fit <- nbsurv(survival::Surv(time, status) ~ x1 + x2, data = dat, cov_structure = "full")
  expect_output(print(fit), "Covariance structure: full")
})

test_that("cv_nbsurv passes cov_structure/shrinkage through to each fold's fit", {
  dat <- simulate_nbsurv_data(n = 200, seed = 24)
  cv_fit <- cv_nbsurv(
    survival::Surv(time, status) ~ x1 + x2, data = dat,
    folds = 3, times = c(1, 2, 3), seed = 1, cov_structure = "full"
  )
  expect_s3_class(cv_fit, "cv_nbsurv")
  expect_true(all(is.finite(cv_fit$summary$brier)))
})

test_that("full covariance is symmetric positive-definite after shrinkage (no mahalanobis errors)", {
  dat <- simulate_nbsurv_data(n = 150, seed = 25)
  fit <- nbsurv(survival::Surv(time, status) ~ x1 + x2, data = dat, cov_structure = "full", shrinkage = 0.05)
  preds <- predict(fit, newdata = dat, times = c(1, 2, 3))
  expect_true(all(is.finite(preds)))
  expect_true(all(preds >= 0 & preds <= 1))
})
