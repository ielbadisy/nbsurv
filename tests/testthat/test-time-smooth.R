test_that("time_smooth = FALSE (default) is unaffected: no grid_stats stored", {
  dat <- simulate_nbsurv_data(n = 100, seed = 30)
  fit <- nbsurv(survival::Surv(time, status) ~ x1 + x2, data = dat)
  expect_null(fit$grid_stats)
  expect_false(isTRUE(fit$time_smooth))
})

test_that("time_smooth = TRUE produces finite, valid predictions (regression guard for the NA bug)", {
  # nw_smooth_matrix() previously propagated NA from grid points with no
  # events/survivors into every prediction, because 0 * NA is still NA in R
  # and colSums() lacked na.rm - every smoothed prediction was NA.
  dat <- simulate_nbsurv_data(n = 150, seed = 31)
  fit <- nbsurv(survival::Surv(time, status) ~ x1 + x2, data = dat, time_smooth = TRUE)
  expect_false(is.null(fit$grid_stats))
  expect_true(is.finite(fit$bandwidth) && fit$bandwidth > 0)

  preds <- predict(fit, newdata = dat[1:10, ], times = c(1, 2, 3))
  expect_true(all(is.finite(preds)))
  expect_true(all(preds >= 0 & preds <= 1))

  m <- evaluate_nbsurv(fit, newdata = dat, times = c(1, 2, 3))
  expect_true(all(is.finite(m$brier)))
  expect_true(all(is.finite(m$concordance)))
})

test_that("an explicit bandwidth is honored and stored on the fit", {
  dat <- simulate_nbsurv_data(n = 120, seed = 32)
  fit <- nbsurv(survival::Surv(time, status) ~ x1 + x2, data = dat, time_smooth = TRUE, bandwidth = 5)
  expect_equal(fit$bandwidth, 5)
})

test_that("print.nbsurv reports time smoothing when active", {
  dat <- simulate_nbsurv_data(n = 100, seed = 33)
  fit <- nbsurv(survival::Surv(time, status) ~ x1 + x2, data = dat, time_smooth = TRUE, bandwidth = 10)
  expect_output(print(fit), "Time smoothing: kernel bandwidth = 10")
})

test_that("time_smooth improves held-out concordance at a sparse extreme horizon vs. diagonal", {
  # At an extreme horizon, one of the two horizon-defined classes (survivors
  # or IPCW-weighted events) becomes small, inflating estimation variance.
  # Smoothing across neighboring horizons should reduce - not necessarily
  # eliminate - that variance. Averaged over repeated splits for stability.
  set.seed(44)
  n <- 400
  z <- rbinom(n, 1, 0.5)
  x1 <- rnorm(n, mean = z * 3, sd = 1.2)
  x2 <- rnorm(n, mean = 0.6 * x1)
  true_time <- ifelse(z == 1, rweibull(n, shape = 1.3, scale = 4), rweibull(n, shape = 1.3, scale = 25))
  censor_time <- rexp(n, rate = 1 / 40)
  dat <- data.frame(
    time = pmin(true_time, censor_time),
    status = as.integer(true_time <= censor_time),
    x1 = x1, x2 = x2
  )

  ipcw_brier <- nbsurv:::ipcw_brier_score
  conc_at <- nbsurv:::concordance_at_horizon
  horizon <- 20  # near the tail, where one class is sparse

  diag_conc <- numeric(8)
  smooth_conc <- numeric(8)
  for (rep in seq_len(8)) {
    set.seed(100 + rep)
    idx <- sample.int(n, floor(0.7 * n))
    train <- dat[idx, ]
    test <- dat[-idx, ]

    fit_diag <- nbsurv(survival::Surv(time, status) ~ x1 + x2, data = train)
    fit_smooth <- nbsurv(survival::Surv(time, status) ~ x1 + x2, data = train, time_smooth = TRUE, bandwidth = 15)

    diag_event <- predict(fit_diag, newdata = test, times = horizon, type = "event")
    smooth_event <- predict(fit_smooth, newdata = test, times = horizon, type = "event")

    diag_conc[rep] <- conc_at(test$time, test$status, diag_event[, 1])
    smooth_conc[rep] <- conc_at(test$time, test$status, smooth_event[, 1])
  }

  expect_gt(mean(smooth_conc, na.rm = TRUE), mean(diag_conc, na.rm = TRUE) - 0.02)
})

test_that("tune_nbsurv can search over time_smooth/bandwidth via param_grid", {
  dat <- simulate_nbsurv_data(n = 150, seed = 34)
  grid <- data.frame(
    scale = TRUE, laplace = 1, min_sd = 0.05,
    time_smooth = c(FALSE, TRUE, TRUE),
    bandwidth = c(NA, 5, 15)
  )
  tuned <- tune_nbsurv(
    survival::Surv(time, status) ~ x1 + x2, data = dat,
    param_grid = grid, folds = 3, times = c(1, 2, 3), seed = 1
  )
  expect_equal(nrow(tuned$results), 3)
  expect_true(all(is.finite(tuned$results$mean_metric)))
})
