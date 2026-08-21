test_that("varimp_nbsurv returns one row per predictor, sorted by decreasing importance", {
  dat <- simulate_nbsurv_data(n = 200, censor_rate = 0.2, seed = 3)
  fit <- nbsurv(survival::Surv(time, status) ~ x1 + x2 + grp, data = dat)
  vi <- varimp_nbsurv(fit, newdata = dat, times = c(1, 2, 3), n_repeats = 3, seed = 1)

  expect_s3_class(vi, "varimp_nbsurv")
  expect_setequal(vi$feature, c("x1", "x2", "grp"))
  expect_true(is.unsorted(-vi$importance) == FALSE)
  expect_true(all(is.finite(vi$importance)))
})

test_that("varimp_nbsurv ranks a strong predictor above pure noise", {
  set.seed(11)
  n <- 500
  z <- rbinom(n, 1, 0.5)
  x1 <- rnorm(n, mean = z * 6, sd = 1)
  noise <- rnorm(n)
  true_time <- ifelse(z == 1, rweibull(n, shape = 1.4, scale = 2), rweibull(n, shape = 1.4, scale = 12))
  censor_time <- rexp(n, rate = 0.05)
  dat <- data.frame(
    time = pmin(true_time, censor_time),
    status = as.integer(true_time <= censor_time),
    x1 = x1,
    noise = noise
  )
  fit <- nbsurv(survival::Surv(time, status) ~ x1 + noise, data = dat)
  vi <- varimp_nbsurv(fit, newdata = dat, times = c(1, 2, 3), n_repeats = 15, seed = 2)

  x1_importance <- vi$importance[vi$feature == "x1"]
  noise_importance <- vi$importance[vi$feature == "noise"]
  expect_gt(x1_importance, noise_importance)
  expect_gt(x1_importance, 0)
})

test_that("varimp_nbsurv works with metric = 'concordance' (higher baseline = better)", {
  dat <- simulate_nbsurv_data(n = 200, censor_rate = 0.2, seed = 4)
  fit <- nbsurv(survival::Surv(time, status) ~ x1 + x2, data = dat)
  vi <- varimp_nbsurv(fit, newdata = dat, times = c(1, 2, 3), metric = "concordance", n_repeats = 3, seed = 1)

  expect_equal(attr(vi, "metric"), "concordance")
  expect_setequal(vi$feature, c("x1", "x2"))
})

test_that("varimp_nbsurv rejects a non-positive n_repeats", {
  dat <- simulate_nbsurv_data(n = 100, seed = 5)
  fit <- nbsurv(survival::Surv(time, status) ~ x1, data = dat)
  expect_error(varimp_nbsurv(fit, newdata = dat, times = 1, n_repeats = 0), "n_repeats")
})

test_that("print.varimp_nbsurv runs without error", {
  dat <- simulate_nbsurv_data(n = 100, seed = 6)
  fit <- nbsurv(survival::Surv(time, status) ~ x1 + x2, data = dat)
  vi <- varimp_nbsurv(fit, newdata = dat, times = c(1, 2), n_repeats = 2, seed = 1)
  expect_output(print(vi), "permutation variable importance")
})

test_that("plot.varimp_nbsurv runs without error", {
  dat <- simulate_nbsurv_data(n = 100, seed = 8)
  fit <- nbsurv(survival::Surv(time, status) ~ x1 + x2, data = dat)
  vi <- varimp_nbsurv(fit, newdata = dat, times = c(1, 2), n_repeats = 2, seed = 1)
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp)
  expect_silent(plot(vi))
  grDevices::dev.off()
  unlink(tmp)
})

test_that("plot.cv_nbsurv runs without error", {
  dat <- simulate_nbsurv_data(n = 150, seed = 9)
  cvfit <- cv_nbsurv(survival::Surv(time, status) ~ x1 + x2, data = dat, folds = 3, times = c(1, 2, 3))
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp)
  expect_silent(plot(cvfit))
  grDevices::dev.off()
  unlink(tmp)
})

test_that("plot.tune_nbsurv runs without error and highlights the best candidate", {
  dat <- simulate_nbsurv_data(n = 150, seed = 10)
  grid <- expand.grid(scale = TRUE, laplace = c(0.5, 1, 2), min_sd = 0.05)
  tf <- tune_nbsurv(survival::Surv(time, status) ~ x1 + x2, data = dat, param_grid = grid, folds = 3, times = c(1, 2))
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp)
  expect_silent(plot(tf))
  grDevices::dev.off()
  unlink(tmp)
  expect_true(tf$best_index >= 1 && tf$best_index <= nrow(tf$results))
})
