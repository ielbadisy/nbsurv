test_that("nbsurv returns valid monotone survival predictions", {
  skip_if_not_installed("survival")
  lung <- make_lung_data()

  fit <- nbsurv(survival::Surv(time, status) ~ age + sex + ph.ecog, data = lung)
  preds <- predict(fit, newdata = lung[1:10, ], times = c(100, 200, 400, 800))

  expect_s3_class(fit, "nbsurv")
  expect_equal(dim(preds), c(10, 4))
  expect_true(all(is.finite(preds)))
  expect_true(all(preds >= 0 & preds <= 1))
  expect_true(all(t(apply(preds, 1, diff)) <= 1e-10))
})

test_that("event predictions complement survival predictions", {
  skip_if_not_installed("survival")
  lung <- make_lung_data()

  fit <- nbsurv(survival::Surv(time, status) ~ age + sex, data = lung)
  surv_preds <- predict(fit, newdata = lung[1:5, ], times = c(200, 500))
  event_preds <- predict(fit, newdata = lung[1:5, ], times = c(200, 500), type = "event")

  expect_equal(event_preds, 1 - surv_preds, tolerance = 1e-10)
  expect_equal(
    predictSurvProb(fit, newdata = lung[1:5, ], times = c(200, 500)),
    surv_preds,
    tolerance = 1e-10
  )
})

test_that("categorical predictors reject unseen levels at prediction time", {
  skip_if_not_installed("survival")
  lung <- make_lung_data()

  fit <- nbsurv(survival::Surv(time, status) ~ age + sex, data = lung)
  newdata <- lung[1:2, ]
  newdata$sex <- factor(c("3", "3"))

  expect_error(
    predict(fit, newdata = newdata, times = 200),
    regexp = "factor"
  )
})

test_that("evaluate_nbsurv returns finite metrics", {
  skip_if_not_installed("survival")
  lung <- make_lung_data()

  train <- lung[1:140, ]
  test <- lung[141:nrow(lung), ]
  fit <- nbsurv(survival::Surv(time, status) ~ age + sex + ph.ecog, data = train)

  metrics <- evaluate_nbsurv(fit, newdata = test, times = c(150, 300, 450))

  expect_equal(names(metrics), c("time", "brier", "concordance"))
  expect_equal(nrow(metrics), 3)
  expect_true(all(is.finite(metrics$brier)))
  expect_true(all(metrics$brier >= 0))
  expect_true(all(is.finite(metrics$concordance)))
  expect_true(all(metrics$concordance >= 0 & metrics$concordance <= 1))
})

test_that("average predicted event risk is in line with observed horizon risk", {
  dat <- simulate_nbsurv_data(n = 300, censor_rate = 0.15, seed = 12)
  fit <- nbsurv(survival::Surv(time, status) ~ x1 + x2 + grp, data = dat)

  horizon <- as.numeric(stats::quantile(dat$time[dat$status == 1], probs = 0.5, type = 2))
  pred_event <- as.numeric(predict(fit, newdata = dat, times = horizon, type = "event")[, 1])

  km <- survival::survfit(survival::Surv(time, status) ~ 1, data = dat)
  observed_event <- 1 - summary(km, times = horizon, extend = TRUE)$surv

  expect_lt(abs(mean(pred_event) - observed_event), 0.15)
})

test_that("cv_nbsurv returns fold-level and summary metrics", {
  skip_if_not_installed("survival")
  lung <- make_lung_data()

  cv_fit <- cv_nbsurv(
    survival::Surv(time, status) ~ age + sex + ph.ecog,
    data = lung,
    folds = 3,
    times = c(150, 300),
    seed = 42
  )

  expect_s3_class(cv_fit, "cv_nbsurv")
  expect_equal(sort(unique(cv_fit$fold_metrics$fold)), 1:3)
  expect_equal(sort(unique(cv_fit$summary$time)), c(150, 300))
  expect_true(all(is.finite(cv_fit$summary$brier)))
  expect_true(all(is.finite(cv_fit$summary$concordance)))
})

test_that("tune_nbsurv ranks parameter settings and returns the best row", {
  skip_if_not_installed("survival")
  lung <- make_lung_data()

  param_grid <- data.frame(
    scale = c(TRUE, FALSE),
    laplace = c(1, 2),
    min_sd = c(0.05, 0.10)
  )
  param_grid$time_grid <- I(list(NULL, NULL))

  tuned <- tune_nbsurv(
    survival::Surv(time, status) ~ age + sex + ph.ecog,
    data = lung,
    param_grid = param_grid,
    folds = 3,
    times = c(150, 300),
    seed = 7
  )

  expect_s3_class(tuned, "tune_nbsurv")
  expect_equal(nrow(tuned$results), 2)
  expect_equal(nrow(tuned$best_params), 1)
  expect_true(is.finite(tuned$results$mean_metric[1]))
})

test_that("fit validation catches unsupported responses and invalid times", {
  skip_if_not_installed("survival")
  lung <- make_lung_data()

  expect_error(
    nbsurv(status ~ age + sex, data = lung),
    regexp = "Surv"
  )

  fit <- nbsurv(survival::Surv(time, status) ~ age + sex, data = lung)
  expect_error(predict(fit, newdata = lung[1:3, ], times = c(100, NA)), regexp = "non-missing")
})

test_that("heavy censoring still yields finite predictions and evaluation metrics", {
  dat <- simulate_nbsurv_data(n = 260, censor_rate = 0.8, seed = 99)
  fit <- nbsurv(survival::Surv(time, status) ~ x1 + x2 + grp, data = dat, eps = 1e-04)

  eval_times <- as.numeric(stats::quantile(dat$time[dat$status == 1], probs = c(0.25, 0.5), type = 2))
  preds <- predict(fit, newdata = dat[1:20, ], times = eval_times)
  metrics <- evaluate_nbsurv(fit, newdata = dat, times = eval_times)

  expect_true(all(is.finite(preds)))
  expect_true(all(is.finite(metrics$brier)))
  expect_true(all(is.finite(metrics$concordance)))
})

test_that("evaluate_nbsurv returns finite IBS when ibs = TRUE", {
  skip_if_not_installed("survival")
  lung <- make_lung_data()

  train <- lung[1:140, ]
  test  <- lung[141:nrow(lung), ]
  fit <- nbsurv(survival::Surv(time, status) ~ age + sex + ph.ecog, data = train)

  res <- evaluate_nbsurv(fit, newdata = test, times = c(150, 300, 450), ibs = TRUE)
  ibs_val <- attr(res, "ibs")

  expect_true(is.finite(ibs_val))
  expect_true(ibs_val >= 0 && ibs_val <= 1)
})

test_that("calibration_plot_nbsurv returns a data frame and plots without error", {
  skip_if_not_installed("survival")
  lung <- make_lung_data()
  fit <- nbsurv(survival::Surv(time, status) ~ age + sex, data = lung)

  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  calib <- calibration_plot_nbsurv(fit, newdata = lung, horizon = 300)

  expect_s3_class(calib, "data.frame")
  expect_true(all(c("mean_pred", "observed", "n") %in% names(calib)))
  expect_true(all(is.finite(calib$mean_pred)))
  expect_true(all(is.finite(calib$observed)))
})

test_that("print and plot methods execute without error", {
  skip_if_not_installed("survival")
  lung <- make_lung_data()
  fit <- nbsurv(survival::Surv(time, status) ~ age + sex, data = lung)

  printed <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(printed, "nbsurv conditional naive Bayes model")

  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  plotted <- plot(fit, times = c(100, 200, 400), n_curves = 3)

  expect_equal(dim(plotted), c(3, 3))
})
