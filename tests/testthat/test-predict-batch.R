test_that("batch predictions equal row-by-row predictions with several continuous predictors", {
  lung <- make_lung_data()
  lung <- lung[stats::complete.cases(lung[, c("age", "ph.karno", "meal.cal", "wt.loss")]), ]
  fit <- nbsurv(survival::Surv(time, status) ~ age + ph.karno + meal.cal + wt.loss, data = lung)
  nd <- lung[1:12, ]
  times <- c(100, 300, 600)
  batch <- predict(fit, newdata = nd, times = times)
  single <- do.call(rbind, lapply(seq_len(nrow(nd)), function(i)
    predict(fit, newdata = nd[i, , drop = FALSE], times = times)))
  expect_equal(unname(as.matrix(batch)), unname(as.matrix(single)), tolerance = 1e-10)
})
