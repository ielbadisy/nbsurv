test_that("numeric 0/1 predictors are categorical and match a factor coding", {
  set.seed(1); n <- 300; trt <- rbinom(n, 1, .5); x <- rnorm(n)
  d <- data.frame(time = rexp(n, exp(-.5 * trt + .5 * x)), status = rbinom(n, 1, .8), trt = trt, x = x)
  f1 <- nbsurv(Surv(time, status) ~ trt + x, d)
  f2 <- nbsurv(Surv(time, status) ~ trt + x, transform(d, trt = factor(trt)))
  expect_identical(unname(f1$feature_info[["trt"]]), "categorical")
  nd <- data.frame(trt = c(0, 1, 1), x = c(0, 0, 1))
  expect_equal(unname(as.matrix(predict(f1, nd, times = c(.5, 1)))),
               unname(as.matrix(predict(f2, transform(nd, trt = factor(trt, 0:1)), times = c(.5, 1)))))
})
