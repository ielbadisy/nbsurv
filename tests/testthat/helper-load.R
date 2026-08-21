if (!exists("nbsurv", mode = "function")) {
  pkgload::load_all(export_all = FALSE, helpers = FALSE, quiet = TRUE)
}

make_lung_data <- function() {
  lung <- survival::lung
  lung$status <- as.integer(lung$status == 2)
  lung$sex <- factor(lung$sex)
  lung
}

simulate_nbsurv_data <- function(n = 250, censor_rate = 0.2, seed = 1) {
  set.seed(seed)

  z <- rbinom(n, size = 1, prob = 0.45)
  x1 <- rnorm(n, mean = z * 0.5, sd = ifelse(z == 1, 1.8, 0.9))
  x2 <- rnorm(n, mean = 0, sd = ifelse(z == 1, 1.5, 0.7))
  grp <- factor(ifelse(z == 1, sample(c("A", "B"), n, TRUE, c(0.3, 0.7)), sample(c("A", "B"), n, TRUE, c(0.7, 0.3))))

  true_time <- ifelse(
    z == 1,
    rweibull(n, shape = 1.4, scale = 2.8),
    rweibull(n, shape = 1.4, scale = 6.5)
  )
  censor_time <- rexp(n, rate = censor_rate)

  data.frame(
    time = pmin(true_time, censor_time),
    status = as.integer(true_time <= censor_time),
    x1 = x1,
    x2 = x2,
    grp = grp
  )
}
