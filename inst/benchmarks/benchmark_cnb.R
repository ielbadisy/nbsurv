simulate_cnb_dgp <- function(n,
                             p_cont = 8,
                             censor_rate = 0.08,
                             horizon_quantile = 0.6,
                             seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  z <- rbinom(n, size = 1, prob = 0.5)

  x_cont <- sapply(seq_len(p_cont), function(j) {
    sd_j <- ifelse(z == 1, 2.2, 0.8)
    rnorm(n, mean = 0, sd = sd_j)
  })
  colnames(x_cont) <- paste0("x", seq_len(p_cont))

  x_cat1 <- ifelse(
    z == 1,
    sample(c("A", "B", "C"), size = n, replace = TRUE, prob = c(0.15, 0.30, 0.55)),
    sample(c("A", "B", "C"), size = n, replace = TRUE, prob = c(0.55, 0.30, 0.15))
  )
  x_cat2 <- ifelse(
    z == 1,
    sample(c("L", "M"), size = n, replace = TRUE, prob = c(0.25, 0.75)),
    sample(c("L", "M"), size = n, replace = TRUE, prob = c(0.75, 0.25))
  )

  true_time <- ifelse(
    z == 1,
    rweibull(n, shape = 1.3, scale = 3.0),
    rweibull(n, shape = 1.3, scale = 8.0)
  )
  censor_time <- rexp(n, rate = censor_rate)

  time <- pmin(true_time, censor_time)
  status <- as.integer(true_time <= censor_time)
  horizon <- unname(stats::quantile(time[status == 1], probs = horizon_quantile))

  data.frame(
    time = time,
    status = status,
    x_cont,
    x_cat1 = factor(x_cat1),
    x_cat2 = factor(x_cat2),
    z = z,
    horizon = horizon
  )
}

survival_to_event <- function(surv_matrix) {
  if (is.null(dim(surv_matrix))) {
    return(1 - surv_matrix)
  }
  1 - surv_matrix[, ncol(surv_matrix)]
}

km_event_rate <- function(time, status, horizon) {
  fit <- survival::survfit(survival::Surv(time, status) ~ 1)
  1 - summary(fit, times = horizon, extend = TRUE)$surv
}

brier_at_horizon <- function(time, status, pred_event, horizon) {
  outcome <- as.integer(time <= horizon & status == 1)
  mean((outcome - pred_event) ^ 2)
}

fit_ranger_survival <- function(formula, data) {
  ranger::ranger(
    formula = formula,
    data = data,
    num.trees = 500,
    mtry = max(2, floor(sqrt(ncol(data) - 2))),
    min.node.size = 10,
    respect.unordered.factors = "order",
    splitrule = "logrank",
    seed = 1
  )
}

predict_ranger_survival <- function(fit, newdata, horizon) {
  pred <- predict(fit, data = newdata)
  eligible <- which(pred$unique.death.times <= horizon)
  if (!length(eligible)) {
    return(rep(1, nrow(newdata)))
  }
  time_index <- max(eligible)
  pred$survival[, time_index]
}

run_one_benchmark <- function(n_train = 400, n_test = 2000, seed = 1) {
  train <- simulate_cnb_dgp(n_train, seed = seed)
  test <- simulate_cnb_dgp(n_test, seed = seed + 1000)
  horizon <- train$horizon[1]

  predictors <- setdiff(names(train), c("time", "status", "z", "horizon"))
  formula <- stats::as.formula(
    paste("survival::Surv(time, status) ~", paste(predictors, collapse = " + "))
  )

  fit_cnb <- nbsurv::nbsurv(formula, data = train)
  fit_cox <- survival::coxph(formula, data = train, x = TRUE, y = TRUE)
  fit_ranger <- fit_ranger_survival(formula, data = train)

  surv_cnb <- stats::predict(fit_cnb, newdata = test, times = horizon)
  surv_cnb_vec <- as.numeric(surv_cnb[, ncol(surv_cnb)])
  event_cnb <- 1 - surv_cnb_vec

  basehaz_cox <- survival::basehaz(fit_cox, centered = FALSE)
  haz_index <- max(which(basehaz_cox$time <= horizon))
  cumhaz_h <- if (length(haz_index)) basehaz_cox$hazard[haz_index] else 0
  lp <- stats::predict(fit_cox, newdata = test, type = "lp")
  surv_cox <- exp(-cumhaz_h * exp(lp))
  event_cox <- 1 - surv_cox

  surv_ranger <- predict_ranger_survival(fit_ranger, newdata = test, horizon = horizon)
  event_ranger <- 1 - surv_ranger

  cnb_cindex <- survival::concordance(survival::Surv(test$time, test$status) ~ surv_cnb_vec)$concordance
  cox_cindex <- survival::concordance(survival::Surv(test$time, test$status) ~ surv_cox)$concordance
  ranger_cindex <- survival::concordance(survival::Surv(test$time, test$status) ~ surv_ranger)$concordance

  data.frame(
    seed = seed,
    horizon = horizon,
    model = c("nbsurv", "CoxPH", "Ranger"),
    c_index = c(cnb_cindex, cox_cindex, ranger_cindex),
    brier = c(
      brier_at_horizon(test$time, test$status, event_cnb, horizon),
      brier_at_horizon(test$time, test$status, event_cox, horizon),
      brier_at_horizon(test$time, test$status, event_ranger, horizon)
    ),
    km_event_rate = km_event_rate(test$time, test$status, horizon)
  )
}

run_benchmark <- function(n_rep = 30, n_train = 400, n_test = 2000) {
  results <- do.call(
    rbind,
    lapply(seq_len(n_rep), function(rep_id) {
      run_one_benchmark(n_train = n_train, n_test = n_test, seed = rep_id)
    })
  )

  summary <- aggregate(
    cbind(c_index, brier) ~ model,
    data = results,
    FUN = mean
  )
  summary$c_index_sd <- aggregate(c_index ~ model, data = results, FUN = stats::sd)$c_index
  summary$brier_sd <- aggregate(brier ~ model, data = results, FUN = stats::sd)$brier

  list(raw = results, summary = summary)
}

if (sys.nframe() == 0) {
  stopifnot(requireNamespace("nbsurv", quietly = TRUE))
  stopifnot(requireNamespace("ranger", quietly = TRUE))
  stopifnot(requireNamespace("survival", quietly = TRUE))

  out <- run_benchmark()
  dir.create("results", recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(out$raw, file = "results/benchmark_raw.csv", row.names = FALSE)
  utils::write.csv(out$summary, file = "results/benchmark_summary.csv", row.names = FALSE)
  print(out$summary)
}
