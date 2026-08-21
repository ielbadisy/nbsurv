nbsurv <- function(formula,
                   data,
                   scale = TRUE,
                   laplace = 1,
                   min_sd = 0.05,
                   time_grid = NULL,
                   eps = 1e-06) {
  mf <- stats::model.frame(formula, data = data, na.action = stats::na.omit)
  y <- stats::model.response(mf)

  validate_surv_response(y)

  x <- mf[, -1, drop = FALSE]
  if (!ncol(x)) {
    stop("The model requires at least one predictor.")
  }

  feature_info <- vapply(x, classify_feature, character(1))
  x <- coerce_features(x, feature_info)

  scaling <- NULL
  cont_vars <- names(feature_info[feature_info == "continuous"])
  if (scale && length(cont_vars) > 0L) {
    means <- colMeans(x[, cont_vars, drop = FALSE])
    sds <- apply(x[, cont_vars, drop = FALSE], 2, stats::sd)
    sds <- floor_sd(sds, min_sd = min_sd)
    x[, cont_vars] <- scale_continuous(x[, cont_vars, drop = FALSE], means, sds)
    scaling <- list(means = means, sds = sds)
  }

  times <- y[, "time"]
  status <- y[, "status"]

  km_event <- survival::survfit(y ~ 1)
  km_censor <- survival::survfit(survival::Surv(times, 1 - status) ~ 1)
  time_grid <- normalise_time_grid(time_grid, times = times, status = status)

  object <- list(
    call = match.call(),
    formula = formula,
    terms = stats::terms(formula, data = data),
    xlevels = lapply(x[, feature_info == "categorical", drop = FALSE], levels),
    feature_info = feature_info,
    training_data = x,
    response = y,
    times = times,
    status = status,
    event_survival = stats::stepfun(km_event$time, c(1, km_event$surv)),
    censor_survival = stats::stepfun(km_censor$time, c(1, km_censor$surv)),
    scaling = scaling,
    laplace = laplace,
    min_sd = min_sd,
    eps = eps,
    time_grid = time_grid
  )

  class(object) <- "nbsurv"
  object
}

cv_nbsurv <- function(formula,
                      data,
                      folds = 5,
                      times = NULL,
                      seed = NULL,
                      scale = TRUE,
                      laplace = 1,
                      min_sd = 0.05,
                      time_grid = NULL,
                      eps = 1e-06) {
  mf <- stats::model.frame(formula, data = data, na.action = stats::na.omit)
  data_complete <- data[rownames(mf), , drop = FALSE]
  y <- stats::model.response(mf)

  validate_surv_response(y)
  validate_folds(folds, nrow(mf))

  eval_times <- normalise_eval_times(times, y)
  fold_id <- make_fold_ids(nrow(mf), folds = folds, seed = seed)

  fold_metrics <- do.call(
    rbind,
    lapply(seq_len(folds), function(fold) {
      train <- data_complete[fold_id != fold, , drop = FALSE]
      test <- data_complete[fold_id == fold, , drop = FALSE]

      fit <- nbsurv(
        formula = formula,
        data = train,
        scale = scale,
        laplace = laplace,
        min_sd = min_sd,
        time_grid = time_grid,
        eps = eps
      )

      metrics <- evaluate_nbsurv(fit, newdata = test, times = eval_times)
      metrics$fold <- fold
      metrics
    })
  )

  summary <- stats::aggregate(
    cbind(brier, concordance) ~ time,
    data = fold_metrics,
    FUN = mean
  )

  out <- list(
    call = match.call(),
    formula = formula,
    folds = folds,
    fold_id = fold_id,
    times = eval_times,
    fold_metrics = fold_metrics,
    summary = summary
  )
  class(out) <- "cv_nbsurv"
  out
}

evaluate_nbsurv <- function(object,
                            newdata,
                            times,
                            metrics = c("brier", "concordance"),
                            ibs = FALSE) {
  metrics <- unique(match.arg(metrics, choices = c("brier", "concordance"), several.ok = TRUE))
  times <- validate_prediction_times(times)

  mf <- stats::model.frame(object$formula, data = newdata, na.action = stats::na.fail)
  y <- stats::model.response(mf)
  validate_surv_response(y)

  preds <- predict(object, newdata = mf, times = times, type = "event")
  obs_time <- y[, "time"]
  obs_status <- y[, "status"]

  results <- data.frame(time = times)
  if ("brier" %in% metrics) {
    results$brier <- vapply(
      seq_along(times),
      function(j) ipcw_brier_score(
        time = obs_time,
        status = obs_status,
        pred_event = preds[, j],
        horizon = times[j],
        eps = object$eps
      ),
      numeric(1)
    )
  }
  if ("concordance" %in% metrics) {
    results$concordance <- vapply(
      seq_along(times),
      function(j) concordance_at_horizon(
        time = obs_time,
        status = obs_status,
        risk = preds[, j]
      ),
      numeric(1)
    )
  }

  if (ibs && "brier" %in% metrics && length(times) > 1L) {
    attr(results, "ibs") <- integrated_brier_score(results$brier, times)
  }

  results
}

calibration_plot_nbsurv <- function(object, newdata, horizon, n_groups = 10, ...) {
  times <- validate_prediction_times(horizon)
  if (length(times) != 1L) {
    stop("`horizon` must be a single prediction time.")
  }

  mf <- stats::model.frame(object$formula, data = newdata, na.action = stats::na.fail)
  y <- stats::model.response(mf)
  validate_surv_response(y)

  pred_event <- as.numeric(predict(object, newdata = mf, times = times, type = "event")[, 1L])
  obs_time <- y[, "time"]
  obs_status <- y[, "status"]

  breaks <- unique(stats::quantile(pred_event, probs = seq(0, 1, length.out = n_groups + 1L), names = FALSE))
  if (length(breaks) < 3L) {
    stop("Too few unique predicted values to form calibration groups.")
  }
  grp <- base::findInterval(pred_event, breaks, rightmost.closed = TRUE)

  calib <- do.call(rbind, lapply(sort(unique(grp)), function(g) {
    idx <- grp == g
    km <- survival::survfit(survival::Surv(obs_time[idx], obs_status[idx]) ~ 1)
    s_hat <- summary(km, times = times[1L], extend = TRUE)$surv
    if (!length(s_hat)) s_hat <- 1
    data.frame(
      mean_pred = mean(pred_event[idx]),
      observed  = 1 - s_hat,
      n         = sum(idx)
    )
  }))

  graphics::plot(
    calib$mean_pred, calib$observed,
    xlab = "Mean predicted event probability",
    ylab = "Observed event probability (1 - KM)",
    pch = 16, xlim = c(0, 1), ylim = c(0, 1), ...
  )
  graphics::abline(0, 1, lty = 2)
  invisible(calib)
}

tune_nbsurv <- function(formula,
                        data,
                        param_grid,
                        folds = 5,
                        times = NULL,
                        metric = c("brier", "concordance"),
                        maximize = NULL,
                        seed = NULL,
                        eps = 1e-06) {
  metric <- match.arg(metric)
  param_grid <- validate_param_grid(param_grid)

  if (is.null(maximize)) {
    maximize <- identical(metric, "concordance")
  }

  results <- do.call(
    rbind,
    lapply(seq_len(nrow(param_grid)), function(i) {
      params <- as.list(param_grid[i, , drop = FALSE])
      cv_fit <- cv_nbsurv(
        formula = formula,
        data = data,
        folds = folds,
        times = times,
        seed = seed,
        scale = params$scale,
        laplace = params$laplace,
        min_sd = params$min_sd,
        time_grid = params$time_grid[[1]],
        eps = eps
      )

      data.frame(
        param_grid[i, , drop = FALSE],
        mean_metric = mean(cv_fit$summary[[metric]], na.rm = TRUE)
      )
    })
  )

  best_index <- if (maximize) {
    which.max(results$mean_metric)
  } else {
    which.min(results$mean_metric)
  }

  out <- list(
    call = match.call(),
    metric = metric,
    maximize = maximize,
    results = results,
    best_index = best_index,
    best_params = results[best_index, , drop = FALSE]
  )
  class(out) <- "tune_nbsurv"
  out
}

predict.nbsurv <- function(object,
                           newdata,
                           times,
                           type = c("survival", "event"),
                           ...) {
  type <- match.arg(type)
  times <- validate_prediction_times(times)
  new_x <- prepare_newdata(object, newdata)

  surv_probs <- matrix(NA_real_, nrow = nrow(new_x), ncol = length(times))
  colnames(surv_probs) <- paste0("t_", format(times, trim = TRUE, scientific = FALSE))
  rownames(surv_probs) <- rownames(new_x)

  train_x <- object$training_data
  cont_vars <- names(object$feature_info[object$feature_info == "continuous"])
  cat_vars <- names(object$feature_info[object$feature_info == "categorical"])

  for (j in seq_along(times)) {
    horizon <- times[j]
    prior_surv <- clip_prob(object$event_survival(horizon), object$eps)
    event_weights <- event_case_weights(object, horizon)
    survivor_index <- object$times >= horizon

    log_surv <- rep(log(prior_surv), nrow(new_x))
    log_event <- rep(log1p(-prior_surv), nrow(new_x))

    if (length(cont_vars) > 0L) {
      cont_stats <- continuous_statistics(
        train_x = train_x[, cont_vars, drop = FALSE],
        event_weights = event_weights,
        survivor_index = survivor_index,
        min_sd = object$min_sd
      )

      if (!is.null(cont_stats)) {
        new_cont <- as.matrix(new_x[, cont_vars, drop = FALSE])
        log_surv <- log_surv + rowSums(stats::dnorm(
          x = new_cont,
          mean = cont_stats$survivor_mean,
          sd = cont_stats$survivor_sd,
          log = TRUE
        ))
        log_event <- log_event + rowSums(stats::dnorm(
          x = new_cont,
          mean = cont_stats$event_mean,
          sd = cont_stats$event_sd,
          log = TRUE
        ))
      }
    }

    if (length(cat_vars) > 0L) {
      for (var in cat_vars) {
        cat_probs <- categorical_statistics(
          values = train_x[[var]],
          event_weights = event_weights,
          survivor_index = survivor_index,
          laplace = object$laplace
        )
        new_values <- as.character(new_x[[var]])
        log_surv <- log_surv + log(clip_prob(unname(cat_probs$survivor[new_values]), object$eps))
        log_event <- log_event + log(clip_prob(unname(cat_probs$event[new_values]), object$eps))
      }
    }

    surv_probs[, j] <- logistic_from_logs(log_surv, log_event)
  }

  surv_probs <- enforce_monotone_survival(surv_probs)
  if (type == "event") {
    return(1 - surv_probs)
  }
  surv_probs
}

predictSurvProb <- function(object, ...) {
  UseMethod("predictSurvProb")
}

predictSurvProb.nbsurv <- function(object, newdata, times, ...) {
  predict.nbsurv(object, newdata = newdata, times = times, type = "survival", ...)
}

print.nbsurv <- function(x, ...) {
  cat("nbsurv conditional naive Bayes model\n")
  cat("Formula: ")
  print(x$formula)
  cat("Training rows:", nrow(x$training_data), "\n")
  cat("Predictors:", paste(names(x$feature_info), collapse = ", "), "\n")
  cat("Prediction grid size:", length(x$time_grid), "\n")
  invisible(x)
}

plot.nbsurv <- function(x, times = NULL, n_curves = 5, ...) {
  if (is.null(times)) {
    times <- x$time_grid
  }
  if (!length(times)) {
    stop("No times available for plotting.")
  }

  curves <- predict.nbsurv(
    x,
    newdata = utils::head(x$training_data, n_curves),
    times = times,
    type = "survival"
  )

  graphics::matplot(
    x = times,
    y = t(curves),
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "Time",
    ylab = "Predicted survival probability",
    ...
  )

  legend_labels <- rownames(curves)
  if (is.null(legend_labels)) {
    legend_labels <- paste("obs", seq_len(nrow(curves)))
  }

  graphics::legend(
    "topright",
    legend = legend_labels,
    col = seq_len(nrow(curves)),
    lty = 1,
    lwd = 2,
    bty = "n"
  )
  invisible(curves)
}

validate_surv_response <- function(y) {
  if (!inherits(y, "Surv")) {
    stop("The left-hand side of `formula` must be a survival::Surv object.")
  }
  if (attr(y, "type") != "right") {
    stop("Only right-censored survival outcomes are supported.")
  }
}

classify_feature <- function(x) {
  if (is.numeric(x) || is.integer(x)) {
    return("continuous")
  }
  "categorical"
}

coerce_features <- function(x, feature_info) {
  for (name in names(feature_info)) {
    if (feature_info[[name]] == "continuous") {
      x[[name]] <- as.numeric(x[[name]])
    } else {
      x[[name]] <- factor(x[[name]])
    }
  }
  x
}

scale_continuous <- function(x, means, sds) {
  scaled <- sweep(as.matrix(x), 2, means, FUN = "-")
  scaled <- sweep(scaled, 2, sds, FUN = "/")
  as.data.frame(scaled, stringsAsFactors = FALSE)
}

floor_sd <- function(sd_vals, min_sd) {
  pmax(sd_vals, min_sd)
}

clip_prob <- function(x, eps) {
  x[is.na(x)] <- eps
  pmin(pmax(x, eps), 1 - eps)
}

normalise_time_grid <- function(time_grid, times, status) {
  if (is.null(time_grid)) {
    time_grid <- sort(unique(times[status == 1]))
  } else {
    time_grid <- sort(unique(as.numeric(unlist(time_grid, use.names = FALSE))))
  }

  if (!length(time_grid)) {
    stop("No event times are available to define the prediction grid.")
  }

  time_grid
}

event_case_weights <- function(object, horizon) {
  at_risk_censor <- clip_prob(object$censor_survival(pmin(object$times, horizon)), object$eps)
  as.numeric(object$times < horizon & object$status == 1) / at_risk_censor
}

continuous_statistics <- function(train_x, event_weights, survivor_index, min_sd) {
  if (!any(event_weights > 0) || !any(survivor_index)) {
    return(NULL)
  }

  train_x <- as.matrix(train_x)
  event_mean <- colSums(train_x * event_weights) / sum(event_weights)
  event_second <- colSums((train_x ^ 2) * event_weights) / sum(event_weights)
  event_sd <- floor_sd(sqrt(pmax(event_second - event_mean ^ 2, 0)), min_sd)

  survivor_mean <- colMeans(train_x[survivor_index, , drop = FALSE])
  survivor_second <- colMeans(train_x[survivor_index, , drop = FALSE] ^ 2)
  survivor_sd <- floor_sd(sqrt(pmax(survivor_second - survivor_mean ^ 2, 0)), min_sd)

  list(
    event_mean = event_mean,
    event_sd = event_sd,
    survivor_mean = survivor_mean,
    survivor_sd = survivor_sd
  )
}

categorical_statistics <- function(values, event_weights, survivor_index, laplace) {
  levels_all <- levels(values)

  survivor_counts <- stats::setNames(rep(0, length(levels_all)), levels_all)
  survivor_tab <- table(values[survivor_index])
  survivor_counts[names(survivor_tab)] <- as.numeric(survivor_tab)

  event_counts <- stats::setNames(rep(0, length(levels_all)), levels_all)
  event_split <- split(event_weights, values)
  event_counts[names(event_split)] <- vapply(event_split, sum, numeric(1))

  list(
    survivor = laplace_smoothing(survivor_counts, laplace),
    event = laplace_smoothing(event_counts, laplace)
  )
}

laplace_smoothing <- function(counts, laplace) {
  (counts + laplace) / (sum(counts) + laplace * length(counts))
}

logistic_from_logs <- function(log_surv, log_event) {
  margin <- log_event - log_surv
  probs <- 1 / (1 + exp(margin))
  pmin(pmax(probs, 0), 1)
}

enforce_monotone_survival <- function(x) {
  if (ncol(x) <= 1L) {
    return(x)
  }
  t(apply(x, 1, cummin))
}

validate_prediction_times <- function(times) {
  if (missing(times) || is.null(times) || !length(times)) {
    stop("`times` must contain at least one prediction horizon.")
  }
  times <- as.numeric(times)
  if (anyNA(times) || any(times < 0)) {
    stop("`times` must be non-missing and non-negative.")
  }
  sort(times)
}

prepare_newdata <- function(object, newdata) {
  terms_x <- stats::delete.response(object$terms)
  mf <- stats::model.frame(
    terms_x,
    data = newdata,
    na.action = stats::na.fail,
    xlev = object$xlevels
  )
  mf <- coerce_features(mf, object$feature_info)

  cont_vars <- names(object$feature_info[object$feature_info == "continuous"])
  if (!is.null(object$scaling) && length(cont_vars) > 0L) {
    mf[, cont_vars] <- scale_continuous(
      mf[, cont_vars, drop = FALSE],
      means = object$scaling$means,
      sds = object$scaling$sds
    )
  }

  mf
}

normalise_eval_times <- function(times, y) {
  if (!is.null(times)) {
    return(validate_prediction_times(times))
  }

  event_times <- y[y[, "status"] == 1, "time"]
  probs <- c(0.25, 0.5, 0.75)
  unique(as.numeric(stats::quantile(event_times, probs = probs, names = FALSE, type = 2)))
}

validate_folds <- function(folds, n) {
  if (length(folds) != 1L || is.na(folds) || folds < 2 || folds > n) {
    stop("`folds` must be an integer between 2 and the number of rows.")
  }
}

make_fold_ids <- function(n, folds, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  sample(rep(seq_len(folds), length.out = n))
}

km_stepfun <- function(time, status) {
  fit <- survival::survfit(survival::Surv(time, status) ~ 1)
  stats::stepfun(fit$time, c(1, fit$surv))
}

ipcw_brier_score <- function(time, status, pred_event, horizon, eps) {
  censor_survival <- km_stepfun(time, 1 - status)
  g_time <- clip_prob(censor_survival(time), eps)
  g_horizon <- clip_prob(censor_survival(horizon), eps)

  outcome <- as.integer(time <= horizon & status == 1)
  weight <- ifelse(
    time <= horizon & status == 1,
    1 / g_time,
    ifelse(time > horizon, 1 / g_horizon, 0)
  )

  sum(weight * (outcome - pred_event) ^ 2) / sum(weight)
}

concordance_at_horizon <- function(time, status, risk) {
  out <- tryCatch(
    survival::concordance(survival::Surv(time, status) ~ risk)$concordance,
    error = function(e) NA_real_
  )
  as.numeric(out)
}

integrated_brier_score <- function(brier, times) {
  ord <- order(times)
  t <- times[ord]
  b <- brier[ord]
  tau <- t[length(t)] - t[1L]
  if (tau <= 0) return(NA_real_)
  sum(diff(t) * (b[-length(b)] + b[-1L]) / 2) / tau
}

validate_param_grid <- function(param_grid) {
  if (!is.data.frame(param_grid) || !nrow(param_grid)) {
    stop("`param_grid` must be a non-empty data.frame.")
  }

  if (!"time_grid" %in% names(param_grid)) {
    param_grid$time_grid <- I(replicate(nrow(param_grid), NULL, simplify = FALSE))
  }

  required <- c("scale", "laplace", "min_sd", "time_grid")
  missing <- setdiff(required, names(param_grid))
  if (length(missing)) {
    stop("`param_grid` must contain columns: ", paste(required, collapse = ", "), ".")
  }

  invisible(param_grid)
}

varimp_nbsurv <- function(object,
                          newdata,
                          times,
                          metric = c("brier", "concordance"),
                          n_repeats = 10,
                          seed = NULL) {
  metric <- match.arg(metric)
  if (!is.numeric(n_repeats) || length(n_repeats) != 1L || n_repeats < 1L) {
    stop("`n_repeats` must be a positive integer.")
  }
  n_repeats <- as.integer(n_repeats)
  if (!is.null(seed)) {
    set.seed(seed)
  }

  feature_names <- names(object$feature_info)
  n <- nrow(newdata)

  baseline <- evaluate_nbsurv(object, newdata = newdata, times = times, metrics = metric)
  baseline_value <- mean(baseline[[metric]])

  importance <- vapply(feature_names, function(feature) {
    permuted_values <- vapply(seq_len(n_repeats), function(rep) {
      permuted <- newdata
      permuted[[feature]] <- permuted[[feature]][sample.int(n)]
      result <- tryCatch(
        evaluate_nbsurv(object, newdata = permuted, times = times, metrics = metric),
        error = function(e) NULL
      )
      if (is.null(result)) NA_real_ else mean(result[[metric]])
    }, numeric(1))
    mean(permuted_values, na.rm = TRUE)
  }, numeric(1))

  score <- if (metric == "brier") {
    importance - baseline_value
  } else {
    baseline_value - importance
  }

  out <- data.frame(
    feature = feature_names,
    baseline = baseline_value,
    permuted = as.numeric(importance),
    importance = as.numeric(score),
    row.names = NULL
  )
  out <- out[order(out$importance, decreasing = TRUE), , drop = FALSE]
  attr(out, "metric") <- metric
  attr(out, "n_repeats") <- n_repeats
  class(out) <- c("varimp_nbsurv", "data.frame")
  out
}

print.varimp_nbsurv <- function(x, ...) {
  cat("nbsurv permutation variable importance\n")
  cat("Metric:", attr(x, "metric"), "| repeats:", attr(x, "n_repeats"), "\n\n")
  print.data.frame(x, row.names = FALSE)
  invisible(x)
}

plot.varimp_nbsurv <- function(x, ...) {
  ord <- order(x$importance)
  graphics::barplot(
    height = x$importance[ord],
    names.arg = x$feature[ord],
    horiz = TRUE,
    las = 1,
    xlab = paste0("Importance (", attr(x, "metric"), " degradation)"),
    ...
  )
  invisible(x)
}

plot.cv_nbsurv <- function(x, metric = c("brier", "concordance"), ...) {
  metric <- match.arg(metric)
  graphics::matplot(
    x = x$summary$time,
    y = x$summary[[metric]],
    type = "b",
    pch = 19,
    lty = 1,
    lwd = 2,
    xlab = "Time",
    ylab = paste("Mean", metric, "across folds"),
    ...
  )
  invisible(x)
}

plot.tune_nbsurv <- function(x, ...) {
  values <- x$results$mean_metric
  cols <- rep("black", length(values))
  cols[x$best_index] <- "firebrick"
  graphics::plot(
    x = seq_along(values),
    y = values,
    pch = 19,
    col = cols,
    xlab = "Candidate",
    ylab = paste("Mean", x$metric),
    ...
  )
  graphics::points(x$best_index, values[x$best_index], pch = 1, cex = 2, col = "firebrick")
  invisible(x)
}
