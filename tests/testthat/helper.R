simulate_data <- function(
    n_obs, n_pred, model = "linear", intercept = NULL,
    coef_true = NULL, design = NULL, seed = NULL, option = list()
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if ((model != "linear")  && !is.null(option$signal_to_noise)) {
    warning(paste(
      "The `signal_to_noise` option is currently unsupported for",
      "non-linear models and will be ignored."
    ))
  }
  if (is.null(coef_true)) {
    coef_true <- rnorm(n_pred, sd = 1 / sqrt(n_pred))
  }
  if (is.null(design)) {
    design <- matrix(rnorm(n_obs * n_pred), nrow = n_obs, ncol = n_pred)
  }
  if (!is.null(intercept)) {
    if (!is.numeric(intercept)) {
      stop("The intercept argument must be numeric.")
    }
    coef_true <- c(intercept, coef_true)
    design <- cbind(rep(1, n_obs), design)
  }
  expected_mean <- as.vector(design %*% coef_true)
  if (model == 'linear') {
    signal_to_noise <- option$signal_to_noise
    if (is.null(signal_to_noise)) {
      signal_to_noise <- 0.1
    }
    noise_magnitude <- sqrt(var(expected_mean) / signal_to_noise^2)
    noise <- noise_magnitude * rnorm(n_obs)
    outcome <- expected_mean + noise
  } else {
    n_trial <- option$n_trial
    prob <- 1 / (1 + exp(-expected_mean))
    # Object type of `outcome` returned by this function is variable. One should
    # in general be careful about introducing this type of inconsistency, but
    # sometimes one might find it the most natural and/or reasonable thing to do.
    if (is.null(n_trial)) {
      outcome <- rbinom(n_obs, 1, prob)
    } else {
      n_success <- rbinom(n_obs, n_trial, prob)
      outcome <- list(n_success = n_success, n_trial = n_trial)
    }
  }
  return(list(design = design, outcome = outcome, coef_true = coef_true))
}

approx_grad <- function(func, x, design, outcome, dx = .Machine$double.eps^(1 / 3)) {
  numerical_grad <- rep(0, length(x))

  for (i in 1:length(x)) {
    add_x <- x
    add_x[i] <- add_x[i] + dx
    minus_x <- x
    minus_x[i] <- minus_x[i] - dx
    add_term <- func(add_x, design, outcome)
    minus_term <- func(minus_x, design, outcome)
    numerical_grad[i] <- (add_term - minus_term) / 2 / dx
  }

  return(numerical_grad)
}
