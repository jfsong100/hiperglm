MLE_pseudo_inverse <- function(design, outcome) {
  return(solve(t(design) %*% design, t(design) %*% outcome))
}

linear_loglik <- function(coef, design, outcome, noise_var = 1) {
  return(-sum((outcome - design %*% coef)^2) / 2 / noise_var)
}

linear_gradient <- function(coef, design, outcome, noise_var = 1) {
  gradient <- -(t(design) %*% design %*% coef - t(design) %*% outcome) / noise_var
  return(gradient)
}

MLE_BFGS_function <- function(loglik_function, gradient_function, design, outcome, ...) {
  initial <- rep(0, ncol(design))
  BFGS_MLE <- stats::optim(initial, loglik_function, gradient_function,
    design = design, outcome = outcome, control = list(fnscale = -1),
    method = "BFGS",...
  )$par
  return(BFGS_MLE)
}

logit_loglik <- function(coef, design, outcome) {

  n <- length(outcome)
  prob <- expit_function(design %*% coef)
  log_lik <- sum(outcome * log(prob) + (1 - outcome) * log(1 - prob))
  return(log_lik)

}

logit_gradient <- function(coef, design, outcome) {

  prob <- expit_function(design %*% coef)
  gradient <- as.vector(t(outcome - prob) %*% design)

  return(gradient)

}


logit_hessian <- function(coef, design, outcome){

  prob <-  expit_function(design %*% coef)
  W <- diag(c(prob * (1 - prob)))
  hessian <- t(design) %*% W %*% design
  return(hessian)
}

logit_newton <- function(design, outcome){

  max_iteration <- 1000
  iter <- 0

  stop_indicator <- FALSE

  coef_old <- rep(0, ncol(design))

  loglik_old <- logit_loglik(coef_old, design, outcome)

  while(!stop_indicator & iter <= max_iteration ){

  prob_old <-  expit_function(design %*% coef_old)
  hessian_old <- logit_hessian(coef_old, design, outcome)

  coef_update <-  coef_old + solve(hessian_old, t(design) %*% (outcome-prob_old))
  loglik_update <- logit_loglik(coef_update, design, outcome)
  stop_indicator <- are_all_close(loglik_update, loglik_old)

  loglik_old <- loglik_update
  coef_old <- coef_update
  iter <- iter+1

  }
  if(iter == max_iteration){
    warning("The max iterations are reached!")
  }
  return(as.vector(coef_update))

}


expit_function <- function(x){
  return(exp(x)/(1+exp(x)))
}

are_all_close <- function(v, w, abs_tol = 1e-6, rel_tol = 1e-6) {
  abs_diff <- abs(v - w)
  are_all_within_atol <- all(abs_diff < abs_tol)
  are_all_within_rtol <- all(abs_diff < rel_tol * pmax(abs(v), abs(w)))
  return(are_all_within_atol && are_all_within_rtol)
}
