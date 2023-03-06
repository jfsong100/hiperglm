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

linear_BFGS <- function(design, outcome, noise_var = 1) {
  initial <- rep(0, ncol(design))
  MLE <- stats::optim(initial, linear_loglik, linear_gradient,
    design = design, outcome = outcome, noise_var = 1,
    control = list(fnscale = -1),
    method = "BFGS"
  )$par
  return(MLE)
}

logit_loglik <- function(coef, design, outcome) {

  n <- length(outcome)
  log_lik <- sum(sapply(1:n, function(i){
    outcome[i] * sum(design[i,] * coef) -
      log(1 + exp(sum(design[i,] * coef)))
    }))
  return(log_lik)

}

logit_gradient <- function(coef, design, outcome) {

  pi <- exp(design %*% coef) / (1 + exp(design %*% coef))
  gradient <- as.vector(t(outcome - pi) %*% design)

  return(gradient)

}

logit_BFGS <- function(design, outcome) {
  initial <- rep(0, ncol(design))
  MLE <- stats::optim(initial, logit_loglik, logit_gradient,
                      design = design, outcome = outcome,
                      control = list(fnscale = -1),
                      method = "BFGS"
  )$par
  return(MLE)
}


logit_hessian <- function(coef, design, outcome){
  eta <- design %*% coef
  pi <- exp(eta)/(1+exp(eta))
  W <- diag(c(pi * (1 - pi)))
  hessian <- t(design) %*% W %*% design
  return(hessian)
}

logit_newton <- function(design, outcome){

  max_iteration <- 1000
  t <- 0

  threshold <- ncol(design)/max_iteration
  diff <- 1

  coef_old <- rep(0, ncol(design))

  loglik_old <- logit_loglik(coef_old, design, outcome)

  while(diff >= threshold & t <= max_iteration ){

  eta_old <- design %*% coef_old
  pi_old <- exp(eta_old)/(1+exp(eta_old))
  hessian_old <- logit_hessian(coef_old, design, outcome)

  coef_update <-  coef_old + solve(hessian_old, t(design) %*% (outcome-pi_old))
  loglik_update <- logit_loglik(coef_update, design, outcome)
  diff <- abs(loglik_update - loglik_old)

  loglik_old <- loglik_update
  coef_old <- coef_update
  t <- t+1

  }

  return(as.vector(coef_update))

}



