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

