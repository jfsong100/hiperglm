#' @export
hiper_glm <- function(design, outcome, model = "linear", option = list()) {
  supported_model <- c("linear", "logit")
  if (!(model %in% supported_model)) {
    stop(sprintf("Model is not available"))
  }
  # TODO: find MLE.
  if (model == "linear") {
    if ((is.null(option$mle_solver)) || (option$mle_solver == "pseudo inverse")) {
      MLE <- MLE_pseudo_inverse(design, outcome)
    } else if (option$mle_solver == "BFGS") {
      MLE <- MLE_BFGS_function(linear_loglik, linear_gradient, design, outcome, noise_var = 1)
    }
  }else{
    if ((is.null(option$mle_solver)) ) {
      MLE <- logit_newton(design, outcome, option)
    } else if (option$mle_solver == "BFGS") {
      MLE <- MLE_BFGS_function(logit_loglik, logit_gradient, design, outcome)
    }
  }

  hglm_out <- list(coef = MLE)
  class(hglm_out) <- "hglm"
  return(hglm_out)
}
