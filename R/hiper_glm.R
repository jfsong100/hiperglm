#' @export
hiper_glm <- function(design, outcome, model = "linear", option = list()) {
  supported_model <- c("linear")
  if (!(model %in% supported_model)) {
    stop(sprintf("Model is not available"))
  }
  # TODO: find MLE.
  if (model == "linear") {
    if ((is.null(option$mle_solver) == TRUE) || (option$mle_solver == "pseudo inverse")) {
      MLE <- MLE_pseudo_inverse(design, outcome)
    } else if (option$mle_solver == "BFGS") {
      MLE <- linear_BFGS(design, outcome, noise_var = 1)
    }
  }

  hglm_out <- list(coef = MLE)
  class(hglm_out) <- "hglm"
  return(hglm_out)
}
