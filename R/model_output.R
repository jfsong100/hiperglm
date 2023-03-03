#' @export
coef.hglm <- function(hglm_out) {
  as.vector(hglm_out$coef)
}

#' @export
vcov.hglm <- function(hglm_out) {
  warning("To be implemented.")
}

#' @export
print.hglm <- function(hglm_out) {
  cat("Output of hiper_glm.")
}
