test_that("linalg and optim least-sq coincide", {
  n_obs <- 32
  n_pred <- 4
  data <- simulate_data(n_obs, n_pred, model = "linear", seed = 1918)
  design <- data$design
  outcome <- data$outcome
  via_linalg_out <- hiper_glm(design, outcome, model = "linear")
  via_bfgs_out <- hiper_glm(
    design, outcome,
    model = "linear", option = list(mle_solver = "BFGS")
  )
  expect_true(are_all_close(
    coef(via_linalg_out), coef(via_bfgs_out),
    abs_tol = 1e-6, rel_tol = 1e-6
  ))
})

test_that("gradient calculation by comparing it against a numerical one linear", {
  n_obs <- 32
  n_pred <- 4
  data <- simulate_data(n_obs, n_pred, model = "linear", seed = 1918)
  design <- data$design
  outcome <- data$outcome
  test_MLE <- rep(1, n_pred)
  expect_true(are_all_close(
    linear_gradient(test_MLE, design, outcome, noise_var = 1),
    approx_grad(linear_loglik, test_MLE, design, outcome),
    abs_tol = 1e-2, rel_tol = 1e-2
  ))
})

test_that("newton and bfgs outputs coincide on logit model", {
  n_obs <- 32; n_pred <- 4
  data <- simulate_data(n_obs, n_pred, model = "logit", seed = 1918)
  design <- data$design; outcome <- data$outcome
  via_newton_out <- hiper_glm(design, outcome, model = "logit")
  via_bfgs_out <- hiper_glm(
    design, outcome, model = "logit", option = list(mle_solver = "BFGS")
  )
  expect_true(are_all_close(
    coef(via_newton_out), coef(via_bfgs_out), abs_tol = 1e-2, rel_tol = 1e-2
  ))
})

test_that("gradient calculation by comparing it against a numerical one logit", {
  n_obs <- 32
  n_pred <- 4
  data <- simulate_data(n_obs, n_pred, model = "logit", seed = 1918)
  design <- data$design
  outcome <- data$outcome
  test_MLE <- rep(1, n_pred)
  expect_true(are_all_close(
    logit_gradient(test_MLE, design, outcome),
    approx_grad(logit_loglik, test_MLE, design, outcome),
    abs_tol = 1e-6, rel_tol = 1e-6
  ))
})

test_that("using QR decomposition and use LU decompisition", {
  n_obs <- 32; n_pred <- 4
  data <- simulate_data(n_obs, n_pred, model = "logit", seed = 1918)
  design <- data$design; outcome <- data$outcome
  via_newton_out_qr <- hiper_glm(design, outcome, model = "logit")
  via_newton_out_lu <- hiper_glm(design, outcome, model = "logit", option= list(solver="lu"))
  expect_true(are_all_close(
    coef(via_newton_out_qr), coef(via_newton_out_lu), abs_tol = 1e-2, rel_tol = 1e-2
  ))
})
