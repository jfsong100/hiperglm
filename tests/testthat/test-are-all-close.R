
test_that("function are_all_close correctly returns TRUE", {
  v <- rep(1, 10)
  u <- rep(1, 10)
  expect_true(are_all_close(
    v, u,
    abs_tol = 1e-2, rel_tol = 1e-2
  ))
})


test_that("function are_all_close correctly returns FALSE because the relative error is above rel_tol", {
  v <- 3
  u <- v + 0.1
  expect_false(are_all_close(
    v, u,
    abs_tol = 1e-2, rel_tol = 1e-16
  ))
})



test_that("function are_all_close correctly returns FALSE because the absolute error is above abs_tol", {
  v <- 1
  u <- v + 10
  expect_false(are_all_close(
    v, u,
    abs_tol = 1e-6, rel_tol = 1e-6
  ))
})
