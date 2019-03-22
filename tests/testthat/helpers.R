# test functions

set.seed(123)

# need to get tf object here
library (tensorflow)

expect_ok <- function (expr)
  expect_error(expr, NA)

# check a greta operation and the equivalent R operation give the same output
# e.g. check_op(sum, randn(100, 3))
check_covariance <- function (kernel, X, X_prime = NULL, expected, tol = 1e-6) {

  tf$reset_default_graph()

  if (is.null(X_prime))
    X_prime <- X
  
  greta_out <- greta::calculate(kernel(X, X_prime))
  difference <- as.vector(abs(expected - greta_out))
  expect_true(all(difference < tol))

}
