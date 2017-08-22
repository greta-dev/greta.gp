# test functions

# need to get tf object here
library (tensorflow)

expect_ok <- function (expr)
  expect_error(expr, NA)

# evaluate a greta_array, node, or tensor
grab <- function (x) {

  if (inherits(x, "node"))
    x <- as.greta_array(x)

  if (inherits(x, "greta_array")) {
    dag <- dag_class$new(list(x))
    x$node$define_tf(dag)
    x <- get(dag$tf_name(x$node),
             envir = dag$tf_environment)
  }

  tf$Session()$run(x)

}

# evaluate the (unadjusted) density of distribution greta array at some data
get_covariance <- function (kernel, X, X_prime = NULL) {

  x <- as_data(X)

  if (is.null(X_prime)) {
    K <- kernel(x)
  } else {
    xp <- as_data(X_prime)
    K <- kernel(x, xp)
  }

  # create dag and define the density
  dag <- greta:::dag_class$new(list(K))
  K$node$define_tf(dag)

  # get the log density as a vector
  tensor_name <- dag$tf_name(K$node)
  tensor <- get(tensor_name, envir = dag$tf_environment)
  grab(tensor)

}

# check a greta operation and the equivalent R operation give the same output
# e.g. check_op(sum, randn(100, 3))
check_covariance <- function (kernel, X, X_prime = NULL, expected, tol = 1e-6) {

  tf$reset_default_graph()

  greta_out <- get_covariance(kernel, X, X_prime)
  difference <- as.vector(abs(expected - greta_out))
  expect_true(all(difference < tol))
}

