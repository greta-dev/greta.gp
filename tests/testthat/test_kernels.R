context('kernels')

test_that('base kernels evaluate self-covariance correctly', {

  source("helpers.R")
  skip_if_not(greta:::check_tf_version())
  skip_if_not(gpflowr::gpflow_available())

  n <- 5
  x_ <- rnorm(n)
  x <- as_data(x_)

  per <- runif(1)
  var <- runif(1)
  len <- runif(1)

  r <- as.matrix(dist(x_ / len))

  # bias
  check_covariance(bias(var),
                   x,
                   expected = matrix(var, n, n))

  # white
  check_covariance(white(var),
                   x,
                   expected = diag(n) * var)

  # linear
  check_covariance(linear(var),
                   x,
                   expected = outer(x_, x_) * var)

  # rbf
  check_covariance(rbf(len, var),
                   x,
                   expected = var * exp(-0.5 * r ^ 2))

  # expo
  check_covariance(expo(len, var),
                   x,
                   expected = var * exp(-0.5 * r))

  # mat12
  check_covariance(mat12(len, var),
                   x,
                   expected = var * exp(-r))

  # mat32
  check_covariance(mat32(len, var),
                   x,
                   expected = var * (1 + sqrt(3) * r) * exp(-sqrt(3) * r))

  # mat52
  check_covariance(mat52(len, var),
                   x,
                   expected = var * (1 + sqrt(5) * r + 5 / 3 * r ^ 2) *
                     exp(-sqrt(5) * r))

  # periodic
  r <- (pi * as.matrix(dist(x_))) / per
  r <- sin(r) / len
  check_covariance(periodic(per, len, var),
                   x,
                   expected = var * exp(-0.5 * r ^ 2), tol = 1e-2)

})

test_that('compound kernels evaluate self-covariance correctly', {

  source("helpers.R")
  skip_if_not(greta:::check_tf_version())
  skip_if_not(gpflowr::gpflow_available())

  n <- 5
  x_ <- rnorm(n)
  x <- as_data(x_)

  var <- runif(1)
  len <- runif(1)

  r <- as.matrix(dist(x_ / len))

  # additive
  check_covariance(linear(var) + rbf(len, var),
                   x,
                   expected = (outer(x_, x_) * var) + (var * exp(-0.5 * r ^ 2)))

  # multiplicative
  check_covariance(linear(var) * rbf(len, var),
                   x,
                   expected = (outer(x_, x_) * var) * (var * exp(-0.5 * r ^ 2)))

})

test_that('compound kernels can act on specific dimensions', {

  source("helpers.R")
  skip_if_not(greta:::check_tf_version())
  skip_if_not(gpflowr::gpflow_available())

  n <- 5
  x_ <- cbind(rnorm(n), runif(n))
  x <- as_data(x_)

  var <- runif(1)
  len <- runif(1)

  x_1 <- x_[, 1]
  r_2 <- as.matrix(dist(x_[, 2] / len))

  # additive
  check_covariance(linear(var, columns = 1) + rbf(len, var, columns = 2),
                   x,
                   expected = (outer(x_1, x_1) * var) + (var * exp(-0.5 * r_2 ^ 2)))

  # multiplicative
  check_covariance(linear(var, columns = 1) * rbf(len, var, columns = 2),
                   x,
                   expected = (outer(x_1, x_1) * var) * (var * exp(-0.5 * r_2 ^ 2)))

})

test_that('kernels error on badly shaped inputs', {

  source("helpers.R")
  skip_if_not(greta:::check_tf_version())
  skip_if_not(gpflowr::gpflow_available())

  kernel <- rbf(1, 1)

  bad_x <- greta_array(1:24, dim = c(2, 3, 4))
  x1 <- greta_array(1:10)
  x2 <- greta_array(1:10, dim = c(2, 5))

  expect_error(kernel(bad_x),
               "X must be a 2D greta array")
  expect_error(kernel(x1, x2),
               "number of columns of X and X_prime do not match")

})

test_that('kernel constructors error on bad columns', {

  source("helpers.R")
  skip_if_not(greta:::check_tf_version())
  skip_if_not(gpflowr::gpflow_available())

  expect_error(rbf(1, 1, columns = 1:2),
               "columns has length 2 but the kernel has dimension 1")

  expect_error(rbf(1, 1, columns = -1),
               "columns must be a vector of positive integers, but was -1")

})

test_that('kernels error if combined with other things', {

  source("helpers.R")
  skip_if_not(greta:::check_tf_version())
  skip_if_not(gpflowr::gpflow_available())

  expect_error(bias(1) + 1,
               "can only combine a greta kernel with another greta kernel")
  expect_error(1 + bias(1),
               "can only combine a greta kernel with another greta kernel")

  expect_error(bias(1) * 1,
               "can only combine a greta kernel with another greta kernel")
  expect_error(1 * bias(1),
               "can only combine a greta kernel with another greta kernel")

})

test_that('kernels print their own names', {

  source("helpers.R")
  skip_if_not(greta:::check_tf_version())
  skip_if_not(gpflowr::gpflow_available())

  expect_output(print(bias(1)),
                "bias kernel")

  expect_output(print(white(1)),
                "white kernel")

  expect_output(print(linear(1)),
                "linear kernel")

  expect_output(print(rbf(1, 1)),
                "radial basis kernel")

  expect_output(print(expo(1, 1)),
                "exponential kernel")

  expect_output(print(mat12(1, 1)),
                "Matern 1/2 kernel")

  expect_output(print(mat32(1, 1)),
                "Matern 3/2 kernel")

  expect_output(print(mat52(1, 1)),
                "Matern 5/2 kernel")

  expect_output(print(periodic(1, 1, 1)),
                "periodic kernel")

  expect_output(print(bias(1) + bias(1)),
                "additive kernel")

  expect_output(print(bias(1) * bias(1)),
                "multiplicative kernel")

})



