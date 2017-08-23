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

