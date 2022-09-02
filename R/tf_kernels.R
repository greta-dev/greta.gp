# tensorflow implementations of common kernels

# given a vector of column numbers of relevant dimensions (in 0-based indexing),
# pull out the corresponding sub-matrix
tf_cols <- function(X, active_dims) {
  X[, , active_dims, drop = FALSE]
}

tf_distance <- function(x1, x2, squared = FALSE) {
  n1 <- dim(x1)[[2]]
  n2 <- dim(x2)[[2]]

  x1 <- tf$tile(tf$expand_dims(x1, 3L), list(1L, 1L, 1L, n2))
  x2 <- tf$transpose(x2, perm = c(0L, 2L, 1L))
  x2 <- tf$tile(tf$expand_dims(x2, 1L), list(1L, n1, 1L, 1L))

  dists <- (x1 - x2)^2
  dist <- tf$reduce_sum(dists, axis = 2L)

  if (!squared) {
    dist <- tf$math$sqrt(dist)
  }

  dist
}

# build a matrix with dimension given by the number of rows in X and the
# number of rows in X_prime, filled with the given *constant* value
tf_empty_along <- function(X, X_prime = NULL, fill = 1) {
  if (is.null(X_prime)) {
    dims_out <- tf$stack(c(tf$shape(X)[0], dim(X)[[2]]))
  } else {
    dims_out <- tf$stack(c(tf$shape(X)[0], dim(X)[[2]], dim(X_prime)[[2]]))
  }

  switch(as.character(fill),
    `1` = tf$ones(dims_out, dtype = tf_float()),
    `0` = tf$zeros(dims_out, dtype = tf_float())
  )
}


# bias (or constant) kernel
# k(x, y) = \sigma^2
tf_bias <- function(X, X_prime, variance, active_dims) {

  # create and return covariance matrix
  tf_empty_along(X, X_prime, 1) * variance
}

# white kernel
# diagonal with specified variance if self-kernel, all 0s otherwise
tf_white <- function(X, X_prime, variance, active_dims) {

  # only non-zero for self-covariance matrices
  if (identical(X, X_prime)) {
    variance <- tf$squeeze(variance, 2L)
    d <- tf_empty_along(X, X_prime = NULL, fill = 1) * variance
    d <- tf$linalg$diag(d)
  } else {
    d <- tf_empty_along(X, X_prime, 0)
  }

  # return constructed covariance matrix
  d
}

tf_iid <- function(X, X_prime, variance, active_dims) {

  # pull out active dimensions
  X <- tf_cols(X, active_dims)
  X_prime <- tf_cols(X_prime, active_dims)

  # find where these values match and assign the variance as a covariance there
  # (else set it to 0)
  distance <- tf_distance(X, X_prime)
  tf_as_float(distance < fl(1e-12)) * variance
}

# squared exponential kernel (RBF)
tf_rbf <- function(X, X_prime, lengthscales, variance, active_dims) {

  # pull out active dimensions
  X <- tf_cols(X, active_dims)
  X_prime <- tf_cols(X_prime, active_dims)

  # calculate squared distances
  r2 <- squared_dist(X, X_prime, lengthscales)

  # construct and return RBF kernel
  variance * tf$math$exp(-r2 / fl(2))
}

# rational_quadratic kernel
tf_rational_quadratic <- function(X, X_prime, lengthscales, variance, alpha, active_dims) {

  # pull out active dimensions
  X <- tf_cols(X, active_dims)
  X_prime <- tf_cols(X_prime, active_dims)

  # calculate squared distances (scaled if needed)
  r2 <- squared_dist(X, X_prime, lengthscales)

  # construct and return rational quadratic kernel
  variance * (fl(1) + r2 / (fl(2) * alpha))^-alpha
}

# linear kernel (base class)
tf_linear <- function(X, X_prime, variances, active_dims) {

  # pull out active dimensions
  X <- tf_cols(X, active_dims)
  X_prime <- tf_cols(X_prime, active_dims)

  # full kernel
  tf$linalg$matmul(
    tf$math$multiply(variances, X),
    X_prime,
    transpose_b = TRUE
  )
}

tf_polynomial <- function(X, X_prime, variances, offset, degree, active_dims) {

  # pull out active dimensions
  X <- tf_cols(X, active_dims)
  X_prime <- tf_cols(X_prime, active_dims)

  # full kernel
  tf$math$pow(
    tf$linalg$matmul(
      tf$math$multiply(variances, X),
      X_prime,
      transpose_b = TRUE
    ) + offset,
    degree
  )
}

# # exponential kernel (stationary class)
tf_exponential <- function(X, X_prime, lengthscales, variance, active_dims) {

  # pull out active dimensions
  X <- tf_cols(X, active_dims)
  X_prime <- tf_cols(X_prime, active_dims)

  # calculate squared distances (scaled if needed)
  r <- absolute_dist(X, X_prime, lengthscales)

  # construct and return exponential kernel
  variance * tf$math$exp(-fl(0.5) * r)
}

# Matern12 kernel (stationary class)
tf_Matern12 <- function(X, X_prime, lengthscales, variance, active_dims) {

  # pull out active dimensions
  X <- tf_cols(X, active_dims)
  X_prime <- tf_cols(X_prime, active_dims)

  # calculate squared distances (scaled if needed)
  r <- absolute_dist(X, X_prime, lengthscales)

  # construct and return Matern12 kernel
  variance * tf$math$exp(-r)
}

# Matern32 kernel (stationary class)
tf_Matern32 <- function(X, X_prime, lengthscales, variance, active_dims) {

  # pull out active dimensions
  X <- tf_cols(X, active_dims)
  X_prime <- tf_cols(X_prime, active_dims)

  # calculate squared distances (scaled if needed)
  r <- absolute_dist(X, X_prime, lengthscales)

  # precalculate root3
  sqrt3 <- fl(sqrt(3))

  # construct and return Matern32 kernel
  variance * (fl(1) + sqrt3 * r) * tf$math$exp(-sqrt3 * r)
}

# Matern52 kernel (stationary class)
tf_Matern52 <- function(X, X_prime, lengthscales, variance, active_dims) {

  # pull out active dimensions
  X <- tf_cols(X, active_dims)
  X_prime <- tf_cols(X_prime, active_dims)

  # calculate squared distances (scaled if needed)
  r <- absolute_dist(X, X_prime, lengthscales)

  # precalculate root5
  sqrt5 <- fl(sqrt(5))

  # construct and return Matern52 kernel
  variance * (fl(1) + sqrt5 * r + fl(5) / fl(3) * tf$math$square(r)) * tf$math$exp(-sqrt5 * r)
}

# cosine kernel (stationary class)
tf_cosine <- function(X, X_prime, lengthscales, variance, active_dims) {

  # pull out active dimensions
  X <- tf_cols(X, active_dims)
  X_prime <- tf_cols(X_prime, active_dims)

  # calculate squared distances (scaled if needed)
  r <- absolute_dist(X, X_prime, lengthscales)

  # construct and return cosine kernel
  variance * tf$math$cos(r)
}

# periodic kernel
tf_periodic <- function(X, X_prime, lengthscale, variance, period) {

  # calculate squared distances (scaled if needed)
  exp_arg <- fl(pi) * absolute_dist(X, X_prime) / period
  exp_arg <- tf$math$sin(exp_arg) / lengthscale

  # construct and return periodic kernel
  variance * tf$math$exp(-fl(0.5) * tf$math$square(exp_arg))
}

tf_Prod <- function(kernel_a, kernel_b) {
  tf$math$multiply(kernel_a, kernel_b)
}

tf_Add <- function(kernel_a, kernel_b) {
  tf$math$add(kernel_a, kernel_b)
}

# rescale, calculate, and return distance
get_dist <- function(X, X_prime, lengthscales = NULL, squared = FALSE) {
  if (!is.null(lengthscales)) {
    X <- X / lengthscales
    X_prime <- X_prime / lengthscales
  }

  tf_distance(X, X_prime, squared = squared)
}

squared_dist <- function(X, X_prime, lengthscales = NULL) {
  get_dist(X, X_prime, lengthscales, squared = TRUE)
}

absolute_dist <- function(X, X_prime, lengthscales = NULL) {
  get_dist(X, X_prime, lengthscales, squared = FALSE)
}

# combine as module for export via internals
# tf_kernels_module <- module(tf_static,
#                             tf_constant,
#                             tf_bias,
#                             tf_squared_exponential,
#                             tf_rational_quadratic,
#                             tf_linear,
#                             tf_polynomial,
#                             tf_exponential,
#                             tf_Matern12,
#                             tf_Matern32,
#                             tf_Matern52,
#                             tf_cosine)
