#' @title Gaussian process kernels
#' @name kernels
#'
#' @description Create and combine Gaussian process kernels (covariance
#'   functions) for use in Gaussian process models.
#'
#' @param variance,variances (scalar/vector) the variance of a Gaussian process
#'   prior in all dimensions (`variance`) or in each dimension
#'   (`variances`)
#' @param lengthscale,lengthscales (scalar/vector) the correlation decay
#'   distance along all dimensions (`lengthscale`) or each dimension
#'   ((`lengthscales`)) of the Gaussian process
#' @param alpha (scalar) additional parameter in rational quadratic kernel
#' @param offset (scalar) offset in polynomial kernel
#' @param degree (scalar) degree of polynomial kernel
#' @param period (scalar) the period of the Gaussian process
#' @param columns (scalar/vector integer, not a greta array) the columns of the
#'   data matrix on which this kernel acts. Must have the same dimensions as
#'   lengthscale parameters.
#'
#' @details The kernel constructor functions each return a *function* (of
#'   class `greta_kernel`) which can be executed on greta arrays to compute
#'   the covariance matrix between points in the space of the Gaussian process.
#'   The `+` and `*` operators can be used to combine kernel functions
#'   to create new kernel functions.
#'
#'   Note that `bias` and `constant` are identical names for the same
#'   underlying kernel.
#'
#'   `iid` is equivalent to `bias` where all entries in `columns`
#'   match (where the absolute euclidean distance is less than
#'   1e-12), and `white` where they don't; i.e. an independent Gaussian
#'   random effect.
#'
#' @return greta kernel with class "greta_kernel"
#'
#' @examples
#' \dontrun{
#' # create a radial basis function kernel on two dimensions
#' k1 <- rbf(lengthscales = c(0.1, 0.2), variance = 0.6)
#'
#' # evaluate it on a greta array to get the variance-covariance matrix
#' x <- greta_array(rnorm(8), dim = c(4, 2))
#' k1(x)
#'
#' # non-symmetric covariance between two sets of points
#' x2 <- greta_array(rnorm(10), dim = c(5, 2))
#' k1(x, x2)
#'
#' # create a bias kernel, with the variance as a variable
#' k2 <- bias(variance = lognormal(0, 1))
#'
#' # combine two kernels and evaluate
#' K <- k1 + k2
#' K(x, x2)
#'
#' # other kernels
#' constant(variance = lognormal(0, 1))
#' white(variance = lognormal(0, 1))
#' iid(variance = lognormal(0, 1))
#' rational_quadratic(lengthscales = c(0.1, 0.2), alpha = 0.5, variance = 0.6)
#' linear(variances = 0.1)
#' polynomial(variances = 0.6, offset = 0.8, degree = 2)
#' expo(lengthscales = 0.6, variance = 0.9)
#' mat12(lengthscales = 0.5, variance = 0.7)
#' mat32(lengthscales = 0.4, variance = 0.8)
#' mat52(lengthscales = 0.3, variance = 0.9)
#' cosine(lengthscales = 0.68, variance = 0.8)
#' periodic(period = 0.71, lengthscale = 0.59, variance = 0.2)
#' }
NULL

#' @rdname kernels
#' @export
bias <- function(variance) {
  greta_kernel("bias",
    tf_name = "tf_bias",
    parameters = list(variance = variance)
  )
}

#' @rdname kernels
#' @export
constant <- function(variance) {
  greta_kernel("constant",
    tf_name = "tf_bias",
    parameters = list(variance = variance)
  )
}

#' @rdname kernels
#' @export
white <- function(variance) {
  greta_kernel("white",
    tf_name = "tf_white",
    parameters = list(variance = variance)
  )
}

#' @rdname kernels
#' @export
iid <- function(variance, columns = 1) {
  greta_kernel("iid",
    tf_name = "tf_iid",
    parameters = list(variance = variance),
    arguments = list(
      active_dims = check_active_dims(
        columns,
        rep(1, length(columns))
      )
    )
  )
}

#' @rdname kernels
#' @export
rbf <- function(lengthscales, variance, columns = seq_along(lengthscales)) {
  greta_kernel("radial basis",
    tf_name = "tf_rbf",
    parameters = list(
      lengthscales = t(lengthscales),
      variance = variance
    ),
    arguments = list(active_dims = check_active_dims(columns, lengthscales))
  )
}

#' @rdname kernels
#' @export
rational_quadratic <- function(lengthscales, variance, alpha, columns = seq_along(lengthscales)) {
  greta_kernel("rational quadratic",
    tf_name = "tf_rational_quadratic",
    parameters = list(
      lengthscales = t(lengthscales),
      variance = variance,
      alpha = alpha
    ),
    arguments = list(active_dims = check_active_dims(columns, lengthscales))
  )
}

#' @rdname kernels
#' @export
linear <- function(variances, columns = seq_along(variances)) {
  greta_kernel("linear",
    tf_name = "tf_linear",
    parameters = list(variance = t(variances)),
    arguments = list(active_dims = check_active_dims(columns, variances))
  )
}

#' @rdname kernels
#' @export
polynomial <- function(variances, offset, degree, columns = seq_along(variances)) {
  greta_kernel("polynomial",
    tf_name = "tf_polynomial",
    parameters = list(
      variance = t(variances),
      offset = offset,
      degree = degree
    ),
    arguments = list(active_dims = check_active_dims(columns, variances))
  )
}

#' @rdname kernels
#' @export
expo <- function(lengthscales, variance, columns = seq_along(lengthscales)) {
  greta_kernel("exponential",
    tf_name = "tf_exponential",
    parameters = list(
      lengthscales = t(lengthscales),
      variance = variance
    ),
    arguments = list(active_dims = check_active_dims(columns, lengthscales))
  )
}

#' @rdname kernels
#' @export
mat12 <- function(lengthscales, variance, columns = seq_along(lengthscales)) {
  greta_kernel("Matern 1/2",
    tf_name = "tf_Matern12",
    parameters = list(
      lengthscales = t(lengthscales),
      variance = variance
    ),
    arguments = list(active_dims = check_active_dims(columns, lengthscales))
  )
}

#' @rdname kernels
#' @export
mat32 <- function(lengthscales, variance, columns = seq_along(lengthscales)) {
  greta_kernel("Matern 3/2",
    tf_name = "tf_Matern32",
    parameters = list(
      lengthscales = t(lengthscales),
      variance = variance
    ),
    arguments = list(active_dims = check_active_dims(columns, lengthscales))
  )
}

#' @rdname kernels
#' @export
mat52 <- function(lengthscales, variance, columns = seq_along(lengthscales)) {
  greta_kernel("Matern 5/2",
    tf_name = "tf_Matern52",
    parameters = list(
      lengthscales = t(lengthscales),
      variance = variance
    ),
    arguments = list(active_dims = check_active_dims(columns, lengthscales))
  )
}

#' @rdname kernels
#' @export
cosine <- function(lengthscales, variance, columns = seq_along(lengthscales)) {
  greta_kernel("cosine",
    tf_name = "tf_cosine",
    parameters = list(
      lengthscales = t(lengthscales),
      variance = variance
    ),
    arguments = list(active_dims = check_active_dims(columns, lengthscales))
  )
}

#' @rdname kernels
#' @export
periodic <- function(period, lengthscale, variance) {
  greta_kernel("periodic",
    tf_name = "tf_periodic",
    parameters = list(
      lengthscale = lengthscale,
      variance = variance,
      period = period
    )
  )
}
