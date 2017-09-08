#' @title Gaussian process kernels
#' @name kernels
#'
#' @description Create and combine Gaussian process kernels (covariance
#'   functions) for use in Gaussian process models.
#'
#' @param variance,variances (scalar/vector) the variance of a Gaussian process
#'   prior in all dimensions (\code{variance}) or in each dimension
#'   (\code{variances})
#' @param lengthscale,lengthscales (scalar/vector) the correlation decay
#'   distance along all dimensions (\code{lengthscale}) or each dimension
#'   ((\code{lengthscales})) of the Gaussian process
#' @param period (scalar) the period of the Gaussian process
#' @param dim (scalar integer, not a greta array) the dimension of the Gaussian
#'   process (number of columns on which it acts)
#'
#' @details The kernel constructor functions each return a \emph{function} (of
#'   class \code{greta_kernel}) which can be executed on greta arrays
#'   to compute the covariance matrix between points in the space of the
#'   Gaussian process. The \code{+} and \code{*} operators can be used to
#'   combine kernel functions to create new kernel functions.
#'
#'   The kernels are imported from the GPflow python package, using the gpflowr
#'   R package. Both of those need to be installed before you can use these
#'   methods. See the \href{gpflow.readthedocs.io}{GPflow website} for details
#'   of the kernels implemented.
#'
#' @examples
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
NULL

#' @rdname kernels
#' @export
bias <- function (variance, dim = 1) {
  greta_kernel("bias",
               gpflow_name = "Bias",
               parameters = list(variance = variance),
               dim = dim)
}

#' @rdname kernels
#' @export
white <- function (variance, dim = 1) {
  greta_kernel("white",
               gpflow_name = "White",
               parameters = list(variance = variance),
               dim = dim)
}

#' @rdname kernels
#' @export
linear <- function (variances) {
  greta_kernel("linear",
               gpflow_name = 'Linear',
               parameters = list(variance = variances),
               dim = length(variances),
               arguments = list(ARD = TRUE))
}

#' @rdname kernels
#' @export
rbf <- function (lengthscales, variance) {
  greta_kernel("radial basis",
               gpflow_name = 'RBF',
               parameters = list(lengthscales = t(lengthscales),
                                 variance = variance),
               dim = length(lengthscales),
               arguments = list(ARD = TRUE))
}

#' @rdname kernels
#' @export
expo <- function (lengthscales, variance) {
  greta_kernel("exponential",
               gpflow_name = 'Exponential',
               parameters = list(lengthscales = t(lengthscales),
                                 variance = variance),
               dim = length(lengthscales),
               arguments = list(ARD = TRUE))
}

#' @rdname kernels
#' @export
mat12 <- function (lengthscales, variance) {
  greta_kernel("Matern 1/2",
               gpflow_name = 'Matern12',
               parameters = list(lengthscales = t(lengthscales),
                                 variance = variance),
               dim = length(lengthscales),
               arguments = list(ARD = TRUE))
}

#' @rdname kernels
#' @export
mat32 <- function (lengthscales, variance) {
  greta_kernel("Matern 3/2",
               gpflow_name = 'Matern32',
               parameters = list(lengthscales = t(lengthscales),
                                 variance = variance),
               dim = length(lengthscales),
               arguments = list(ARD = TRUE))
}

#' @rdname kernels
#' @export
mat52 <- function (lengthscales, variance) {
  greta_kernel("Matern 5/2",
               gpflow_name = 'Matern52',
               parameters = list(lengthscales = t(lengthscales),
                                 variance = variance),
               dim = length(lengthscales),
               arguments = list(ARD = TRUE))
}

#' @rdname kernels
#' @export
periodic <- function (period, lengthscale, variance, dim = 1) {
  greta_kernel("periodic",
               gpflow_name = 'PeriodicKernel',
               parameters = list(period = period,
                                 lengthscales = lengthscale,
                                 variance = variance),
               dim = dim)
}
