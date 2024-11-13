# create a zero-mean Gaussian process with control points at x, and the specified kernel

#' @title Define a Gaussian process
#' @name gp
#'
#' @description Define Gaussian processes, and project them to new coordinates.
#'
#' @param x,x_new greta array giving the coordinates at which to evaluate the
#'   Gaussian process
#' @param kernel a kernel function created using one of the
#'   [kernel()][greta.gp::kernels] methods
#' @param inducing an optional greta array giving the coordinates of inducing
#'   points in a sparse (reduced rank) Gaussian process model
#' @param n the number of independent Gaussian processes to define with
#'   the same kernel
#' @param tol a numerical tolerance parameter, added to the diagonal of the
#'   self-covariance matrix when computing the cholesky decomposition. If the
#'   sampler is hitting a lot of numerical errors, increasing this parameter
#'   could help
#' @param f a greta array created with `gp$gp` representing the values of
#'   one or more Gaussian processes
#'
#' @details `gp()` returns a greta array representing the values of the
#'   Gaussian process(es) evaluated at `x`. This Gaussian process can be
#'   made sparse (via a reduced-rank representation of the covariance) by
#'   providing an additional set of inducing point coordinates `inducing`.
#'   `project()` evaluates the values of an existing Gaussian process
#'   (created with `gp()`) to new data.
#'
#' @return A greta array
#'
#' @examples
#' \dontrun{
#' # build a kernel function on two dimensions
#' k1 <- rbf(lengthscales = c(0.1, 0.2), variance = 0.6)
#' k2 <- bias(variance = lognormal(0, 1))
#' K <- k1 + k2
#'
#' # use this kernel in a full-rank Gaussian process
#' f <- gp(1:10, K)
#'
#' # or in sparse Gaussian process
#' f_sparse <- gp(1:10, K, inducing = c(2, 5, 8))
#'
#' # project the values of the GP to new coordinates
#' f_new <- project(f, 11:15)
#'
#' # or project with a different kernel (e.g. a sub-kernel)
#' f_new_bias <- project(f, 11:15, k2)
#' }
NULL

#' @rdname gp
#' @importFrom greta normal %*% forwardsolve
#' @export
gp <- function(x, kernel, inducing = NULL, n = 1, tol = 1e-4) {
  sparse <- !is.null(inducing)

  x <- as.greta_array(x)

  if (!sparse) {
    inducing <- x
  } else {
    inducing <- as.greta_array(inducing)
  }

  # calculate key objects
  m <- nrow(inducing)
  v <- normal(0, 1, dim = c(m, n))
  Kmm <- kernel(inducing)

  if (!identical(tol, 0)) {
    Kmm <- Kmm + diag(m) * tol
  }

  Lm <- t(chol(Kmm))

  # evaluate gp at x
  if (sparse) {
    Kmn <- kernel(inducing, x)
    A <- forwardsolve(Lm, Kmn)
    f <- t(A) %*% v
  } else {
    f <- Lm %*% v
  }

  # add the info to the greta array
  attr(f, "gp_info") <- list(
    kernel = kernel,
    inducing = inducing,
    v = v,
    Lm = Lm
  )
  f
}


#' @rdname gp
#' @export
project <- function(f, x_new, kernel = NULL) {
  # get the gp information and project to x_new
  info <- attr(f, "gp_info")

  if (is.null(info)) {
    cli::cli_abort(
      "Can only project from greta arrays created with {.code greta.gp::gp}"
    )
  }

  # this should be the same as the is.null pattern below but it doesn't work?!
  # kernel <- kernel %||% info$kernal

  if (is.null(kernel)){
    kernel <- info$kernel
  }

  Kmn <- kernel(info$inducing, x_new)
  A <- forwardsolve(info$Lm, Kmn)
  t(A) %*% info$v
}
