#' @title Gaussian process modelling in greta
#' @name greta.gp
#'
#' @description A greta module to create and combine covariance functions and
#'   use them to build Gaussian process models in greta. See
#'   \code{\link{kernels}} and \code{\link{gp}}
#'
#' @docType package
#' @importFrom gpflowr gpflow
#' @importFrom greta .internals
#'
NULL


#' @importFrom reticulate py_set_attr
#' @importFrom utils capture.output
check_gpflowr <- function () {

  # create/modify a gpflow config file whilst loading gpflow
  config <- './.gpflowrc'
  extant <- file.exists(config)

  cleanup <- function () {
    file.remove(config)
    if (extant)
      file.rename(config_temp, config)
  }
  on.exit(cleanup())

  if (extant)
    file.rename(config, (config_temp <- tempfile()))

  file.create(config)

  # get float type
  type <- capture.output(options()$greta_tf_float)
  float_text <- switch(type,
                       "<dtype: 'float32'>" = "float32",
                       "<dtype: 'float64'>" = "float64")

  text <- sprintf("[dtypes]\nfloat_type = %s\nint_type = int32\n", float_text)
  cat(text, file = config)


  gpflowr_available <- requireNamespace('gpflowr', quietly = TRUE)

  if (gpflowr_available)
    gpflow_available <- gpflowr::gpflow_available()

  if (!gpflowr_available || !gpflow_available) {

    stop ('the GPflow python package and the gpflowr R package must be installed to plot greta models',
          call. = FALSE)

  } else {

    assign('gpflow',
           gpflowr::gpflow,
           envir = globalenv())

  }

}
