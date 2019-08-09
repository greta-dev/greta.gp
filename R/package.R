#' @title Gaussian process modelling in greta
#' @name greta.gp
#'
#' @description A greta module to create and combine covariance functions and
#'   use them to build Gaussian process models in greta. See
#'   \code{\link{kernels}} and \code{\link{gp}}
#'
#' @docType package
#'
#' @importFrom greta .internals
#' @importFrom tensorflow tf
#'
NULL

# need some internal greta functions accessible
as.greta_array <- greta::.internals$greta_array$as.greta_array
op <- greta::.internals$nodes$constructors$op
fl <- greta::.internals$utils$misc$fl
tf_float <- greta::.internals$utils$misc$tf_float
