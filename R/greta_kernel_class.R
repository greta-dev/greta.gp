# import greta internals
op <- greta::.internals$nodes$constructors$op
as.greta_array <- greta::.internals$greta_arrays$as.greta_array

# create a greta kernel function (to create ops)
greta_kernel <- function (kernel_name,
                          gpflow_name,
                          parameters,
                          components = NULL,
                          arguments = list()) {

  kernel_name <- paste(kernel_name, "kernel")

  parameters <- lapply(parameters, as.greta_array)

  kernel <- list(name = kernel_name,
                 parameters = parameters,
                 gpflow_method = gpflow_name,
                 components = components,
                 arguments = arguments)

  # check and get the dimension of a target matrix
  get_dim <- function (x, name = 'X') {

    x_dim <- dim(x)

    if (length(x_dim) != 2) {
      stop (name, " must be a 2D greta array",
            call. = FALSE)
    }

    x_dim

  }

  # return a function here, acting on either one or two datasets
  kernel_function <- function (X, X_prime = NULL) {

    X <- as.greta_array(X)

    if (is.null(X_prime)) {

      op_data_list <- list(operation = 'self-covariance matrix',
                           X = X)
      tf_op <- tf_self_K

      dimfun <- function (elem_list) {
        X_dim <- get_dim(elem_list[[1]], 'X')
        rep(X_dim[1], 2)
      }

    } else {

      X_prime <- as.greta_array(X_prime)

      op_data_list <- list(operation = 'covariance matrix',
                           X = X,
                           X_prime = X_prime)
      tf_op <- tf_K

      dimfun <- function (elem_list) {

        X_dim <- get_dim(elem_list[[1]], 'X')
        X_prime_dim <- get_dim(elem_list[[2]], 'X_prime')

        if (X_dim[2] != X_prime_dim[2]) {
          stop ('number of columns of X and X_prime do not match',
                call. = FALSE)
        }

        c(X_dim[1], X_prime_dim[1])
      }

    }

    # kernel parameters (as greta arrays) are getting fetched here anyway so
    # just need method to fetch/assign parameters across more complex kernels

    args <- c(op_data_list,
              kernel$parameters,
              list(dimfun = dimfun,
                   operation_args = list(greta_kernel = kernel),
                   tf_operation = tf_op))

    do.call(op, args)

  }

  # give it a class and return
  class(kernel_function) <- c('greta_kernel', class(kernel_function))
  kernel_function

}

#' @export
print.greta_kernel <- function (x, ...)
  cat(environment(x)$kernel$name, "\n")

is.greta_kernel <- function (x)
  inherits(x, "greta_kernel")

# combine greta kernel function objects
combine_greta_kernel <- function (a, b,
                                  combine = c('additive', 'multiplicative')) {

  combine <- match.arg(combine)

  if (!is.greta_kernel(a) | !is.greta_kernel(b)) {
    stop ("can only combine a greta kernel with another greta kernel",
          call. = FALSE)
  }

  kernel_a <- environment(a)$kernel
  kernel_b <- environment(b)$kernel

  gpflow_name <- switch(combine,
                        additive = 'Add',
                        multiplicative = 'Prod')

  greta_kernel(kernel_name = combine,
               gpflow_name = gpflow_name,
               parameters = c(kernel_a$parameters, kernel_b$parameters),
               components = list(kernel_a, kernel_b))

}

#' @export
`+.greta_kernel` <- function (e1, e2)
  combine_greta_kernel(e1, e2, 'additive')

#' @export
`*.greta_kernel` <- function (e1, e2)
  combine_greta_kernel(e1, e2, 'multiplicative')

# recursively iterate through nested greta kernels, creating corresponding
# gpflow kernels and replacing their parameters with tensors
recurse_kernel <- function (greta_kernel, tf_parameters, counter) {

  gpflow_fun <- gpflow$kernels[[greta_kernel$gpflow_method]]

  # if it's compound, recursively call this function on the components then
  # combine them
  if (!is.null(greta_kernel$components)) {

    a <- recurse_kernel(greta_kernel$components[[1]],
                        tf_parameters,
                        counter)

    b <- recurse_kernel(greta_kernel$components[[2]],
                        tf_parameters,
                        counter)

    gpflow_kernel <- gpflow_fun(list(a, b))

  } else {

    # get gpflow version of the basis kernel
    gpflow_kernel <- do.call(gpflow_fun,
                             greta_kernel$arguments)

    # find the relevant tensors
    n_param <- length(greta_kernel$parameters)
    previous <- counter$count
    counter$count <- counter$count + n_param
    idx <- previous + seq_len(n_param)
    tf_parameters <- tf_parameters[idx]

    # put tensors in the gpflow kernel object
    parameter_names <- names(greta_kernel$parameters)
    for (i in seq_along(tf_parameters)) {
      name <- parameter_names[i]
      gpflow_kernel[[name]] <- tf_parameters[[i]]
    }

  }

  gpflow_kernel

}

# function to create gpflow kernel from a greta kernel; called when compiling
# the tf graph
compile_gpflow_kernel <- function (greta_kernel, tf_parameters) {

  # check GPflow is available and load with correct settings
  check_gpflowr()
  counter <- new.env()
  counter$count <- 0
  recurse_kernel(greta_kernel, tf_parameters, counter)

}

# create gpflow kernel and evaluate with tensors
tf_K <- function (X, X_prime, ..., greta_kernel) {

  tf_parameters <- list(...)
  gpflow_kernel <- compile_gpflow_kernel(greta_kernel, tf_parameters)
  gpflow_kernel$K(X, X_prime)

}

tf_self_K <- function (X, ..., greta_kernel) {

  tf_parameters <- list(...)
  gpflow_kernel <- compile_gpflow_kernel(greta_kernel, tf_parameters)
  gpflow_kernel$K(X)

}
