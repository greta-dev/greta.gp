# create a greta kernel function (to create ops)
greta_kernel <- function(kernel_name,
                         tf_name,
                         parameters,
                         components = NULL,
                         arguments = list()) {
  kernel_name <- glue::glue("{kernel_name} kernel")

  parameters <- lapply(parameters, as.greta_array)

  kernel <- list(
    name = kernel_name,
    parameters = parameters,
    tf_name = tf_name,
    components = components,
    arguments = arguments
  )

  # check and get the dimension of a target matrix
  get_dim <- function(x, name = "X") {
    x_dim <- dim(x)

    if (length(x_dim) != 2) {
      msg <- cli::format_error("{.var {name}} must be a 2D greta array")
      stop(
        msg,
        call. = FALSE
      )
    }

    x_dim
  }

  # return a function here, acting on either one or two datasets
  kernel_function <- function(X, X_prime = NULL) {
    X <- as.greta_array(X)

    if (is.null(X_prime)) {
      op_data_list <- list(
        operation = "self-covariance matrix",
        X = X,
        X_prime = X
      )

      X_dim <- get_dim(X, "X")
      dim <- rep(X_dim[1], 2)
    } else {
      X_prime <- as.greta_array(X_prime)

      op_data_list <- list(
        operation = "covariance matrix",
        X = X,
        X_prime = X_prime
      )

      X_dim <- get_dim(X, "X")
      X_prime_dim <- get_dim(X_prime, "X_prime")

      if (X_dim[2] != X_prime_dim[2]) {
        msg <- cli::format_error(
          c(
            "Number of columns of {.var X} and {.var X_prime} do not match",
            "i" = "{.var X} has {.val {X_dim[2]}} columns",
            "i" = "{.var X_prime} has {.val {X_prime_dim[2]}} columns"
          )
        )
        stop(
          msg,
          call. = FALSE
        )
      }

      dim <- c(X_dim[1], X_prime_dim[1])
    }

    args <- c(op_data_list,
      greta_kernel = list(kernel),
      out_dim = list(dim)
    )

    do.call("tf_K", args)
  }

  # give it a class and return
  class(kernel_function) <- c("greta_kernel", class(kernel_function))
  kernel_function
}

#' @export
print.greta_kernel <- function(x, ...) {
  cat(environment(x)$kernel_name, "\n")
}

is.greta_kernel <- function(x) {
  inherits(x, "greta_kernel")
}

# combine greta kernel function objects
combine_greta_kernel <- function(a, b,
                                 combine = c("additive", "multiplicative")) {
  combine <- match.arg(combine)

  if (!is.greta_kernel(a) | !is.greta_kernel(b)) {
    msg <- cli::format_error(
      "Can only combine a greta kernel with another greta kernel"
      )
    stop(
      msg,
      call. = FALSE
    )
  }

  kernel_a <- environment(a)
  kernel_b <- environment(b)

  tf_name <- switch(combine,
    additive = "tf_Add",
    multiplicative = "tf_Prod"
  )

  par_idx <- c(length(kernel_a$parameters), length(kernel_b$parameters))
  par_idx <- list(1:par_idx[1], par_idx[1] + 1:par_idx[2])

  greta_kernel(
    kernel_name = glue::glue("{combine} kernel"),
    tf_name = tf_name,
    parameters = c(kernel_a$parameters, kernel_b$parameters),
    components = list(kernel_a, kernel_b),
    arguments = list(parameter_idx = par_idx)
  )
}

#' @export
`+.greta_kernel` <- function(e1, e2) {
  combine_greta_kernel(e1, e2, "additive")
}

#' @export
`*.greta_kernel` <- function(e1, e2) {
  combine_greta_kernel(e1, e2, "multiplicative")
}

# recursively iterate through nested greta kernels, creating corresponding
# kernels and replacing their parameters with tensors
recurse_kernel <- function(operation, X, X_prime, greta_kernel, greta_parameters, out_dim) {

  # if it's compound, recursively call this function on the components then
  # combine them
  if (!is.null(greta_kernel$components)) {

    # parameters are in a big list for both kernels; need to
    # pull out correct pars for each kernel
    par_idx <- greta_kernel$arguments$parameter_idx

    # recursively call and calculate component kernels
    a <- do.call(recurse_kernel, list(
      operation = operation,
      X = X, X_prime = X_prime,
      greta_kernel = greta_kernel$components[[1]],
      greta_parameters = greta_kernel$parameters[par_idx[[1]]],
      out_dim = out_dim
    ))
    b <- do.call(recurse_kernel, list(
      operation = operation,
      X = X, X_prime = X_prime,
      greta_kernel = greta_kernel$components[[2]],
      greta_parameters = greta_kernel$parameters[par_idx[[2]]],
      out_dim = out_dim
    ))

    # combine evaluated base kernels
    tf_kernel <- op(operation,
      kernel_a = a, kernel_b = b,
      tf_operation = greta_kernel$tf_name,
      dim = out_dim
    )
  } else {

    # get tf version of the basis kernel
    tf_kernel <- do.call(op, c(
      list(
        operation = operation,
        X = X, X_prime = X_prime
      ),
      greta_parameters,
      list(
        operation_args = greta_kernel$arguments,
        tf_operation = greta_kernel$tf_name,
        dim = out_dim
      )
    ))
  }

  tf_kernel
}

# internal function to calculate kernel K
#  could remove this and call recurse_kernel directly
tf_K <- function(operation, X, X_prime, greta_kernel, out_dim) {
  tf_kernel <- recurse_kernel(
    operation, X, X_prime,
    greta_kernel, greta_kernel$parameters,
    out_dim
  )
}
