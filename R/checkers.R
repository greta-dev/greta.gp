check_active_dims <- function(columns, lengthscales) {
  columns <- as.integer(columns)

  if (!all(columns >= 1)) {
    msg <- cli::format_error(
      c(
        "Columns must be a vector of positive integers",
        "But we see columns as {columns}"
      )
    )
    stop(
      msg,
      call. = FALSE
    )
  }

  if (length(columns) != length(lengthscales)) {
    msg <- cli::format_error(
      c(
        "Columns has length: {.val {length(columns)}}, but the kernel has \\
         dimension: {.val {length(lengthscales)}}"
      )
    )
    stop(
      msg,
      call. = FALSE
    )
  }

  columns - 1L
}
