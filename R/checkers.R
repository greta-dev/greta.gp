check_active_dims <- function(columns, lengthscales) {
  columns <- as.integer(columns)

  columns_are_positive_int <- all(columns >= 1)
  if (!columns_are_positive_int) {
    cli::cli_abort(
      c(
        "Columns must be a vector of positive integers",
        "But we see columns as {columns}"
      )
    )
  }

  column_length_match_kernel_length <- length(columns) != length(lengthscales)
  if (column_length_match_kernel_length) {
    cli::cli_abort(
      c(
        "Columns has length: {.val {length(columns)}}, but the kernel has \\
         dimension: {.val {length(lengthscales)}}"
      )
    )
  }

  columns - 1L
}
