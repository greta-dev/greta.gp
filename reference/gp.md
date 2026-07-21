# Define a Gaussian process

Define Gaussian processes, and project them to new coordinates.

## Usage

``` r
gp(x, kernel, inducing = NULL, n = 1, tol = 1e-04)

project(f, x_new, kernel = NULL)
```

## Arguments

- x, x_new:

  greta array giving the coordinates at which to evaluate the Gaussian
  process

- kernel:

  a kernel function created using one of the
  [kernel()](https://greta-dev.github.io/greta.gp/reference/kernels.md)
  methods

- inducing:

  an optional greta array giving the coordinates of inducing points in a
  sparse (reduced rank) Gaussian process model

- n:

  the number of independent Gaussian processes to define with the same
  kernel

- tol:

  a numerical tolerance parameter, added to the diagonal of the
  self-covariance matrix when computing the cholesky decomposition. If
  the sampler is hitting a lot of numerical errors, increasing this
  parameter could help

- f:

  a greta array created with `gp$gp` representing the values of one or
  more Gaussian processes

## Value

A greta array

## Details

`gp()` returns a greta array representing the values of the Gaussian
process (or processes) evaluated at `x`. This Gaussian process can be
made sparse (via a reduced-rank representation of the covariance) by
providing an additional set of inducing point coordinates `inducing`.
`project()` evaluates the values of an existing Gaussian process
(created with `gp()`) to new data.

## Examples

``` r
if (FALSE) { # \dontrun{
# build a kernel function on two dimensions
k1 <- rbf(lengthscales = c(0.1, 0.2), variance = 0.6)
k2 <- bias(variance = lognormal(0, 1))
K <- k1 + k2

# use this kernel in a full-rank Gaussian process
f <- gp(1:10, K)

# or in sparse Gaussian process
f_sparse <- gp(1:10, K, inducing = c(2, 5, 8))

# project the values of the GP to new coordinates
f_new <- project(f, 11:15)

# or project with a different kernel (e.g. a sub-kernel)
f_new_bias <- project(f, 11:15, k2)
} # }
```
