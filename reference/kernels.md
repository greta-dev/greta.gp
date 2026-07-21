# Gaussian process kernels

Create and combine Gaussian process kernels (covariance functions) for
use in Gaussian process models.

## Usage

``` r
bias(variance)

constant(variance)

white(variance)

iid(variance, columns = 1)

rbf(lengthscales, variance, columns = seq_along(lengthscales))

rational_quadratic(
  lengthscales,
  variance,
  alpha,
  columns = seq_along(lengthscales)
)

linear(variances, columns = seq_along(variances))

polynomial(variances, offset, degree, columns = seq_along(variances))

expo(lengthscales, variance, columns = seq_along(lengthscales))

mat12(lengthscales, variance, columns = seq_along(lengthscales))

mat32(lengthscales, variance, columns = seq_along(lengthscales))

mat52(lengthscales, variance, columns = seq_along(lengthscales))

cosine(lengthscales, variance, columns = seq_along(lengthscales))

periodic(period, lengthscale, variance)
```

## Arguments

- variance, variances:

  (scalar/vector) the variance of a Gaussian process prior in all
  dimensions (`variance`) or in each dimension (`variances`)

- columns:

  (scalar/vector integer, not a greta array) the columns of the data
  matrix on which this kernel acts. Must have the same dimensions as
  lengthscale parameters.

- alpha:

  (scalar) additional parameter in rational quadratic kernel

- offset:

  (scalar) offset in polynomial kernel

- degree:

  (scalar) degree of polynomial kernel

- period:

  (scalar) the period of the Gaussian process

- lengthscale, lengthscales:

  (scalar/vector) the correlation decay distance along all dimensions
  (`lengthscale`) or each dimension ((`lengthscales`)) of the Gaussian
  process

## Value

greta kernel with class "greta_kernel"

## Details

The kernel constructor functions each return a *function* (of class
`greta_kernel`) which can be executed on greta arrays to compute the
covariance matrix between points in the space of the Gaussian process.
The `+` and `*` operators can be used to combine kernel functions to
create new kernel functions.

Note that `bias` and `constant` are identical names for the same
underlying kernel.

`iid` is equivalent to `bias` where all entries in `columns` match
(where the absolute euclidean distance is less than 1e-12), and `white`
where they don't; i.e. an independent Gaussian random effect.

## Examples

``` r
if (FALSE) { # \dontrun{
# create a radial basis function kernel on two dimensions
k1 <- rbf(lengthscales = c(0.1, 0.2), variance = 0.6)

# evaluate it on a greta array to get the variance-covariance matrix
x <- greta_array(rnorm(8), dim = c(4, 2))
k1(x)

# non-symmetric covariance between two sets of points
x2 <- greta_array(rnorm(10), dim = c(5, 2))
k1(x, x2)

# create a bias kernel, with the variance as a variable
k2 <- bias(variance = lognormal(0, 1))

# combine two kernels and evaluate
K <- k1 + k2
K(x, x2)

# other kernels
constant(variance = lognormal(0, 1))
white(variance = lognormal(0, 1))
iid(variance = lognormal(0, 1))
rational_quadratic(lengthscales = c(0.1, 0.2), alpha = 0.5, variance = 0.6)
linear(variances = 0.1)
polynomial(variances = 0.6, offset = 0.8, degree = 2)
expo(lengthscales = 0.6, variance = 0.9)
mat12(lengthscales = 0.5, variance = 0.7)
mat32(lengthscales = 0.4, variance = 0.8)
mat52(lengthscales = 0.3, variance = 0.9)
cosine(lengthscales = 0.68, variance = 0.8)
periodic(period = 0.71, lengthscale = 0.59, variance = 0.2)
} # }
```
