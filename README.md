greta.gp
========

### Gaussian processes in greta

`greta.gp` extends greta to let you define Gaussian processes as part of
your model. It provides a syntax to create and combine GP kernels, and
use them to define either full rank or sparse Gaussian processes.

[![build
status](https://travis-ci.org/greta-dev/greta.gp.svg?branch=master)](https://travis-ci.org/greta-dev/greta.gp)
[![codecov.io](https://codecov.io/github/greta-dev/greta.gp/coverage.svg?branch=master)](https://codecov.io/github/greta-dev/greta.gp?branch=master)

#### Example

``` r
# simulate data
x <- runif(20, 0, 10)
y <- sin(x) + rnorm(20, 0, 0.5)
x_plot <- seq(-1, 11, length.out = 200)
```

``` r
library(greta)
library(greta.gp)

# hyperparameters
rbf_var <- lognormal(0, 1)
rbf_len <- lognormal(0, 1)
obs_sd <- lognormal(0, 1)

# kernel & GP
kernel <- rbf(rbf_len, rbf_var) + bias(1)
f <- gp(x, kernel)

# likelihood
distribution(y) <- normal(f, obs_sd)

# prediction
f_plot <- project(f, x_plot)
```

``` r
# fit the model by Hamiltonian Monte Carlo
m <- model(f_plot)
draws <- mcmc(m)
```

``` r
# plot 200 posterior samples
plot(y ~ x, pch = 16, col = grey(0.4), xlim = c(0, 10), ylim = c(-2.5, 2.5))
for (i in 1:200) {
  lines(draws[[1]][i, ] ~ x_plot,
        lwd = 2,
        col = rgb(0.7, 0.1, 0.4, 0.1))  
}
```

![](README_files/figure-markdown_github/plotting-1.png)

#### Installation

This package depends on version 0.3.0 or higher of
[greta](https://github.com/greta-dev/greta).

You can install the R packages from within R using
[devtools](https://CRAN.R-project.org/package=devtools):

``` r
library (devtools)
install_github("greta-dev/greta")
install_github("greta-dev/greta.gp")
```
