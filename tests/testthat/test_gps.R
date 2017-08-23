context('gaussian processes')

test_that('gaussian processes can be defined in models and sampled from', {

  source("helpers.R")
  skip_if_not(greta:::check_tf_version())
  skip_if_not(gpflowr::gpflow_available())

  k <- rbf(1, 1)

  # full
  f = gp(1:10, k)
  expect_ok(m <- model(f))
  expect_ok(mcmc(m, warmup = 2, n_samples = 2))

  # sparse
  f = gp(1:10, k,
         inducing = c(2, 4, 6, 8))
  expect_ok(m <- model(f))
  expect_ok(mcmc(m, warmup = 2, n_samples = 2))

})

