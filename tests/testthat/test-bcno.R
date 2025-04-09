# Density ------------------------------------------------------------------------------------------
test_that("dBCS works", {
  x <- 2
  mu <- stats::runif(1, 3, 5)
  sigma <- stats::runif(1, 0.2, 0.8)
  lambda <- stats::runif(1, -1, 1)
  family <- "NO"

  # Inconsistent arguments
  expect_equal(dBCS(-x, mu, sigma, lambda, family = family), 0)
  expect_warning(expect_true(is.nan(dBCS(x, -mu, sigma, lambda, family = family))))
  expect_warning(expect_true(is.nan(dBCS(x, mu, -sigma, lambda, family = family))))
  expect_error(dBCS(x, mu, sigma, lambda, family = "Other"))
  expect_error(dBCS(x, mu, sigma, lambda, family = "ST"))

  # Vectorization
  x <- matrix(stats::runif(3 * 4, 3, 5), ncol = 4)
  expect_equal(sum(dBCS(x, mu, sigma, lambda, family = family) -
    matrix(c(
      dBCS(x[1, 1], mu, sigma, lambda, family = family),
      dBCS(x[2, 1], mu, sigma, lambda, family = family),
      dBCS(x[3, 1], mu, sigma, lambda, family = family),
      dBCS(x[1, 2], mu, sigma, lambda, family = family),
      dBCS(x[2, 2], mu, sigma, lambda, family = family),
      dBCS(x[3, 2], mu, sigma, lambda, family = family),
      dBCS(x[1, 3], mu, sigma, lambda, family = family),
      dBCS(x[2, 3], mu, sigma, lambda, family = family),
      dBCS(x[3, 3], mu, sigma, lambda, family = family),
      dBCS(x[1, 4], mu, sigma, lambda, family = family),
      dBCS(x[2, 4], mu, sigma, lambda, family = family),
      dBCS(x[3, 4], mu, sigma, lambda, family = family)
    ), ncol = 4)), 0)

  mu <- matrix(stats::runif(3 * 4, 3, 5), ncol = 4)
  expect_equal(sum(dBCS(x, mu, sigma, lambda, family = family) -
    matrix(c(
      dBCS(x[1, 1], mu[1, 1], sigma, lambda, family = family),
      dBCS(x[2, 1], mu[2, 1], sigma, lambda, family = family),
      dBCS(x[3, 1], mu[3, 1], sigma, lambda, family = family),
      dBCS(x[1, 2], mu[1, 2], sigma, lambda, family = family),
      dBCS(x[2, 2], mu[2, 2], sigma, lambda, family = family),
      dBCS(x[3, 2], mu[3, 2], sigma, lambda, family = family),
      dBCS(x[1, 3], mu[1, 3], sigma, lambda, family = family),
      dBCS(x[2, 3], mu[2, 3], sigma, lambda, family = family),
      dBCS(x[3, 3], mu[3, 3], sigma, lambda, family = family),
      dBCS(x[1, 4], mu[1, 4], sigma, lambda, family = family),
      dBCS(x[2, 4], mu[2, 4], sigma, lambda, family = family),
      dBCS(x[3, 4], mu[3, 4], sigma, lambda, family = family)
    ), ncol = 4)), 0)


  sigma <- matrix(stats::runif(3 * 4, 1, 2), ncol = 4)
  expect_equal(sum(dBCS(x, mu, sigma, lambda, family = family) -
    matrix(c(
      dBCS(x[1, 1], mu[1, 1], sigma[1, 1], lambda, family = family),
      dBCS(x[2, 1], mu[2, 1], sigma[2, 1], lambda, family = family),
      dBCS(x[3, 1], mu[3, 1], sigma[3, 1], lambda, family = family),
      dBCS(x[1, 2], mu[1, 2], sigma[1, 2], lambda, family = family),
      dBCS(x[2, 2], mu[2, 2], sigma[2, 2], lambda, family = family),
      dBCS(x[3, 2], mu[3, 2], sigma[3, 2], lambda, family = family),
      dBCS(x[1, 3], mu[1, 3], sigma[1, 3], lambda, family = family),
      dBCS(x[2, 3], mu[2, 3], sigma[2, 3], lambda, family = family),
      dBCS(x[3, 3], mu[3, 3], sigma[3, 3], lambda, family = family),
      dBCS(x[1, 4], mu[1, 4], sigma[1, 4], lambda, family = family),
      dBCS(x[2, 4], mu[2, 4], sigma[2, 4], lambda, family = family),
      dBCS(x[3, 4], mu[3, 4], sigma[3, 4], lambda, family = family)
    ), ncol = 4)), 0)

  lambda <- matrix(stats::runif(3 * 4, -1, 1), ncol = 4)
  expect_equal(sum(dBCS(x, mu, sigma, lambda, family = family) -
    matrix(c(
      dBCS(x[1, 1], mu[1, 1], sigma[1, 1], lambda[1, 1], family = family),
      dBCS(x[2, 1], mu[2, 1], sigma[2, 1], lambda[2, 1], family = family),
      dBCS(x[3, 1], mu[3, 1], sigma[3, 1], lambda[3, 1], family = family),
      dBCS(x[1, 2], mu[1, 2], sigma[1, 2], lambda[1, 2], family = family),
      dBCS(x[2, 2], mu[2, 2], sigma[2, 2], lambda[2, 2], family = family),
      dBCS(x[3, 2], mu[3, 2], sigma[3, 2], lambda[3, 2], family = family),
      dBCS(x[1, 3], mu[1, 3], sigma[1, 3], lambda[1, 3], family = family),
      dBCS(x[2, 3], mu[2, 3], sigma[2, 3], lambda[2, 3], family = family),
      dBCS(x[3, 3], mu[3, 3], sigma[3, 3], lambda[3, 3], family = family),
      dBCS(x[1, 4], mu[1, 4], sigma[1, 4], lambda[1, 4], family = family),
      dBCS(x[2, 4], mu[2, 4], sigma[2, 4], lambda[2, 4], family = family),
      dBCS(x[3, 4], mu[3, 4], sigma[3, 4], lambda[3, 4], family = family)
    ), ncol = 4)), 0)


  # Inconsistent arguments with vectorization
  x[1, 1] <- -2
  mu[1, 2] <- -2
  mu[2, 3] <- -1
  sigma[2, 1] <- -2
  sigma[3, 4] <- -2

  expect_warning(expect_equal(sum(is.nan(dBCS(x, mu, sigma, lambda, family = family))), 4))
  expect_warning(expect_equal(sum(dBCS(x, mu, sigma, lambda, family = family) == 0, na.rm = TRUE), 1))
})

# CDF ----------------------------------------------------------------------------------------------
test_that("pBCS works", {
  x <- 2
  mu <- stats::runif(1, 3, 5)
  sigma <- stats::runif(1, 2, 3)
  lambda <- stats::runif(1, -1, 1)
  family <- "NO"

  # Inconsistent arguments
  expect_equal(pBCS(-x, mu, sigma, lambda, family = family), 0)
  expect_warning(expect_true(is.nan(pBCS(x, -mu, sigma, lambda, family = family))))
  expect_warning(expect_true(is.nan(pBCS(x, mu, -sigma, lambda, family = family))))
  expect_error(pBCS(x, mu, sigma, lambda, family = "Other"))
  expect_error(pBCS(x, mu, sigma, lambda, family = "ST"))

  # Vectorization
  x <- matrix(stats::runif(3 * 4, 3, 5), ncol = 4)
  expect_equal(sum(pBCS(x, mu, sigma, lambda, family = family) -
    matrix(c(
      pBCS(x[1, 1], mu, sigma, lambda, family = family),
      pBCS(x[2, 1], mu, sigma, lambda, family = family),
      pBCS(x[3, 1], mu, sigma, lambda, family = family),
      pBCS(x[1, 2], mu, sigma, lambda, family = family),
      pBCS(x[2, 2], mu, sigma, lambda, family = family),
      pBCS(x[3, 2], mu, sigma, lambda, family = family),
      pBCS(x[1, 3], mu, sigma, lambda, family = family),
      pBCS(x[2, 3], mu, sigma, lambda, family = family),
      pBCS(x[3, 3], mu, sigma, lambda, family = family),
      pBCS(x[1, 4], mu, sigma, lambda, family = family),
      pBCS(x[2, 4], mu, sigma, lambda, family = family),
      pBCS(x[3, 4], mu, sigma, lambda, family = family)
    ), ncol = 4)), 0)

  mu <- matrix(stats::runif(3 * 4, 3, 5), ncol = 4)
  expect_equal(sum(pBCS(x, mu, sigma, lambda, family = family) -
    matrix(c(
      pBCS(x[1, 1], mu[1, 1], sigma, lambda, family = family),
      pBCS(x[2, 1], mu[2, 1], sigma, lambda, family = family),
      pBCS(x[3, 1], mu[3, 1], sigma, lambda, family = family),
      pBCS(x[1, 2], mu[1, 2], sigma, lambda, family = family),
      pBCS(x[2, 2], mu[2, 2], sigma, lambda, family = family),
      pBCS(x[3, 2], mu[3, 2], sigma, lambda, family = family),
      pBCS(x[1, 3], mu[1, 3], sigma, lambda, family = family),
      pBCS(x[2, 3], mu[2, 3], sigma, lambda, family = family),
      pBCS(x[3, 3], mu[3, 3], sigma, lambda, family = family),
      pBCS(x[1, 4], mu[1, 4], sigma, lambda, family = family),
      pBCS(x[2, 4], mu[2, 4], sigma, lambda, family = family),
      pBCS(x[3, 4], mu[3, 4], sigma, lambda, family = family)
    ), ncol = 4)), 0)


  sigma <- matrix(stats::runif(3 * 4, 1, 2), ncol = 4)
  expect_equal(sum(pBCS(x, mu, sigma, lambda, family = family) -
    matrix(c(
      pBCS(x[1, 1], mu[1, 1], sigma[1, 1], lambda, family = family),
      pBCS(x[2, 1], mu[2, 1], sigma[2, 1], lambda, family = family),
      pBCS(x[3, 1], mu[3, 1], sigma[3, 1], lambda, family = family),
      pBCS(x[1, 2], mu[1, 2], sigma[1, 2], lambda, family = family),
      pBCS(x[2, 2], mu[2, 2], sigma[2, 2], lambda, family = family),
      pBCS(x[3, 2], mu[3, 2], sigma[3, 2], lambda, family = family),
      pBCS(x[1, 3], mu[1, 3], sigma[1, 3], lambda, family = family),
      pBCS(x[2, 3], mu[2, 3], sigma[2, 3], lambda, family = family),
      pBCS(x[3, 3], mu[3, 3], sigma[3, 3], lambda, family = family),
      pBCS(x[1, 4], mu[1, 4], sigma[1, 4], lambda, family = family),
      pBCS(x[2, 4], mu[2, 4], sigma[2, 4], lambda, family = family),
      pBCS(x[3, 4], mu[3, 4], sigma[3, 4], lambda, family = family)
    ), ncol = 4)), 0)

  lambda <- matrix(stats::runif(3 * 4, -1, 1), ncol = 4)
  expect_equal(sum(pBCS(x, mu, sigma, lambda, family = family) -
    matrix(c(
      pBCS(x[1, 1], mu[1, 1], sigma[1, 1], lambda[1, 1], family = family),
      pBCS(x[2, 1], mu[2, 1], sigma[2, 1], lambda[2, 1], family = family),
      pBCS(x[3, 1], mu[3, 1], sigma[3, 1], lambda[3, 1], family = family),
      pBCS(x[1, 2], mu[1, 2], sigma[1, 2], lambda[1, 2], family = family),
      pBCS(x[2, 2], mu[2, 2], sigma[2, 2], lambda[2, 2], family = family),
      pBCS(x[3, 2], mu[3, 2], sigma[3, 2], lambda[3, 2], family = family),
      pBCS(x[1, 3], mu[1, 3], sigma[1, 3], lambda[1, 3], family = family),
      pBCS(x[2, 3], mu[2, 3], sigma[2, 3], lambda[2, 3], family = family),
      pBCS(x[3, 3], mu[3, 3], sigma[3, 3], lambda[3, 3], family = family),
      pBCS(x[1, 4], mu[1, 4], sigma[1, 4], lambda[1, 4], family = family),
      pBCS(x[2, 4], mu[2, 4], sigma[2, 4], lambda[2, 4], family = family),
      pBCS(x[3, 4], mu[3, 4], sigma[3, 4], lambda[3, 4], family = family)
    ), ncol = 4)), 0)


  # Inconsistent arguments and vectorization
  x[1, 1] <- -2
  mu[1, 2] <- -2
  mu[2, 3] <- -1
  sigma[2, 1] <- -2
  sigma[3, 4] <- -2
  expect_warning(expect_equal(sum(is.nan(pBCS(x, mu, sigma, lambda, family = family))), 4))
  expect_warning(expect_equal(sum(pBCS(x, mu, sigma, lambda, family = family) == 0, na.rm = TRUE), 1))
})

# Quantile function --------------------------------------------------------------------------------
test_that("qBCS works", {
  x <- stats::runif(1, 0, 1)
  mu <- stats::runif(1, 3, 5)
  sigma <- stats::runif(1, 2, 3)
  lambda <- stats::runif(1, -1, 1)
  family <- "NO"

  # Inconsistent arguments
  expect_warning(expect_true(is.nan(qBCS(stats::runif(1, -100, 0), mu, sigma, lambda, family = family))))
  expect_warning(expect_true(is.nan(qBCS(stats::runif(1, 1, 100), mu, sigma, lambda, family = family))))
  expect_warning(expect_true(is.nan(qBCS(x, -mu, sigma, lambda, family = family))))
  expect_warning(expect_true(is.nan(qBCS(x, mu, -sigma, lambda, family = family))))

  # Vectorization
  x <- matrix(stats::runif(3 * 4, 0, 1), ncol = 4)
  expect_equal(sum(qBCS(x, mu, sigma, lambda, family = family) -
    matrix(c(
      qBCS(x[1, 1], mu, sigma, lambda, family = family),
      qBCS(x[2, 1], mu, sigma, lambda, family = family),
      qBCS(x[3, 1], mu, sigma, lambda, family = family),
      qBCS(x[1, 2], mu, sigma, lambda, family = family),
      qBCS(x[2, 2], mu, sigma, lambda, family = family),
      qBCS(x[3, 2], mu, sigma, lambda, family = family),
      qBCS(x[1, 3], mu, sigma, lambda, family = family),
      qBCS(x[2, 3], mu, sigma, lambda, family = family),
      qBCS(x[3, 3], mu, sigma, lambda, family = family),
      qBCS(x[1, 4], mu, sigma, lambda, family = family),
      qBCS(x[2, 4], mu, sigma, lambda, family = family),
      qBCS(x[3, 4], mu, sigma, lambda, family = family)
    ), ncol = 4)), 0)

  mu <- matrix(stats::runif(3 * 4, 3, 5), ncol = 4)
  expect_equal(sum(qBCS(x, mu, sigma, lambda, family = family) -
    matrix(c(
      qBCS(x[1, 1], mu[1, 1], sigma, lambda, family = family),
      qBCS(x[2, 1], mu[2, 1], sigma, lambda, family = family),
      qBCS(x[3, 1], mu[3, 1], sigma, lambda, family = family),
      qBCS(x[1, 2], mu[1, 2], sigma, lambda, family = family),
      qBCS(x[2, 2], mu[2, 2], sigma, lambda, family = family),
      qBCS(x[3, 2], mu[3, 2], sigma, lambda, family = family),
      qBCS(x[1, 3], mu[1, 3], sigma, lambda, family = family),
      qBCS(x[2, 3], mu[2, 3], sigma, lambda, family = family),
      qBCS(x[3, 3], mu[3, 3], sigma, lambda, family = family),
      qBCS(x[1, 4], mu[1, 4], sigma, lambda, family = family),
      qBCS(x[2, 4], mu[2, 4], sigma, lambda, family = family),
      qBCS(x[3, 4], mu[3, 4], sigma, lambda, family = family)
    ), ncol = 4)), 0)


  sigma <- matrix(stats::runif(3 * 4, 1, 2), ncol = 4)
  expect_equal(sum(qBCS(x, mu, sigma, lambda, family = family) -
    matrix(c(
      qBCS(x[1, 1], mu[1, 1], sigma[1, 1], lambda, family = family),
      qBCS(x[2, 1], mu[2, 1], sigma[2, 1], lambda, family = family),
      qBCS(x[3, 1], mu[3, 1], sigma[3, 1], lambda, family = family),
      qBCS(x[1, 2], mu[1, 2], sigma[1, 2], lambda, family = family),
      qBCS(x[2, 2], mu[2, 2], sigma[2, 2], lambda, family = family),
      qBCS(x[3, 2], mu[3, 2], sigma[3, 2], lambda, family = family),
      qBCS(x[1, 3], mu[1, 3], sigma[1, 3], lambda, family = family),
      qBCS(x[2, 3], mu[2, 3], sigma[2, 3], lambda, family = family),
      qBCS(x[3, 3], mu[3, 3], sigma[3, 3], lambda, family = family),
      qBCS(x[1, 4], mu[1, 4], sigma[1, 4], lambda, family = family),
      qBCS(x[2, 4], mu[2, 4], sigma[2, 4], lambda, family = family),
      qBCS(x[3, 4], mu[3, 4], sigma[3, 4], lambda, family = family)
    ), ncol = 4)), 0)

  lambda <- matrix(stats::runif(3 * 4, -1, 1), ncol = 4)
  expect_equal(sum(qBCS(x, mu, sigma, lambda, family = family) -
    matrix(c(
      qBCS(x[1, 1], mu[1, 1], sigma[1, 1], lambda[1, 1], family = family),
      qBCS(x[2, 1], mu[2, 1], sigma[2, 1], lambda[2, 1], family = family),
      qBCS(x[3, 1], mu[3, 1], sigma[3, 1], lambda[3, 1], family = family),
      qBCS(x[1, 2], mu[1, 2], sigma[1, 2], lambda[1, 2], family = family),
      qBCS(x[2, 2], mu[2, 2], sigma[2, 2], lambda[2, 2], family = family),
      qBCS(x[3, 2], mu[3, 2], sigma[3, 2], lambda[3, 2], family = family),
      qBCS(x[1, 3], mu[1, 3], sigma[1, 3], lambda[1, 3], family = family),
      qBCS(x[2, 3], mu[2, 3], sigma[2, 3], lambda[2, 3], family = family),
      qBCS(x[3, 3], mu[3, 3], sigma[3, 3], lambda[3, 3], family = family),
      qBCS(x[1, 4], mu[1, 4], sigma[1, 4], lambda[1, 4], family = family),
      qBCS(x[2, 4], mu[2, 4], sigma[2, 4], lambda[2, 4], family = family),
      qBCS(x[3, 4], mu[3, 4], sigma[3, 4], lambda[3, 4], family = family)
    ), ncol = 4)), 0)


  # Inconsistent arguments and vectorization
  x[1, 1] <- -2
  mu[1, 2] <- -2
  mu[2, 3] <- -1
  sigma[2, 1] <- -2
  sigma[3, 4] <- -2
  expect_warning(expect_equal(sum(is.nan(qBCS(x, mu, sigma, lambda, family = family))), 5))
})
