# gen.weight() returns nodes and normalized weights for a distribution

test_that("gen.weight() returns normalized weights with the requested moments", {
  w <- gen.weight(n = 41, dist = "norm", mu = 0.5, sigma = 1.3)
  expect_equal(sum(w$weight), 1)
  expect_equal(sum(w$theta * w$weight), 0.5, tolerance = 1e-10)
  expect_equal(sum(w$theta^2 * w$weight) - 0.5^2, 1.3^2, tolerance = 1e-10)
  expect_equal(gen.weight(n = 5, dist = "unif", l = -2, u = 2)$weight, rep(0.2, 5))
  expect_equal(gen.weight(dist = "emp", theta = c(1, 2, 2))$weight, rep(1 / 3, 3))
})

test_that("gen.weight() normalizes the normal density over a given theta grid", {
  theta <- seq(-4, 4, 0.5)
  w <- gen.weight(dist = "norm", mu = 0, sigma = 1, theta = theta)
  expect_equal(sum(w$weight), 1)
  expect_equal(w$weight, dnorm(theta) / sum(dnorm(theta)))
})

test_that("gen.weight() matches dist without regard to case", {
  expect_identical(gen.weight(n = 5, dist = "NORM"), gen.weight(n = 5, dist = "norm"))
})

test_that("gen.weight() errors on an unsupported dist", {
  expect_error(gen.weight(dist = "normal"), "must be one of")
  expect_error(gen.weight(n = 5, dist = c("norm", "unif")), "must be one of")
})

test_that("gen.weight() errors when theta is given with dist = 'unif'", {
  expect_error(gen.weight(dist = "unif", theta = c(-1, 0, 1)), "cannot be used")
})

test_that("gen.weight() errors when dist = 'emp' has no theta", {
  expect_error(gen.weight(dist = "emp"), "required")
})
