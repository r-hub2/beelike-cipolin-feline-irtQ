# simdat() simulates item responses from item metadata or from item parameter vectors

test_that("simdat() treats NA guessing parameters as zeros when parameters are given directly", {
  th <- rnorm(200)
  set.seed(1)
  res_na <- simdat(
    theta = th, a.drm = c(1, 1, 1), b.drm = c(0, 0, 0),
    g.drm = c(0.2, NA, NA), cats = c(2, 2, 2), D = 1
  )
  set.seed(1)
  res_zero <- simdat(
    theta = th, a.drm = c(1, 1, 1), b.drm = c(0, 0, 0),
    g.drm = c(0.2, 0, 0), cats = c(2, 2, 2), D = 1
  )
  # no response is missing and NA behaves as a zero guessing parameter
  expect_false(anyNA(res_na))
  expect_identical(res_na, res_zero)
  expect_true(all(res_na %in% c(0, 1)))
})

test_that("simdat() with NA guessing parameters matches the item metadata input", {
  th <- rnorm(200)
  x <- shape_df(
    par.drm = list(a = c(1, 1.2, 0.8), b = c(0, 0.5, -0.5), g = c(0.2, NA, NA)),
    cats = c(2, 2, 2), model = c("3PLM", "2PLM", "1PLM")
  )
  set.seed(5)
  res_meta <- simdat(x = x, theta = th, D = 1)
  set.seed(5)
  res_par <- simdat(
    theta = th, a.drm = c(1, 1.2, 0.8), b.drm = c(0, 0.5, -0.5),
    g.drm = c(0.2, NA, NA), cats = c(2, 2, 2), D = 1
  )
  expect_identical(res_par, res_meta)
})

test_that("simdat() errors when dichotomous parameter vectors do not match cats", {
  th <- c(-1, 0, 1)
  # two slopes for four dichotomous items
  expect_error(
    simdat(theta = th, a.drm = c(1, 1), b.drm = c(0, 0), cats = c(2, 2, 2, 2)),
    "one value per item"
  )
  # g vector of the wrong length (a single value is repeated)
  expect_error(
    simdat(theta = th, a.drm = c(1, 1), b.drm = c(0, 0), g.drm = c(0.2, 0.2, 0.2), cats = c(2, 2)),
    "one value per item"
  )
})

test_that("simdat() errors when polytomous parameters do not match cats", {
  th <- c(-1, 0, 1)
  # three thresholds for an item with three categories
  expect_error(
    simdat(theta = th, a.prm = 1, d.prm = list(c(-1, 0, 1)), cats = 3, pr.model = "GRM"),
    "cats - 1 thresholds"
  )
  # two polytomous items with a single slope
  expect_error(
    simdat(
      theta = th, a.prm = 1, d.prm = list(c(0, 1), c(0, 1)), cats = c(3, 3),
      pr.model = "GRM"
    ),
    "cats - 1 thresholds"
  )
})

test_that("simdat() repeats a single pr.model for every polytomous item", {
  th <- rnorm(100)
  set.seed(3)
  res1 <- simdat(
    theta = th, a.prm = c(1, 1.2), d.prm = list(c(-1, 0.5), c(0, 1)),
    cats = c(3, 3), pr.model = "GPCM"
  )
  set.seed(3)
  res2 <- simdat(
    theta = th, a.prm = c(1, 1.2), d.prm = list(c(-1, 0.5), c(0, 1)),
    cats = c(3, 3), pr.model = c("GPCM", "GPCM")
  )
  expect_identical(res1, res2)
  expect_true(all(res1 %in% 0:2))
})

test_that("simdat() reproduces the category probabilities", {
  skip_on_cran()
  meta <- data.frame(
    id = paste0("I", 1:4), cats = c(2, 4, 2, 5), model = c("3PLM", "GRM", "2PLM", "GPCM"),
    par.1 = c(1.2, 1.4, 0.8, 0.9), par.2 = c(0.3, -1, -0.5, -1), par.3 = c(0.2, 0, NA, 0),
    par.4 = c(NA, 1, NA, 0.7), par.5 = c(NA, NA, NA, 1.6)
  )
  n <- 200000
  set.seed(11)
  resp <- simdat(meta, theta = rep(0.4, n), D = 1.702)
  pr <- traceline(meta, 0.4, D = 1.702)$prob.cats

  # observed proportions are within five standard errors of the model probabilities
  for (j in 1:4) {
    obs <- tabulate(resp[, j] + 1, meta$cats[j]) / n
    expect_lt(max(abs(obs - c(pr[[j]]))), 5 * 0.5 / sqrt(n))
  }

  # the result is always a matrix, whatever the numbers of examinees and items
  expect_equal(dim(simdat(meta, theta = 0.1)), c(1L, 4L))
  expect_equal(dim(simdat(meta[2, ], theta = c(-1, 1))), c(2L, 1L))
})

test_that("simdat() repeats a single guessing value for every dichotomous item", {
  th <- rnorm(100)
  set.seed(8)
  res1 <- simdat(theta = th, a.drm = c(1, 1.2, 0.8), b.drm = c(0, 0.5, -0.5), g.drm = 0.2, cats = c(2, 2, 2))
  set.seed(8)
  res2 <- simdat(
    theta = th, a.drm = c(1, 1.2, 0.8), b.drm = c(0, 0.5, -0.5), g.drm = rep(0.2, 3),
    cats = c(2, 2, 2)
  )
  expect_identical(res1, res2)
})
