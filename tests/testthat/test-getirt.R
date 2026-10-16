# getirt() extracts every element of est_irt, est_mg, and est_item objects

test_that("every element of an est_irt object can be extracted", {
  set.seed(5)
  x <- shape_df(par.drm = list(a = rlnorm(5, 0, .3), b = rnorm(5), g = NA), cats = 2, model = "2PLM")
  d <- simdat(x = x, theta = rnorm(300), D = 1)
  fit <- est_irt(data = d, D = 1, model = "2PLM", cats = 2, verbose = FALSE, se = FALSE)
  for (e in names(fit)) {
    expect_identical(getirt(fit, what = e), fit[[e]], info = e)
  }
  expect_identical(getirt(fit, what = "call"), fit$call)
  expect_error(getirt(fit, what = "nothing"), "Could not extract")
})

test_that("every element of an est_mg object can be extracted", {
  M <- make_mg_data(seed = 11)
  # the run stops after three cycles, so it warns about convergence
  fit <- suppressWarnings(est_mg(x = M$x[1:2], data = M$data[1:2], group.name = c("G1", "G2"), free.group = 2,
                                 MaxE = 3, verbose = FALSE))
  for (e in names(fit)) {
    expect_identical(getirt(fit, what = e), fit[[e]], info = e)
  }
  expect_identical(getirt(fit, what = "call"), fit$call)
  expect_identical(names(getirt(fit, what = "npar.detail")), c("item", "dist", "weight"))
  expect_error(getirt(fit, what = "nothing"), "Could not extract")
})

test_that("every element of an est_item object can be extracted", {
  set.seed(3)
  x <- shape_df(par.drm = list(a = rlnorm(5, 0, .3), b = rnorm(5), g = NA), cats = 2, model = "2PLM")
  th <- rnorm(300)
  d <- simdat(x = x, theta = th, D = 1)
  fit <- est_item(x = x, data = d, score = th, verbose = FALSE)
  for (e in names(fit)) {
    expect_identical(getirt(fit, what = e), fit[[e]], info = e)
  }
  expect_identical(getirt(fit, what = "group.par"), c(mu = mean(th), sigma = stats::sd(th)))
  expect_identical(getirt(fit, what = "call"), fit$call)
})
