# convergence of the EM algorithm is judged by the largest absolute change in the estimated item parameters
# (the mean and variance of the latent distributions only when all items are fixed in FIPC)

# largest absolute change between two fits in the item parameters and in the listed latent distribution
# moments; 'items' returns the item parameter columns and 'mom' the moments as a numeric vector
em_change <- function(f1, f0, items, mom) {
  d_item <- abs(as.matrix(items(f1)) - as.matrix(items(f0)))
  max(c(d_item, abs(mom(f1) - mom(f0))), na.rm = TRUE)
}

irt_items <- function(f) f$par.est[, grepl("^par\\.", names(f$par.est)), drop = FALSE]
irt_mom <- function(f) unlist(f$group.par[1, c("mu", "sigma2")])
no_mom <- function(f) numeric(0)
mg_items <- function(f) f$par.est$overall[, grepl("^par\\.", names(f$par.est$overall)), drop = FALSE]
# mean and variance of the weights of every group, also those of a rescaled reference histogram
mg_w_mom <- function(f) {
  unlist(lapply(f$weights, function(w) {
    m <- sum(w[, 1] * w[, 2])
    c(m, sum((w[, 1] - m)^2 * w[, 2]))
  }))
}

test_that("est_irt with FIPC stops on the change in the item parameters that are estimated", {
  skip_on_cran()
  set.seed(21)
  x <- shape_df(par.drm = list(a = rlnorm(10, 0, .25), b = rnorm(10), g = NA), cats = 2, model = "2PLM")
  d <- simdat(x = x, theta = rnorm(400, 0.3, 1.2), D = 1)
  for (emp in c(FALSE, TRUE)) {
    fit <- function(MaxE) {
      # a fit that is stopped before convergence warns about it
      suppressWarnings(est_irt(x = x, data = d, D = 1, fipc = TRUE, fix.loc = 1:5, EmpHist = emp,
                               MaxE = MaxE, se = FALSE, verbose = FALSE))
    }
    f <- fit(500)
    expect_lte(f$maxpar.diff, f$Etol)
    f0 <- fit(f$niter - 1)
    expect_equal(em_change(f, f0, irt_items, no_mom), f$maxpar.diff, tolerance = 1e-8)
  }
})

test_that("est_irt with FIPC and an empirical histogram converges with the default Etol and MaxE", {
  skip_on_cran()
  set.seed(21)
  x <- shape_df(par.drm = list(a = rlnorm(10, 0, .25), b = rnorm(10), g = NA), cats = 2, model = "2PLM")
  d <- simdat(x = x, theta = rnorm(400, 0.3, 1.2), D = 1)
  # no convergence warning is issued and the cycles stop before MaxE
  expect_no_warning(
    f <- est_irt(x = x, data = d, D = 1, fipc = TRUE, fix.loc = 1:5, EmpHist = TRUE, se = FALSE, verbose = FALSE)
  )
  expect_lt(f$niter, f$MaxE)
  expect_lte(f$maxpar.diff, f$Etol)
})

test_that("est_irt without FIPC stops on the change in the item parameters", {
  skip_on_cran()
  set.seed(22)
  x <- shape_df(par.drm = list(a = rlnorm(8, 0, .25), b = rnorm(8), g = NA), cats = 2, model = "2PLM")
  d <- simdat(x = x, theta = rnorm(500), D = 1)
  fit <- function(MaxE) {
    suppressWarnings(est_irt(data = d, D = 1, model = "2PLM", cats = 2, MaxE = MaxE, se = FALSE, verbose = FALSE))
  }
  f <- fit(500)
  expect_lte(f$maxpar.diff, f$Etol)
  f0 <- fit(f$niter - 1)
  expect_equal(em_change(f, f0, irt_items, no_mom), f$maxpar.diff, tolerance = 1e-8)
})

test_that("est_mg stops on the change in the item parameters", {
  skip_on_cran()
  M <- make_mg_data(seed = 11)
  gn <- c("G1", "G2")
  for (emp in c(FALSE, TRUE)) {
    fit <- function(MaxE) {
      suppressWarnings(est_mg(x = M$x[1:2], data = M$data[1:2], group.name = gn, free.group = 2, EmpHist = emp,
                              MaxE = MaxE, se = FALSE, verbose = FALSE))
    }
    f <- fit(500)
    expect_lte(f$maxpar.diff, f$Etol)
    f0 <- fit(f$niter - 1)
    expect_equal(em_change(f, f0, mg_items, no_mom), f$maxpar.diff, tolerance = 1e-8)
  }
})

test_that("est_mg with FIPC stops on the change in the free items", {
  skip_on_cran()
  M <- make_mg_data(seed = 11)
  gn <- c("G1", "G2")
  fit <- function(MaxE) {
    suppressWarnings(est_mg(x = M$x[1:2], data = M$data[1:2], group.name = gn, free.group = 2, fipc = TRUE,
                            fix.id = c("C1", "C2"), MaxE = MaxE, se = FALSE, verbose = FALSE))
  }
  f <- fit(500)
  expect_lte(f$maxpar.diff, f$Etol)
  f0 <- fit(f$niter - 1)
  expect_equal(em_change(f, f0, mg_items, no_mom), f$maxpar.diff, tolerance = 1e-8)
})

test_that("est_irt with FIPC and all items fixed stops on the change in the mean and variance", {
  skip_on_cran()
  set.seed(23)
  x <- shape_df(par.drm = list(a = rlnorm(10, 0, .25), b = rnorm(10), g = NA), cats = 2, model = "2PLM")
  d <- simdat(x = x, theta = rnorm(400, 0.3, 1.2), D = 1)
  fit <- function(MaxE) {
    suppressWarnings(est_irt(x = x, data = d, D = 1, fipc = TRUE, fix.loc = 1:10, MaxE = MaxE, se = FALSE,
                             verbose = FALSE))
  }
  f <- fit(500)
  expect_lte(f$maxpar.diff, f$Etol)
  f0 <- fit(f$niter - 1)
  expect_equal(max(abs(irt_mom(f) - irt_mom(f0))), f$maxpar.diff, tolerance = 1e-8)
})

test_that("est_mg with FIPC, all items fixed, and an empirical histogram stops on the change in the moments", {
  skip_on_cran()
  M <- make_mg_data(seed = 11)
  gn <- c("G1", "G2")
  fit <- function(MaxE) {
    suppressWarnings(est_mg(x = M$x[1:2], data = M$data[1:2], group.name = gn, free.group = 2, fipc = TRUE,
                            fix.id = M$xall$id[1:14], EmpHist = TRUE, MaxE = MaxE, se = FALSE, verbose = FALSE))
  }
  f <- fit(500)
  expect_lte(f$maxpar.diff, f$Etol)
  f0 <- fit(f$niter - 1)
  # the moments of all groups are compared, including the rescaled histogram of the reference group
  expect_equal(max(abs(mg_w_mom(f) - mg_w_mom(f0))), f$maxpar.diff, tolerance = 1e-8)
})
