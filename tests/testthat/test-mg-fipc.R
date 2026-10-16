# fixed item parameter calibration in est_mg() and est_irt(): input checks and edge cases

M <- make_mg_data(seed = 11)
gn <- c("G1", "G2", "G3")
# fits stop after two cycles, so calls that do not test the convergence warning suppress it call by call
q <- function(...) est_mg(..., verbose = FALSE, MaxE = 2, se = FALSE)

test_that("item.id replaces the IDs in x in MG-FIPC", {
  skip_on_cran()
  ids2 <- lapply(M$x, function(z) paste0(z$id, "_n"))
  ff <- suppressWarnings(q(x = M$x, data = M$data, group.name = gn, free.group = 1:3, item.id = ids2,
                           fipc = TRUE, fix.loc = list(1:4, 1:4, NULL)))
  expect_identical(ff$par.est$group$G3$id, ids2[[3]])
  # fix.id refers to the new IDs
  fn <- suppressWarnings(q(x = M$x, data = M$data, group.name = gn, free.group = 1:3, item.id = ids2,
                           fipc = TRUE, fix.id = c("C1_n", "C2_n")))
  expect_identical(sort(fn$par.est$overall$id[fn$fix.loc$overall]), c("C1_n", "C2_n"))
  # the same IDs given in item.id reproduce the run without item.id
  id <- lapply(M$x, function(z) z$id)
  f0 <- suppressWarnings(q(x = M$x, data = M$data, group.name = gn, free.group = 1:3, fipc = TRUE, fix.id = c("C1", "C2")))
  f1 <- suppressWarnings(q(x = M$x, data = M$data, group.name = gn, free.group = 1:3, fipc = TRUE, fix.id = c("C1", "C2"),
                           item.id = id))
  expect_identical(f1$par.est, f0$par.est)
})

test_that("FIPC stops when the responses to all fixed items show no variation", {
  skip_on_cran()
  set.seed(1)
  x <- shape_df(par.drm = list(a = rlnorm(6, 0, .3), b = rnorm(6), g = NA), cats = 2, model = "2PLM")
  d <- simdat(x = x, theta = rnorm(400), D = 1)
  d[, 2] <- 1
  d[, 4] <- 0
  d[sample(400, 30), 4] <- NA
  for (m in c("MEM", "OEM")) {
    expect_error(est_irt(x = x, data = d, fipc = TRUE, fipc.method = m, fix.loc = c(2, 4), verbose = FALSE),
                 "responses to all fixed items \\(V2, V4\\) show no variation")
  }
  expect_error(est_irt(x = x[c(2, 4), ], data = d[, c(2, 4)], fipc = TRUE, fix.loc = 1:2, verbose = FALSE), "show no variation")
  expect_error(est_irt(x = x, data = d, fipc = TRUE, fipc.method = "oem", fix.loc = 6, verbose = FALSE), "'fipc.method' must be")
  # one fixed item with varying responses is enough
  f <- suppressWarnings(est_irt(x = x, data = d, fipc = TRUE, fix.loc = c(2, 4, 6), MaxE = 3, se = FALSE, verbose = FALSE))
  expect_s3_class(f, "est_irt")
  # the responses to the fixed items are pooled over the groups in est_mg()
  dd <- M$data
  dd[[1]][, 1] <- 1
  dd[[2]][, 1] <- 1
  dd[[1]][, 2] <- 0
  dd[[2]][, 2] <- 0
  expect_error(q(x = M$x, data = dd, group.name = gn, free.group = 1:3, fipc = TRUE, fix.id = c("C1", "C2")), "show no variation")
  dd[[2]][, 1] <- rep(0:1, length.out = nrow(dd[[2]]))
  f2 <- suppressWarnings(q(x = M$x, data = dd, group.name = gn, free.group = 1:3, fipc = TRUE, fix.id = c("C1", "C2")))
  expect_s3_class(f2, "est_mg")
})

test_that("a fixed item with different parameters across groups gives a warning", {
  skip_on_cran()
  xf <- M$x
  xf[[2]]$par.1[1] <- xf[[2]]$par.1[1] + 0.5
  r <- with_warnings(q(x = xf, data = M$data, group.name = gn, free.group = 1:3, fipc = TRUE, fix.loc = list(1:4, 1:4, NULL)))
  expect_true(any(grepl("fixed item\\(s\\) C1 have different parameter values", r$warnings)))
  expect_equal(r$value$par.est$overall$par.1[1], M$x[[1]]$par.1[1])
  # no warning when the values agree
  r0 <- with_warnings(q(x = M$x, data = M$data, group.name = gn, free.group = 1:3, fipc = TRUE, fix.loc = list(1:4, 1:4, NULL)))
  expect_false(any(grepl("different parameter values", r0$warnings)))
  # a difference at the level of floating point noise gives no warning
  xn <- M$x
  xn[[2]]$par.1[1] <- xn[[2]]$par.1[1] * (1 + 1e-12)
  rn <- with_warnings(q(x = xn, data = M$data, group.name = gn, free.group = 1:3, fipc = TRUE, fix.loc = list(1:4, 1:4, NULL)))
  expect_false(any(grepl("different parameter values", rn$warnings)))
  # a difference of 1e-4 gives the warning
  xm <- M$x
  xm[[2]]$par.1[1] <- xm[[2]]$par.1[1] + 1e-4
  rw <- with_warnings(q(x = xm, data = M$data, group.name = gn, free.group = 1:3, fipc = TRUE, fix.loc = list(1:4, 1:4, NULL)))
  expect_true(any(grepl("fixed item\\(s\\) C1 have different parameter values", rw$warnings)))
})

test_that("MG-FIPC with all items fixed reports that the information matrix is not computed", {
  f <- suppressWarnings(q(x = M$x, data = M$data, group.name = gn, free.group = 1:3, fipc = TRUE, fix.id = M$xall$id))
  expect_identical(f$test.2, "No item parameters are estimated, so the information matrix is not computed.")
  expect_identical(unname(f$npar.detail), c(0L, 6L, 0L))
})

test_that("an item keeps its previous values and is marked as not converged when its last search fails", {
  quadpt <- seq(-4, 4, length.out = 21)
  nan_freq <- rep(NaN, 21)
  # 2PLM item; the NaN frequencies make every search evaluate the objective at NaN
  r <- with_warnings(irtQ:::estimation2(
    f_i = nan_freq, r_i = nan_freq, s_i = nan_freq, quadpt = quadpt, mod = "2PLM", D = 1, n.quad = 21,
    control = list(eval.max = 50, iter.max = 20), startval = c(1.2, 0.3), lower = c(0.01, -6), upper = c(5, 6), iter = 5
  ))
  est <- r$value
  expect_true(length(r$warnings) > 0L)
  expect_true(all(grepl("NA/NaN function evaluation", r$warnings)))
  expect_equal(est$pars, c(1.2, 0.3))
  expect_equal(est$convergence, 1L)
  # 1PLM items that share one slope
  r1 <- with_warnings(irtQ:::estimation2(
    f_i = nan_freq, r_i = nan_freq, s_i = nan_freq, quadpt = quadpt, mod = "1PLM", D = 1, n.quad = 21,
    fix.a.1pl = FALSE, n.1PLM = 2, control = list(eval.max = 50, iter.max = 20), startval = c(1, 0.1, -0.2),
    lower = c(0.01, -6, -6), upper = c(5, 6, 6), iter = 5
  ))
  est1 <- r1$value
  expect_true(all(grepl("NA/NaN function evaluation", r1$warnings)))
  expect_equal(est1$pars, c(1, 0.1, -0.2))
  expect_equal(est1$convergence, 1L)
})

test_that("a failed search in the first EM cycles keeps the previous values", {
  quadpt <- seq(-4, 4, length.out = 21)
  f_i <- rep(20, 21)
  r_i <- round(f_i * plogis(quadpt))
  # a gradient that cannot be evaluated makes every gradient-based search fail
  fns <- irtQ:::make_drm_optim_fns
  local_mocked_bindings(make_drm_optim_fns = function(...) {
    out <- fns(...)
    out$gradient <- function(p) rep(NaN, length(p))
    out
  })
  for (it in 1:2) {
    est <- irtQ:::estimation2(
      f_i = f_i, r_i = r_i, s_i = f_i - r_i, quadpt = quadpt, mod = "2PLM", D = 1, n.quad = 21,
      control = list(eval.max = 50, iter.max = 20), startval = c(1.2, 0.3), lower = c(0.01, -6), upper = c(5, 6),
      iter = it
    )
    expect_equal(est$pars, c(1.2, 0.3))
    expect_equal(est$convergence, 1L)
  }
})
