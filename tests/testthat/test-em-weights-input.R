# User-supplied quadrature weights of est_irt() and est_mg()

xf <- bring.flexmirt(file = system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ"), "par")$Group1$full_df
x12 <- xf[c(1:10, 39:40), ]
set.seed(1)
d12 <- simdat(x = x12, theta = rnorm(500), D = 1)

sum_msg <- "The weights do not sum to 1 and were rescaled to sum to 1."

test_that("weights that do not sum to one are rescaled with a warning and give the normalized results", {
  skip_on_cran()
  w1 <- gen.weight(dist = "norm", theta = seq(-6, 6, length.out = 49))
  w3 <- w1
  w3$weight <- 3 * w3$weight
  for (eh in c(FALSE, TRUE)) {
    r1 <- with_warnings(est_irt(x = x12, data = d12, D = 1, weights = w1, EmpHist = eh, MaxE = 20, se = FALSE, verbose = FALSE))
    r3 <- with_warnings(est_irt(x = x12, data = d12, D = 1, weights = w3, EmpHist = eh, MaxE = 20, se = FALSE, verbose = FALSE))
    f1 <- r1$value
    f3 <- r3$value
    expect_false(sum_msg %in% r1$warnings)
    expect_true(sum_msg %in% r3$warnings)
    expect_equal(f3$loglikelihood, f1$loglikelihood, tolerance = 1e-10)
    expect_equal(f3$group.par, f1$group.par, tolerance = 1e-10)
    expect_equal(f3$par.est, f1$par.est, tolerance = 1e-8)
  }
  # the reported log-likelihood is computed with weights that sum to one
  ref <- indep_marg(f1$par.est, d12, f1$weights, D = 1)
  expect_equal(f3$loglikelihood, ref$loglik, tolerance = 1e-10)
})

test_that("the rescaling of the weights also applies to FIPC and to est_mg()", {
  skip_on_cran()
  w1 <- gen.weight(dist = "norm", theta = seq(-6, 6, length.out = 49))
  w2 <- w1
  w2$weight <- 2 * w2$weight
  r2 <- with_warnings(est_irt(x = x12, data = d12, D = 1, fipc = TRUE, fix.loc = 1:6, weights = w2, MaxE = 5,
                              se = FALSE, verbose = FALSE))
  expect_true(sum_msg %in% r2$warnings)
  r1 <- with_warnings(est_irt(x = x12, data = d12, D = 1, fipc = TRUE, fix.loc = 1:6, weights = w1, MaxE = 5,
                              se = FALSE, verbose = FALSE))
  expect_equal(r2$value$loglikelihood, r1$value$loglikelihood, tolerance = 1e-10)

  set.seed(2)
  dat <- list(g1 = simdat(x = x12, theta = rnorm(300), D = 1), g2 = simdat(x = x12, theta = rnorm(300, 0.5), D = 1))
  args <- list(x = list(x12, x12), data = dat, group.name = c("G1", "G2"), free.group = "G2", D = 1, MaxE = 5,
               se = FALSE, verbose = FALSE)
  m1 <- with_warnings(do.call(est_mg, c(args, list(weights = w1))))
  m2 <- with_warnings(do.call(est_mg, c(args, list(weights = w2))))
  expect_false(sum_msg %in% m1$warnings)
  expect_true(sum_msg %in% m2$warnings)
  expect_equal(m2$value$loglikelihood, m1$value$loglikelihood, tolerance = 1e-10)
})

test_that("weights that are negative, not finite, or sum to zero stop the estimation", {
  w <- gen.weight(dist = "norm", theta = seq(-6, 6, length.out = 49))
  msg <- "nonnegative finite values with a positive sum"
  wneg <- w
  wneg$weight[1] <- -0.1
  expect_error(est_irt(x = x12, data = d12, D = 1, weights = wneg, verbose = FALSE), msg)
  wna <- w
  wna$weight[2] <- NA
  expect_error(est_irt(x = x12, data = d12, D = 1, weights = wna, verbose = FALSE), msg)
  wzero <- w
  wzero$weight <- 0
  expect_error(est_irt(x = x12, data = d12, D = 1, weights = wzero, verbose = FALSE), msg)
})

test_that("the empirical histogram works with weights whose length differs from Quadrature[1]", {
  for (nq in c(41, 61)) {
    w <- gen.weight(dist = "norm", theta = seq(-6, 6, length.out = nq))
    f <- suppressWarnings(est_irt(x = x12, data = d12, D = 1, weights = w, EmpHist = TRUE, MaxE = 20, se = FALSE,
                                  verbose = FALSE))
    ww <- f$weights[, 2]
    q <- f$weights[, 1]
    expect_equal(sum(ww), 1, tolerance = 1e-12)
    # the histogram is rescaled to mean 0 and variance 1 (up to the interpolation)
    expect_lt(abs(sum(q * ww)), 0.01)
    expect_lt(abs(sum(q^2 * ww) - sum(q * ww)^2 - 1), 0.01)
    # the extrapolated tails decrease towards the ends of the grid
    expect_true(all(diff(utils::tail(ww, 4)) < 0))
    expect_true(all(diff(utils::head(ww, 4)) > 0))
  }
})
