# est_item() checked against independent formulas

set.seed(7)
xi <- shape_df(par.drm = list(a = c(1, 1, 1, rlnorm(5, 0, .3)), b = rnorm(8), g = c(rep(NA, 5), .15, .2, .1)),
               par.prm = list(a = rlnorm(2, 0, .3), d = list(sort(rnorm(2)), sort(rnorm(3)))),
               cats = c(rep(2, 8), 3, 4), model = c(rep("1PLM", 3), "2PLM", "2PLM", rep("3PLM", 3), "GRM", "GPCM"))
th <- rnorm(600)
di <- simdat(x = xi, theta = th, D = 1.702)
di[sample(length(di), 0.08 * length(di))] <- NA

# independent log-likelihood of the responses given the abilities
indep_item_ll <- function(par, data, theta, D) {
  sum(sapply(seq_len(nrow(par)), function(j) {
    o <- which(!is.na(data[, j]))
    p <- indep_prob(theta[o], par$model[j], as.numeric(unlist(par[j, -(1:3)])), par$cats[j], D)
    sum(log(p[cbind(seq_along(o), data[o, j] + 1)]))
  }))
}

test_that("the log-likelihood excludes the log prior densities", {
  skip_on_cran()
  f <- est_item(x = xi, data = di, score = th, D = 1.702, verbose = FALSE)
  expect_equal(f$loglikelihood, indep_item_ll(f$par.est, di, th, 1.702), tolerance = 1e-8)
  fp <- est_item(x = xi, data = di, score = th, D = 1.702, use.aprior = TRUE, use.bprior = TRUE, verbose = FALSE)
  expect_equal(fp$loglikelihood, indep_item_ll(fp$par.est, di, th, 1.702), tolerance = 1e-8)
  # without prior distributions, the log-likelihood is the maximized objective
  fn <- est_item(x = xi, data = di, score = th, D = 1.702, use.gprior = FALSE, verbose = FALSE)
  expect_equal(fn$loglikelihood, indep_item_ll(fn$par.est, di, th, 1.702), tolerance = 1e-8)
})

test_that("the standard errors match the inverse of an independent observed information with the g prior", {
  skip_on_cran()
  f <- est_item(x = xi, data = di, score = th, D = 1.702, verbose = FALSE)
  pos <- as.matrix(f$pos.par[, -(1:3)])
  free <- which(!is.na(pos), arr.ind = TRUE)
  pv <- as.matrix(f$par.est[, -(1:3)])
  lpost <- function(v) {
    pm <- pv
    pm[free] <- v
    # the three 1PLM items share one slope
    pm[1:3, 1] <- pm[1, 1]
    s <- 0
    for (j in seq_len(nrow(xi))) {
      o <- which(!is.na(di[, j]))
      p <- indep_prob(th[o], xi$model[j], pm[j, !is.na(pm[j, ])], xi$cats[j], 1.702)
      s <- s + sum(log(p[cbind(seq_along(o), di[o, j] + 1)]))
    }
    s + sum(dbeta(pm[6:8, 3], 5, 16, log = TRUE))
  }
  H <- stats::optimHess(pv[free], lpost)
  cov_ind <- solve(-H)
  ord <- pos[free]
  expect_equal(as.matrix(f$se.est[, -(1:3)])[free], sqrt(diag(cov_ind)), tolerance = 1e-3)
  expect_equal(unname(f$covariance[ord, ord]), unname(cov_ind), tolerance = 1e-3)
})
