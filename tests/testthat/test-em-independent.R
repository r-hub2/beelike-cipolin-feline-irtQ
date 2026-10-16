# Independent checks of the EM results: marginal log-likelihood, posterior,
# expected counts, AIC/BIC, and the FIPC (MEM) fixed point

set.seed(21)
xm <- shape_df(
  par.drm = list(a = c(rep(1, 3), rlnorm(3, 0, .3), rlnorm(4, 0, .3)), b = rnorm(10), g = c(rep(NA, 6), runif(4, .1, .25))),
  par.prm = list(a = rlnorm(3, 0, .3), d = list(sort(rnorm(2)), sort(rnorm(3)), c(-1, .5, 0.2))),
  cats = c(rep(2, 10), 3, 4, 4), model = c(rep("1PLM", 3), rep("2PLM", 3), rep("3PLM", 4), "GRM", "GRM", "GPCM"))
dm <- simdat(x = xm, theta = rnorm(600), D = 1.702)
dm[sample(length(dm), 0.1 * length(dm))] <- NA
dm[1, ] <- NA

test_that("est_irt() log-likelihood, posterior, expected counts, and AIC/BIC match independent values", {
  skip_on_cran()
  fit <- suppressWarnings(est_irt(x = xm, data = dm, D = 1.702, se = FALSE, verbose = FALSE))
  ref <- indep_marg(fit$par.est, dm, fit$weights, D = 1.702)
  expect_equal(fit$loglikelihood, ref$loglik, tolerance = 1e-12)
  expect_equal(fit$posterior.dist, ref$post, tolerance = 1e-12)
  # free parameters: common 1PLM slope + 3 b, 3 x 2PLM, 4 x 3PLM, GRM 3 + 4, GPCM 4
  npar <- 1 + 3 + 6 + 12 + 3 + 4 + 4
  expect_equal(fit$npar.est, npar)
  expect_equal(fit$aic, 2 * npar - 2 * ref$loglik, tolerance = 1e-12)
  expect_equal(fit$bic, npar * log(nrow(dm)) - 2 * ref$loglik, tolerance = 1e-12)
  # expected counts of the E-step at the final estimates
  elm <- breakdown(fit$par.est)
  idx <- idxfinder(elm)
  dl <- divide_data(dm, idx, build_freqcat(dm, xm$cats))
  es <- Estep(elm, idx$idx.drm, idx$idx.prm, cbind(dl$data_drm_q, dl$data_drm_p), dl$data_prm, dl$data_all, fit$weights, D = 1.702)
  cols <- cols4item(nrow(xm), xm$cats, NULL)$cols.all
  for (j in seq_len(nrow(xm))) {
    rj <- sapply(0:(xm$cats[j] - 1), function(k) colSums(ref$post * (!is.na(dm[, j]) & dm[, j] == k)))
    expect_equal(unname(es$freq.exp[, cols[[j]]]), unname(rj), tolerance = 1e-10)
  }
})

test_that("est_irt(fipc = TRUE, fipc.method = 'MEM') reaches the fixed point of Kim (2006)", {
  skip_on_cran()
  fix <- c(1:2, 7:9, 11)
  fit <- suppressWarnings(est_irt(x = xm, data = dm, D = 1.702, fipc = TRUE, fix.loc = fix, Etol = 1e-6, se = FALSE, verbose = FALSE))
  # fixed items are unchanged
  expect_equal(unname(as.matrix(fit$par.est[fix, 4:5])), unname(as.matrix(xm[fix, 4:5])))
  # the prior is the normal distribution with the moments of the mean posterior
  ref <- indep_marg(fit$par.est, dm, fit$weights, D = 1.702)
  q <- fit$weights[, 1]
  pw <- colMeans(ref$post)
  mu <- sum(q * pw)
  s2 <- sum(q^2 * pw) - mu^2
  w <- stats::dnorm(q, mu, sqrt(s2))
  expect_equal(fit$weights[, 2], w / sum(w), tolerance = 1e-6)
  expect_equal(fit$loglikelihood, ref$loglik, tolerance = 1e-12)
})
