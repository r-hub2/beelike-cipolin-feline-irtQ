# Long tests: the likelihood of an examinee can be far below the smallest
# double; the E-step, the marginal log-likelihood, and pcd2() must stay finite

# 1200 low-discrimination 2PLM items: most examinees have log-likelihoods below -745
make_long <- function(n = 30, nit = 1200, seed = 41) {
  set.seed(seed)
  x <- shape_df(par.drm = list(a = rep(0.3, nit), b = rnorm(nit), g = rep(NA, nit)), cats = 2, model = "2PLM")
  list(x = x, data = simdat(x = x, theta = rnorm(n), D = 1))
}

test_that("likelihood() rescales underflowing rows and posterior() stays finite", {
  lt <- make_long()
  w <- gen.weight(dist = "norm", theta = seq(-6, 6, length.out = 49))
  ref <- indep_marg(lt$x, lt$data, w, D = 1)
  # the design must actually underflow
  expect_gt(sum(ref$loglik_i < -745), 10)
  elm <- breakdown(lt$x)
  idx <- idxfinder(elm)
  dl <- divide_data(lt$data, idx, build_freqcat(lt$data, lt$x$cats))
  lk <- likelihood(elm, idx.drm = idx$idx.drm, idx.prm = idx$idx.prm,
                   data_drm = cbind(dl$data_drm_q, dl$data_drm_p), data_prm = dl$data_prm,
                   theta = w[, 1], D = 1)
  post <- posterior(lk$L, w)
  expect_false(any(is.nan(post)))
  expect_equal(post, ref$post, tolerance = 1e-10)
  lli <- log(as.vector(lk$L %*% w[, 2])) + lk$offset
  expect_equal(lli, ref$loglik_i, tolerance = 1e-12)
})

test_that("est_irt() and est_irt(fipc = TRUE) report the correct log-likelihood for a long test", {
  skip_on_cran()
  lt <- make_long()
  fit <- suppressWarnings(est_irt(x = lt$x, data = lt$data, D = 1, MaxE = 1, se = FALSE, verbose = FALSE))
  ref <- indep_marg(fit$par.est, lt$data, fit$weights, D = 1)
  expect_true(is.finite(fit$loglikelihood))
  expect_equal(fit$loglikelihood, ref$loglik, tolerance = 1e-12)
  expect_equal(fit$posterior.dist, ref$post, tolerance = 1e-10)
  ff <- suppressWarnings(est_irt(x = lt$x, data = lt$data, D = 1, fipc = TRUE, fix.loc = 1:1195,
                                 MaxE = 2, se = FALSE, verbose = FALSE))
  ref2 <- indep_marg(ff$par.est, lt$data, ff$weights, D = 1)
  expect_equal(ff$loglikelihood, ref2$loglik, tolerance = 1e-12)
  expect_true(all(is.finite(unlist(ff$group.par[1, ]))))
})

test_that("pcd2() gives finite statistics for a long test", {
  lt <- make_long()
  pc <- pcd2(x = lt$x, data = lt$data, D = 1)$no_purify$ipd_stat$pcd2
  expect_false(any(is.nan(pc)))
  # independent pcd2 for the first three items (points without posterior mass add nothing)
  w <- gen.weight(dist = "norm", theta = seq(-6, 6, length.out = 49))
  post <- indep_marg(lt$x, lt$data, w, D = 1)$post
  ref <- sapply(1:3, function(j) {
    n <- colSums(post)
    r <- colSums(post * lt$data[, j])
    p <- indep_prob(w[, 1], "2PLM", c(0.3, lt$x$par.2[j]), 2, 1)[, 2]
    ok <- n > 0
    sum((n[ok] / nrow(lt$data)) * (r[ok] / n[ok] - p[ok])^2)
  })
  expect_equal(pc[1:3], ref, tolerance = 1e-10)
})
