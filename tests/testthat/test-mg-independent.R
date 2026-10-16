# est_mg() checked against an independent multiple-group marginal likelihood

M <- make_mg_data(seed = 11)
gn <- c("G1", "G2", "G3")

test_that("the group and overall log-likelihoods match an independent multiple-group computation", {
  fit <- suppressWarnings(est_mg(x = M$x, data = M$data, group.name = gn, free.group = 2:3, verbose = FALSE))
  ll <- indep_mg_loglik(fit$par.est$overall, M$data, fit$weights, D = 1)
  expect_equal(unname(unlist(fit$loglikelihood$group)), ll, tolerance = 1e-9)
  expect_equal(fit$loglikelihood$overall, sum(ll), tolerance = 1e-9)
  expect_equal(fit$aic, -2 * sum(ll) + 2 * fit$npar.est, tolerance = 1e-9)
  expect_equal(fit$bic, -2 * sum(ll) + fit$npar.est * log(900), tolerance = 1e-9)
  # the reference group keeps the normal N(0, 1) weights
  th <- fit$weights[[1]][, 1]
  expect_equal(fit$weights[[1]][, 2], indep_normw(th, 0, 1), tolerance = 1e-12)
  # a free group has normal weights at its reported mean and variance
  expect_equal(fit$weights[[2]][, 2], indep_normw(th, fit$group.par$G2[1, 1], fit$group.par$G2[1, 2]), tolerance = 1e-4)
  # npar: free item parameters plus two moments of each free group
  expect_equal(fit$npar.est, sum(fit$npar.detail))
  expect_equal(unname(fit$npar.detail["dist"]), 4L)
})

test_that("the standard errors match an independent cross-product information with the g prior", {
  skip_on_cran()
  fit <- suppressWarnings(est_mg(x = M$x[1:2], data = M$data[1:2], group.name = gn[1:2], free.group = 2, verbose = FALSE))
  xo <- fit$par.est$overall
  pos <- as.matrix(fit$pos.par[, -(1:3)])
  free <- which(!is.na(pos), arr.ind = TRUE)
  pv <- as.matrix(xo[, -(1:3)])
  th <- fit$weights[[1]][, 1]
  ll_i <- function(v) {
    pm <- pv
    pm[free] <- v
    xx <- xo
    xx[, -(1:3)] <- pm
    unlist(lapply(1:2, function(k) {
      indep_marg_i(xx[match(colnames(M$data[[k]]), xx$id), ], M$data[[k]], th, fit$weights[[k]][, 2], 1)
    }))
  }
  v0 <- pv[free]
  h <- 1e-5
  S <- sapply(seq_along(v0), function(k) {
    e <- replace(numeric(length(v0)), k, h)
    (ll_i(v0 + e) - ll_i(v0 - e)) / (2 * h)
  })
  info <- crossprod(S)
  for (j in which(xo$model == "3PLM")) {
    k <- which(free[, 1] == j & free[, 2] == 3)
    g <- pv[j, 3]
    info[k, k] <- info[k, k] + 4 / g^2 + 15 / (1 - g)^2
  }
  ord <- pos[free]
  se_ind <- numeric(length(ord))
  se_ind[ord] <- sqrt(diag(solve(info)))
  expect_equal(as.matrix(fit$se.est$overall[, -(1:3)])[free], se_ind[ord], tolerance = 1e-4)
})

test_that("two reference groups give the est_irt() results of the pooled data", {
  set.seed(5)
  x1 <- M$x[[1]]
  dA <- simdat(x = x1, theta = rnorm(300), D = 1)
  dB <- simdat(x = x1, theta = rnorm(300), D = 1)
  mg <- suppressWarnings(est_mg(x = list(x1, x1), data = list(dA, dB), group.name = c("A", "B"), Etol = 1e-4, verbose = FALSE))
  ir <- suppressWarnings(est_irt(x = x1, data = rbind(dA, dB), Etol = 1e-4, verbose = FALSE))
  expect_equal(mg$par.est$overall, ir$par.est, tolerance = 1e-10)
  expect_equal(mg$se.est$overall, ir$se.est, tolerance = 1e-8)
  expect_equal(mg$loglikelihood$overall, ir$loglikelihood, tolerance = 1e-9)
  expect_equal(mg$npar.est, ir$npar.est)
})

test_that("MG-FIPC keeps the fixed items and reports the independent log-likelihood", {
  skip_on_cran()
  fit <- suppressWarnings(est_mg(x = M$x, data = M$data, group.name = gn, free.group = 1:3, fipc = TRUE,
                                 fix.id = c("C1", "C2", "C3", "C5"), verbose = FALSE))
  xo <- fit$par.est$overall
  fx <- match(c("C1", "C2", "C3", "C5"), xo$id)
  ref <- M$xall[match(xo$id[fx], M$xall$id), ]
  expect_equal(unname(as.matrix(xo[fx, 4:5])), unname(as.matrix(ref[, 4:5])), tolerance = 1e-12)
  ll <- indep_mg_loglik(xo, M$data, fit$weights, D = 1)
  expect_equal(unname(unlist(fit$loglikelihood$group)), ll, tolerance = 1e-9)
})

test_that("MG-FIPC standard errors of the new items match an independent cross-product information", {
  skip_on_cran()
  fix <- c("C1", "C2", "C3", "C5")
  fit <- suppressWarnings(est_mg(x = M$x[1:2], data = M$data[1:2], group.name = gn[1:2], free.group = 1:2, fipc = TRUE,
                                 fix.id = fix, verbose = FALSE))
  xo <- fit$par.est$overall
  pos <- as.matrix(fit$pos.par[, -(1:3)])
  free <- which(!is.na(pos), arr.ind = TRUE)
  expect_false(any(xo$id[free[, 1]] %in% fix))
  pv <- as.matrix(xo[, -(1:3)])
  th <- fit$weights[[1]][, 1]
  ll_i <- function(v) {
    pm <- pv
    pm[free] <- v
    xx <- xo
    xx[, -(1:3)] <- pm
    unlist(lapply(1:2, function(k) {
      indep_marg_i(xx[match(colnames(M$data[[k]]), xx$id), ], M$data[[k]], th, fit$weights[[k]][, 2], 1)
    }))
  }
  v0 <- pv[free]
  h <- 1e-5
  S <- sapply(seq_along(v0), function(k) {
    e <- replace(numeric(length(v0)), k, h)
    (ll_i(v0 + e) - ll_i(v0 - e)) / (2 * h)
  })
  info <- crossprod(S)
  for (j in intersect(which(xo$model == "3PLM"), unique(free[, 1]))) {
    k <- which(free[, 1] == j & free[, 2] == 3)
    g <- pv[j, 3]
    info[k, k] <- info[k, k] + 4 / g^2 + 15 / (1 - g)^2
  }
  ord <- pos[free]
  se_ind <- numeric(length(ord))
  se_ind[ord] <- sqrt(diag(solve(info)))
  expect_equal(as.matrix(fit$se.est$overall[, -(1:3)])[free], se_ind[ord], tolerance = 1e-4)
})
