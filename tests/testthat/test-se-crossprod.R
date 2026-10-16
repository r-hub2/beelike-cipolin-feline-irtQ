# Standard errors of est_irt() and est_item() compared with information
# matrices built from closed-form likelihoods and finite differences

prob_ref <- function(theta, model, par, D) {
  if (model %in% c("1PLM", "2PLM", "3PLM")) {
    g <- if (length(par) >= 3) par[3] else 0
    p <- g + (1 - g) / (1 + exp(-D * par[1] * (theta - par[2])))
    return(cbind(1 - p, p))
  }
  if (model == "GRM") {
    ps <- matrix(sapply(par[-1], function(dk) 1 / (1 + exp(-D * par[1] * (theta - dk)))), nrow = length(theta))
    return(cbind(1, ps) - cbind(ps, 0))
  }
  lz <- t(apply(outer(theta, c(0, par[-1]), "-") * D * par[1], 1, cumsum))
  e <- exp(lz - apply(lz, 1, max))
  e / rowSums(e)
}

# per-examinee log marginal likelihood
marg_ll_ref <- function(items, data, nodes, wts, D) {
  LL <- matrix(0, nrow(data), length(nodes))
  for (j in seq_along(items)) {
    lp <- log(prob_ref(nodes, items[[j]]$model, items[[j]]$par, D))
    ok <- !is.na(data[, j])
    LL[ok, ] <- LL[ok, ] + t(lp[, data[ok, j] + 1, drop = FALSE])
  }
  mx <- apply(LL, 1, max)
  log(as.vector(exp(LL - mx) %*% wts)) + mx
}

# cross-product information of the per-examinee marginal scores
xpd_ref <- function(items, free, data, nodes, wts, D, h = 1e-6) {
  v0 <- mapply(function(j, k) items[[j]]$par[k], free$item, free$k)
  setv <- function(v) {
    for (r in seq_along(v)) items[[free$item[r]]]$par[free$k[r]] <- v[r]
    items
  }
  S <- sapply(seq_along(v0), function(r) {
    e <- rep(0, length(v0))
    e[r] <- h
    (marg_ll_ref(setv(v0 + e), data, nodes, wts, D) - marg_ll_ref(setv(v0 - e), data, nodes, wts, D)) / (2 * h)
  })
  crossprod(S)
}

xt_se <- shape_df(
  par.drm = list(a = c(1.2, 0.8, 1.3), b = c(-0.5, 0.3, 0.1), g = c(0.2, 0.15, NA)),
  par.prm = list(a = c(1.1, 0.9), d = list(c(-1, 0, 1.1), c(-0.4, 0.6))),
  item.id = paste0("I", 1:5), cats = c(2, 2, 2, 4, 3), model = c("3PLM", "3PLM", "2PLM", "GRM", "GPCM")
)

se_check <- function(D, fix.g, g.val) {
  set.seed(5)
  resp <- simdat(xt_se, theta = rnorm(1000), D = D)
  fit <- est_irt(
    data = resp, D = D, model = xt_se$model, cats = xt_se$cats, fix.g = fix.g, g.val = g.val,
    use.gprior = TRUE, verbose = FALSE
  )
  pe <- fit$par.est
  items <- lapply(1:5, function(j) {
    v <- as.numeric(pe[j, 4:ncol(pe)])
    list(model = pe$model[j], par = v[!is.na(v)])
  })
  if (fix.g) for (j in 1:2) items[[j]]$par[3] <- g.val
  free <- do.call(rbind, lapply(1:5, function(j) {
    k <- switch(pe$model[j], "3PLM" = if (fix.g) 1:2 else 1:3, "2PLM" = 1:2, GRM = 1:4, GPCM = 1:3)
    data.frame(item = j, k = k)
  }))
  I <- xpd_ref(items, free, as.matrix(resp), fit$weights[, 1], fit$weights[, 2], D)
  # information of the Beta(5, 16) prior of the estimated guessing parameters
  if (!fix.g) {
    for (r in which(free$k == 3 & free$item <= 2)) {
      g <- items[[free$item[r]]]$par[3]
      I[r, r] <- I[r, r] + 4 / g^2 + 15 / (1 - g)^2
    }
  }
  list(
    pkg = mapply(function(j, k) fit$se.est[j, 3 + k], free$item, free$k),
    ref = sqrt(diag(solve(I)))
  )
}

test_that("est_irt() cross-product SEs match marginal scores with D = 1.702", {
  skip_on_cran()
  rst <- se_check(1.702, FALSE, 0.2)
  expect_equal(unname(rst$pkg), unname(rst$ref), tolerance = 1e-4)
})

test_that("est_irt() cross-product SEs use the fixed guessing value", {
  rst <- se_check(1, TRUE, 0.25)
  expect_equal(unname(rst$pkg), unname(rst$ref), tolerance = 1e-4)
})

test_that("est_irt() SEs are invariant to the scaling constant for rescaled slopes", {
  set.seed(8)
  resp <- simdat(xt_se, theta = rnorm(1500), D = 1)
  ap <- list(dist = "lnorm", params = c(0, 0.5))
  fit_d1 <- est_irt(
    data = resp, D = 1, model = xt_se$model, cats = xt_se$cats, use.gprior = TRUE, use.aprior = TRUE,
    aprior = ap, Etol = 1e-8, MaxE = 500, verbose = FALSE
  )
  fit_d17 <- est_irt(
    data = resp, D = 1.702, model = xt_se$model, cats = xt_se$cats, use.gprior = TRUE, use.aprior = TRUE,
    aprior = ap, Etol = 1e-8, MaxE = 500, verbose = FALSE
  )

  # the slope estimates and their SEs shrink by the factor 1.702, so the rescaled values agree with D = 1
  par_d17 <- as.matrix(fit_d17$par.est[, -(1:3)])
  par_d17[, 1] <- par_d17[, 1] * 1.702
  expect_equal(par_d17, as.matrix(fit_d1$par.est[, -(1:3)]), tolerance = 1e-5, ignore_attr = TRUE)
  se_d17 <- as.matrix(fit_d17$se.est[, -(1:3)])
  se_d17[, 1] <- se_d17[, 1] * 1.702
  expect_equal(se_d17, as.matrix(fit_d1$se.est[, -(1:3)]), tolerance = 1e-5, ignore_attr = TRUE)
})

test_that("est_item() SEs of 1PLM items with a common slope match the observed information", {
  x0 <- shape_df(par.drm = list(a = rep(1, 5), b = rep(0, 5), g = NA), cats = 2, model = "2PLM")
  sc <- est_score(x0, LSAT6, D = 1, method = "EAP")$est.theta
  for (D in c(1, 1.702)) {
    fit <- est_item(data = LSAT6, score = sc, D = D, model = "1PLM", cats = 2, fix.a.1pl = FALSE, verbose = FALSE)
    p <- c(fit$par.est$par.1[1], fit$par.est$par.2)
    nll <- function(v) {
      -sum(sapply(1:5, function(j) {
        P <- prob_ref(sc, "2PLM", c(v[1], v[j + 1]), D)[, 2]
        sum(LSAT6[, j] * log(P) + (1 - LSAT6[, j]) * log(1 - P))
      }))
    }

    # central finite-difference Hessian of the negative log-likelihood
    h <- 1e-4
    H <- matrix(0, 6, 6)
    for (i in 1:6) {
      for (j in 1:6) {
        ei <- rep(0, 6)
        ei[i] <- h
        ej <- rep(0, 6)
        ej[j] <- h
        H[i, j] <- (nll(p + ei + ej) - nll(p + ei - ej) - nll(p - ei + ej) + nll(p - ei - ej)) / (4 * h^2)
      }
    }
    expect_equal(c(fit$se.est$par.1[1], fit$se.est$par.2), sqrt(diag(solve(H))), tolerance = 1e-4)
  }
})
