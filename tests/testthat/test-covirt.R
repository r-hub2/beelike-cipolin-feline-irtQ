# covirt() compared with the expected Fisher information computed from
# numerical derivatives of closed-form category probabilities

cat_prob_ref <- function(theta, model, par, D) {
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

# expected information sum_q w_q sum_k (dP_k)(dP_k)' / P_k for the free parameters
exp_info_ref <- function(model, par, free, D, nodes, wts, h = 1e-6) {
  P <- cat_prob_ref(nodes, model, par, D)
  dP <- lapply(free, function(j) {
    e <- rep(0, length(par))
    e[j] <- h
    (cat_prob_ref(nodes, model, par + e, D) - cat_prob_ref(nodes, model, par - e, D)) / (2 * h)
  })
  I <- matrix(0, length(free), length(free))
  for (r in seq_along(free)) {
    for (s in seq_along(free)) {
      I[r, s] <- sum(wts * rowSums(dP[[r]] * dP[[s]] / P))
    }
  }
  I
}

wq <- gen.weight(n = 41, dist = "norm", mu = 0, sigma = 1)

test_that("covirt() matches the expected information for all models", {
  cases <- list(
    list(m = "3PLM", par = c(1.3, 0.4, 0.18), free = 1:3),
    list(m = "2PLM", par = c(0.9, -0.6), free = 1:2),
    list(m = "GRM", par = c(1.2, -0.8, 0.9), free = 1:3),
    list(m = "GRM", par = c(1.2, -1.0, 0.2, 1.1), free = 1:4),
    list(m = "GRM", par = c(0.8, -1.5, -0.4, 0.5, 1.6), free = 1:5),
    list(m = "GPCM", par = c(1.1, -0.5, 0.7), free = 1:3),
    list(m = "GPCM", par = c(1.2, -1.0, 0.2, 1.1), free = 1:4),
    list(m = "GPCM", par = c(0.7, 0.4, -0.6, 1.0, 0.3), free = 1:5)
  )
  for (D in c(1, 1.702)) {
    for (cs in cases) {
      x <- data.frame(id = "I1", cats = max(2, length(cs$par)), model = cs$m)
      for (k in 1:5) x[[paste0("par.", k)]] <- if (k <= length(cs$par)) cs$par[k] else NA
      if (cs$m %in% c("2PLM", "3PLM")) x$cats <- 2
      ref <- solve(1500 * exp_info_ref(cs$m, cs$par, cs$free, D, wq[, 1], wq[, 2]))
      rst <- covirt(x, D = D, nstd = 1500)
      expect_equal(unname(rst$cov[[1]]), ref, tolerance = 1e-6)
      expect_equal(unname(rst$se[[1]]), sqrt(diag(ref)), tolerance = 1e-6)
    }
  }
})

test_that("covirt() handles PCM items through pcm.loc", {
  x <- data.frame(id = "P1", cats = 4, model = "GPCM", par.1 = 1, par.2 = -0.5, par.3 = 0.4, par.4 = 1.0)
  for (D in c(1, 1.702)) {
    ref <- solve(1000 * exp_info_ref("GPCM", c(1, -0.5, 0.4, 1.0), 2:4, D, wq[, 1], wq[, 2]))
    expect_equal(unname(covirt(x, D = D, nstd = 1000, pcm.loc = 1)$cov[[1]]), ref, tolerance = 1e-6)
  }
})

test_that("covirt() treats GRM and GPCM items with two categories as dichotomous items", {
  ref2 <- solve(1000 * exp_info_ref("2PLM", c(1.1, 0.3), 1:2, 1, wq[, 1], wq[, 2]))
  x <- data.frame(id = "G2", cats = 2, model = "GRM", par.1 = 1.1, par.2 = 0.3, par.3 = NA)
  expect_equal(unname(covirt(x, D = 1, nstd = 1000)$cov[[1]]), ref2, tolerance = 1e-6)
  x$model <- "GPCM"
  expect_equal(unname(covirt(x, D = 1, nstd = 1000)$cov[[1]]), ref2, tolerance = 1e-6)

  # a two-category PCM item has the location parameter only
  ref1 <- solve(1000 * exp_info_ref("2PLM", c(1, 0.4), 2, 1, wq[, 1], wq[, 2]))
  xp <- data.frame(id = "P2", cats = 2, model = "GPCM", par.1 = 1, par.2 = 0.4, par.3 = NA)
  expect_equal(unname(covirt(xp, D = 1, nstd = 1000, pcm.loc = 1)$cov[[1]]), ref1, tolerance = 1e-6)

  # mixed with a polytomous item
  x2 <- data.frame(
    id = c("G2", "T3"), cats = c(2, 3), model = c("GRM", "GRM"),
    par.1 = c(1.1, 0.9), par.2 = c(0.3, -0.5), par.3 = c(NA, 0.6)
  )
  rst <- covirt(x2, D = 1, nstd = 1000)
  expect_equal(unname(rst$cov[[1]]), ref2, tolerance = 1e-6)
  expect_equal(dim(rst$cov[[2]]), c(3L, 3L))
})

test_that("covirt() warns and returns named NA values when an information matrix is singular", {
  # a 2PLM item with a zero slope has no information about its difficulty
  x <- data.frame(
    id = c("S1", "S2"), cats = c(2, 2), model = c("2PLM", "2PLM"),
    par.1 = c(0, 1.1), par.2 = c(0.3, 0.2), par.3 = c(NA, NA)
  )
  expect_warning(rst <- covirt(x, D = 1, nstd = 1000), "The information matrix of item S1 is singular")
  expect_true(all(is.na(rst$cov$S1)))
  expect_identical(dimnames(rst$cov$S1), list(c("par.1", "par.2"), c("par.1", "par.2")))
  expect_identical(names(rst$se$S1), c("par.1", "par.2"))

  # the other item is not affected
  expect_false(anyNA(rst$cov$S2))
  expect_identical(dimnames(rst$cov$S2), dimnames(rst$cov$S1))
})
