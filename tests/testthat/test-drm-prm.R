test_that("drm() returns a matrix with correct dimensions", {
  P <- drm(theta = c(-1, 0, 1), a = c(1, 1.5), b = c(0, 0.5), g = c(0, 0.2), D = 1)
  expect_true(is.matrix(P))
  expect_equal(dim(P), c(3L, 2L))
})

test_that("drm() probabilities are bounded between g and 1", {
  g <- c(0.1, 0.2)
  P <- drm(theta = seq(-4, 4, by = 0.5), a = c(1, 1.5), b = c(0, 0.5), g = g, D = 1)
  expect_true(all(P > 0))
  expect_true(all(P < 1))
  # each column's minimum must be >= its guessing parameter
  expect_true(all(P[, 1] >= g[1]))
  expect_true(all(P[, 2] >= g[2]))
})

test_that("drm() 2PL: g=0 gives monotone increasing probabilities", {
  theta <- c(-2, -1, 0, 1, 2)
  P <- drm(theta = theta, a = 1.2, b = 0, D = 1)
  expect_true(all(diff(P[, 1]) > 0))
})

test_that("drm() 1PL: P(theta=b) is 0.5 when g=0", {
  P <- drm(theta = 0.5, a = 1, b = 0.5, g = 0, D = 1)
  expect_equal(as.numeric(P), 0.5, tolerance = 1e-6)
})

test_that("drm() 3PL: P(theta=b) equals (1+g)/2", {
  g_val <- 0.2
  b_val <- 1.0
  P <- drm(theta = b_val, a = 1, b = b_val, g = g_val, D = 1)
  expect_equal(as.numeric(P), (1 + g_val) / 2, tolerance = 1e-6)
})

test_that("drm() with NULL g defaults to 0 (2PL behaviour)", {
  P_null_g <- drm(theta = 0, a = 1, b = 0, D = 1)
  P_zero_g <- drm(theta = 0, a = 1, b = 0, g = 0, D = 1)
  expect_equal(as.numeric(P_null_g), as.numeric(P_zero_g), tolerance = 1e-10)
})

test_that("drm() single theta scalar input works", {
  P <- drm(theta = 0, a = 1, b = 0, g = 0.2, D = 1)
  expect_equal(dim(P), c(1L, 1L))
})

# ---- prm() ----------------------------------------------------------------------

test_that("prm() GRM: category probabilities sum to 1", {
  P <- prm(theta = c(-1, 0, 1), a = 1.2, d = c(-1, 0, 1), D = 1, pr.model = "GRM")
  expect_equal(dim(P), c(3L, 4L))
  row_sums <- rowSums(P)
  expect_equal(row_sums, rep(1, 3), tolerance = 1e-8)
})

test_that("prm() GPCM: category probabilities sum to 1", {
  P <- prm(theta = c(-2, 0, 2), a = 1.4, d = c(-0.2, 0, 0.5), D = 1, pr.model = "GPCM")
  expect_equal(dim(P), c(3L, 4L))
  row_sums <- rowSums(P)
  expect_equal(row_sums, rep(1, 3), tolerance = 1e-8)
})

test_that("prm() GRM: all probabilities strictly positive", {
  P <- prm(theta = seq(-3, 3, by = 1), a = 1, d = c(-1.5, 0, 1.5), D = 1, pr.model = "GRM")
  expect_true(all(P > 0))
})

test_that("prm() GPCM: all probabilities strictly positive", {
  P <- prm(theta = seq(-3, 3, by = 1), a = 1, d = c(-0.5, 0.5), D = 1, pr.model = "GPCM")
  expect_true(all(P > 0))
})

test_that("prm() GRM five-category item returns 5 columns", {
  P <- prm(theta = 0, a = 1.2, d = c(-1.5, -0.5, 0.5, 1.5), D = 1, pr.model = "GRM")
  expect_equal(ncol(P), 5L)
})

test_that("prm() GPCM and GRM give different results for same parameters", {
  theta <- c(-1, 0, 1)
  d     <- c(-0.5, 0.5)
  P_grm  <- prm(theta, a = 1, d = d, D = 1, pr.model = "GRM")
  P_gpcm <- prm(theta, a = 1, d = d, D = 1, pr.model = "GPCM")
  expect_false(isTRUE(all.equal(P_grm, P_gpcm)))
})

test_that("drm() known probability: 2PL D=1.702 matches hand calculation", {
  # P = 1 / (1 + exp(-1.702 * 1 * (0 - 0))) = 0.5
  P <- drm(theta = 0, a = 1, b = 0, g = 0, D = 1.702)
  expect_equal(as.numeric(P), 0.5, tolerance = 1e-6)
})

test_that("prm() uses the GRM when pr.model is not given", {
  P1 <- prm(theta = c(-1, 0, 1), a = 1.2, d = c(-1, 0.5), D = 1)
  P2 <- prm(theta = c(-1, 0, 1), a = 1.2, d = c(-1, 0.5), D = 1, pr.model = "GRM")
  expect_identical(P1, P2)

  # model names are matched without regard to case
  P3 <- prm(theta = c(-1, 0, 1), a = 1.2, d = c(-1, 0.5), D = 1, pr.model = "gpcm")
  P4 <- prm(theta = c(-1, 0, 1), a = 1.2, d = c(-1, 0.5), D = 1, pr.model = "GPCM")
  expect_identical(P3, P4)
})

test_that("prm() errors on an unknown pr.model", {
  expect_error(prm(theta = 0, a = 1, d = c(0, 1), D = 1, pr.model = "PCM"))
  expect_error(prm(theta = 0, a = 1, d = c(0, 1), D = 1, pr.model = "xyz"))
})

test_that("prm() GPCM probabilities at one theta do not depend on the other theta values", {
  d <- c(-1, 0, 1, 2)
  p1 <- prm(theta = 0, a = 2, d = d, D = 1.702, pr.model = "GPCM")
  p2 <- prm(theta = c(0, 50), a = 2, d = d, D = 1.702, pr.model = "GPCM")[1, , drop = FALSE]
  expect_equal(p2, p1, tolerance = 1e-12)
  expect_equal(unname(rowSums(prm(c(-50, 50), 1, d, 1, "GPCM"))), c(1, 1))
})

test_that("prm() GPCM probabilities stay above zero for very large slopes", {
  th <- seq(-6, 6, length.out = 49)
  P <- prm(th, a = 100, d = c(-1, 0, 1), D = 1.702, pr.model = "GPCM")
  expect_true(all(P > 0))
  expect_true(all(is.finite(log(P))))

  # the log-likelihood from these probabilities is finite
  ll <- irtQ:::loglike_prm(
    item_par = c(100, -1, 0, 1), r_i = matrix(1, 49, 4), theta = th,
    pr.mod = "GPCM", D = 1.702
  )
  expect_true(is.finite(ll))
})

test_that("prm() GPCM matches the closed-form probabilities", {
  th <- seq(-4, 4, 0.5)
  d <- c(-1, 0.2, 1.5)
  a <- 1.3
  # numerators exp(sum_{v <= k} D a (theta - b_v)) with b_0 = 0
  num <- sapply(0:3, function(k) {
    exp(vapply(th, function(t) sum(1.702 * a * (t - c(0, d)[seq_len(k + 1)])), numeric(1)))
  })
  ref <- num / rowSums(num)
  expect_equal(prm(th, a, d, 1.702, "GPCM"), unname(ref), tolerance = 1e-12)
})

test_that("drm() applies the lower bound with each item's own guessing parameter", {
  # a 2PLM item next to a 3PLM item keeps the exact logistic probabilities
  P <- drm(theta = c(-3, 0), a = c(1, 1), b = c(0, 0), g = c(0.2, 0), D = 1)
  expect_lt(max(abs(P[, 2] - plogis(c(-3, 0)))), 1e-12)

  # a 3PLM item far below its location stays above its guessing parameter
  P <- drm(theta = c(0, -50), a = c(1, 1), b = c(0, 0), g = c(0.2, 0), D = 1)
  expect_gt(P[2, 1], 0.2)
})

test_that("drm() matches the 1PLM, 2PLM, and 3PLM formulas on a grid", {
  th <- seq(-4, 4, 0.25)
  a <- c(0.8, 1.5, 2)
  b <- c(-1, 0, 1.2)
  g <- c(0.2, 0, 0.15)
  for (D in c(1, 1.702)) {
    ref <- sapply(1:3, function(j) g[j] + (1 - g[j]) / (1 + exp(-D * a[j] * (th - b[j]))))
    expect_lt(max(abs(drm(th, a, b, g, D) - ref)), 1e-12)
  }
})

test_that("drm() treats NA guessing parameters as zeros", {
  P <- drm(theta = c(-1, 0, 1), a = c(1, 1), b = c(0, 0), g = c(NA, 0.2), D = 1)
  expect_equal(P[, 1], plogis(c(-1, 0, 1)), tolerance = 1e-12)
  expect_equal(P[, 2], 0.2 + 0.8 * plogis(c(-1, 0, 1)), tolerance = 1e-12)
  expect_equal(as.numeric(drm(0, a = 1, b = 0, g = NA, D = 1)), 0.5)
})

test_that("prm() matches the GRM and GPCM formulas (b_v = beta - tau_v)", {
  # hand formulas used as references
  grm_ref <- function(th, a, d, D) {
    s <- matrix(sapply(d, function(dk) plogis(D * a * (th - dk))), nrow = length(th))
    cbind(1, s) - cbind(s, 0)
  }
  gpcm_ref <- function(th, a, d, D) {
    num <- sapply(0:length(d), function(k) if (k == 0) rep(0, length(th)) else D * a * (k * th - sum(d[1:k])))
    num <- matrix(num, nrow = length(th))
    num <- exp(num - apply(num, 1, max))
    num / rowSums(num)
  }
  th <- seq(-4, 4, 0.5)
  for (D in c(1, 1.702)) {
    expect_lt(max(abs(prm(th, 1.3, c(-1, 0.2, 1.5), D, "GRM") - grm_ref(th, 1.3, c(-1, 0.2, 1.5), D))), 1e-12)
    expect_lt(max(abs(prm(th, 1.3, c(-1, 0.2, 1.5), D, "GPCM") - gpcm_ref(th, 1.3, c(-1, 0.2, 1.5), D))), 1e-12)
  }

  # symmetric steps around beta give a symmetric distribution at theta = beta
  P <- prm(0.5, 1, 0.5 - c(0.8, 0, -0.8), 1, "GPCM")
  expect_equal(c(P), c(0.1550128, 0.3449872, 0.3449872, 0.1550128), tolerance = 1e-6)
})
