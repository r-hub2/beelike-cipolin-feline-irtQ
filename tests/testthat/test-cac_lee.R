# cac_lee() computes Lee's (2010) classification accuracy and consistency;
# expected values come from an independent Lord-Wingersky recursion
# (helper-mst-reference.R) or are worked out by hand; irtQ bounds item
# probabilities away from 0 and 1 by 1e-10, so the reference values agree
# to about 1e-9

# fixtures: a mixed-format test with 3PLM, GRM, and GPCM items
set.seed(101)
x_lee <- shape_df(
  par.drm = list(a = runif(12, 0.7, 1.8), b = rnorm(12), g = rep(0.2, 12)),
  par.prm = list(a = c(1.1, 0.9, 1.3, 0.8),
                 d = list(c(-1, 0.2, 1.1), c(-0.5, 0.6), c(-1.2, 0, 0.9), c(-0.3, 0.8))),
  cats = c(rep(2L, 12), 4L, 3L, 4L, 3L),
  model = c(rep("3PLM", 12), "GRM", "GPCM", "GPCM", "GRM")
)
wts_lee <- gen.weight(n = 31, dist = "norm", mu = 0, sigma = 1)
set.seed(102)
theta_lee <- rnorm(80)


test_that("cac_lee() matches a hand-computed two-item example", {
  # two 2PLM items with b = 0 at theta = 0: the summed scores 0, 1, 2 have
  # probabilities 0.25, 0.5, 0.25 and the true score is 1
  x2 <- shape_df(par.drm = list(a = c(1, 1), b = c(0, 0), g = c(0, 0)),
                 cats = 2, model = "2PLM")
  res <- cac_lee(x = x2, cutscore = 1, theta = 0)
  # a true score equal to the cut score belongs to the higher level
  expect_equal(as.integer(res$conditional$level), 2L)
  expect_equal(unname(unlist(res$prob.level[, c("p.level.1", "p.level.2")])), c(0.25, 0.75))
  expect_equal(res$conditional$accuracy, 0.75)
  expect_equal(res$conditional$consistency, 0.25^2 + 0.75^2)
})

test_that("cac_lee() with quadrature weights matches the reference for D = 1 and D = 1.702", {
  for (D in c(1, 1.702)) {
    res <- cac_lee(x = x_lee, cutscore = c(5, 10, 15), weights = wts_lee, D = D)
    ref <- ref_lee(x_lee, c(5, 10, 15), wts_lee[, 1], wts_lee[, 2], D = D)
    tot <- res$marginal[res$marginal$level == "marginal", ]
    expect_equal(tot$accuracy, ref$acc, tolerance = 1e-8)
    expect_equal(tot$consistency, ref$con, tolerance = 1e-8)
    expect_equal(res$conditional$accuracy, ref$cond_acc, tolerance = 1e-8)
    expect_equal(res$conditional$consistency, ref$cond_con, tolerance = 1e-8)
    expect_equal(unname(as.matrix(res$prob.level[, paste0("p.level.", 1:4)])), ref$P,
                 tolerance = 1e-8)
  }
})

test_that("cac_lee() with ability estimates and theta cut scores matches the reference", {
  res <- cac_lee(x = x_lee, cutscore = c(-1, 0, 1), theta = theta_lee, cut.obs = FALSE)
  # theta cut scores are converted to expected summed scores with the TCC
  cut_obs <- colSums(ref_lw(x_lee, c(-1, 0, 1)) * (0:(sum(x_lee$cats - 1))))
  expect_equal(res$cutscore, cut_obs, tolerance = 1e-8)
  ref <- ref_lee(x_lee, cut_obs, theta_lee, rep(1 / 80, 80))
  tot <- res$marginal[res$marginal$level == "marginal", ]
  expect_equal(tot$accuracy, ref$acc, tolerance = 1e-8)
  expect_equal(tot$consistency, ref$con, tolerance = 1e-8)
})

test_that("cac_lee() confusion matrix adds up to the weights of the true levels", {
  res <- cac_lee(x = x_lee, cutscore = c(5, 10, 15), weights = wts_lee)
  lev_mass <- tapply(wts_lee[, 2], res$conditional$level, sum)
  expect_equal(unname(rowSums(res$confusion)), as.vector(lev_mass), tolerance = 1e-6)
  expect_equal(sum(res$confusion), 1, tolerance = 1e-6)
  # the accuracy of each level is the diagonal of the confusion matrix
  expect_equal(res$marginal$accuracy[1:4], unname(diag(res$confusion)), tolerance = 1e-6)
})

test_that("cac_lee() keeps a level that contains no observed summed score", {
  # no integer score lies in [10.2, 10.7), so level 2 cannot be observed
  res <- cac_lee(x = x_lee, cutscore = c(10.2, 10.7), weights = wts_lee)
  expect_equal(dim(res$confusion), c(3L, 3L))
  expect_true(all(res$prob.level$p.level.2 == 0))
  expect_equal(unname(rowSums(res$prob.level[, paste0("p.level.", 1:3)])), rep(1, 31),
               tolerance = 1e-8)
  # a cut score above the maximum score leaves the top level empty
  res2 <- cac_lee(x = x_lee, cutscore = c(10, 40), weights = wts_lee)
  expect_true(all(res2$prob.level$p.level.3 == 0))
})

test_that("cac_lee() excludes missing ability estimates with a warning", {
  expect_warning(
    res_na <- cac_lee(x = x_lee, cutscore = c(5, 10, 15), theta = c(theta_lee, NA)),
    "missing"
  )
  res <- cac_lee(x = x_lee, cutscore = c(5, 10, 15), theta = theta_lee)
  expect_equal(res_na$marginal, res$marginal)
  expect_error(cac_lee(x = x_lee, cutscore = 5, theta = c(NA_real_, NA_real_)), "missing")
})

test_that("cac_lee() stops on cut scores that are not in ascending order", {
  expect_error(cac_lee(x = x_lee, cutscore = c(15, 5), weights = wts_lee), "ascending")
  expect_error(cac_lee(x = x_lee, cutscore = c(5, 5), weights = wts_lee), "ascending")
})

test_that("cac_lee() rescales weights that do not sum to one", {
  expect_warning(
    res2 <- cac_lee(x = x_lee, cutscore = c(5, 10), weights = cbind(wts_lee[, 1], 2 * wts_lee[, 2])),
    "rescaled"
  )
  res1 <- cac_lee(x = x_lee, cutscore = c(5, 10), weights = wts_lee)
  expect_equal(res2$marginal, res1$marginal)
})

test_that("cac_lee() does not change weights that sum to one", {
  expect_no_warning(res <- cac_lee(x = x_lee, cutscore = c(5, 10), weights = wts_lee))
  expect_equal(res$conditional$weights, wts_lee[, 2])
  expect_error(cac_lee(x = x_lee, cutscore = 5, weights = cbind(wts_lee[, 1], NA_real_)),
               "weights")
  # a negative weight stops even when the weights sum to one
  wts_neg <- wts_lee
  wts_neg[1, 2] <- -0.1
  wts_neg[2, 2] <- wts_neg[2, 2] + sum(wts_lee[, 2]) - sum(wts_neg[, 2])
  expect_equal(sum(wts_neg[, 2]), 1)
  expect_error(cac_lee(x = x_lee, cutscore = 5, weights = wts_neg), "nonnegative")
})
