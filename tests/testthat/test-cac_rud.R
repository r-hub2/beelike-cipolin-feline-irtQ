# cac_rud() computes classification accuracy and consistency from ability estimates and standard errors

set.seed(11)
theta_cac <- rnorm(50)
se_cac <- rep(0.3, 50)

test_that("cac_rud() labels the total row of the marginal table", {
  res <- cac_rud(cutscore = c(-0.5, 0.8), theta = theta_cac, se = se_cac)
  expect_equal(nrow(res$marginal), 4L)
  expect_identical(as.character(res$marginal$level), c("1", "2", "3", "marginal"))
  expect_false(anyNA(res$marginal$level))
})

test_that("cac_rud() marginal accuracy and consistency are the weighted sums of the conditional values", {
  res <- cac_rud(cutscore = c(-0.5, 0.8), theta = theta_cac, se = se_cac)
  total <- res$marginal[res$marginal$level == "marginal", ]
  expect_equal(total$accuracy, mean(res$conditional$accuracy))
  expect_equal(total$consistency, mean(res$conditional$consistency))
})

test_that("cac_rud() handles a performance level that no examinee reaches", {
  res <- cac_rud(cutscore = c(-0.5, 10), theta = theta_cac, se = se_cac)
  # the empty level has zero accuracy and consistency
  expect_equal(nrow(res$marginal), 4L)
  expect_equal(res$marginal$accuracy[3], 0)
  expect_equal(res$marginal$consistency[3], 0)
  # the confusion matrix keeps one row and one column per level
  expect_equal(dim(res$confusion), c(3L, 3L))
  expect_equal(unname(res$confusion[3, ]), c(0, 0, 0))
  # the marginal values still add up across the occupied levels
  expect_equal(res$marginal$accuracy[4], sum(res$marginal$accuracy[1:3]))
})

test_that("cac_rud() handles an empty level when quadrature weights are supplied", {
  wts <- gen.weight(n = 41, dist = "norm", mu = 0, sigma = 1)
  res <- cac_rud(cutscore = c(-0.5, 10), weights = wts, se = rep(0.3, 41))
  expect_equal(dim(res$confusion), c(3L, 3L))
  expect_equal(nrow(res$marginal), 4L)
  expect_lt(res$marginal$accuracy[3], 1e-10)
})

test_that("cac_rud() matches hand-computed normal probabilities", {
  # theta = 0 and se = 1 with cut scores -1 and 1
  res <- cac_rud(cutscore = c(-1, 1), theta = 0, se = 1)
  p <- c(stats::pnorm(-1), stats::pnorm(1) - stats::pnorm(-1), 1 - stats::pnorm(1))
  expect_equal(unname(unlist(res$prob.level[, paste0("p.level.", 1:3)])), p)
  expect_equal(res$conditional$accuracy, p[2])
  expect_equal(res$conditional$consistency, sum(p^2))
  # a theta equal to a cut score belongs to the higher level
  res_tie <- cac_rud(cutscore = c(-1, 1), theta = 1, se = 0.5)
  expect_equal(res_tie$conditional$level, 3L)
  expect_equal(res_tie$conditional$accuracy, 0.5)
})

test_that("cac_rud() with item metadata matches the reference", {
  x <- simMST$item_bank[1:20, ]
  wts <- gen.weight(n = 21, dist = "norm", mu = 0, sigma = 1)
  res <- cac_rud(x = x, cutscore = c(-0.5, 0.8), weights = wts, D = 1.702)
  se <- 1 / sqrt(info(x = x, theta = wts[, 1], D = 1.702, tif = TRUE)$tif)
  ref <- ref_rud(c(-0.5, 0.8), wts[, 1], se, wts[, 2])
  tot <- res$marginal[res$marginal$level == "marginal", ]
  expect_equal(tot$accuracy, ref$acc, tolerance = 1e-12)
  expect_equal(tot$consistency, ref$con, tolerance = 1e-12)
})

test_that("cac_rud() excludes examinees with a missing theta or se with a warning", {
  set.seed(11)
  th <- rnorm(20)
  se <- runif(20, 0.2, 0.4)
  res <- cac_rud(cutscore = c(-0.5, 0.8), theta = th, se = se)
  expect_warning(
    res_na <- cac_rud(cutscore = c(-0.5, 0.8), theta = c(th, NA, 0.1), se = c(se, 0.3, NA)),
    "2 examinee"
  )
  expect_equal(res_na$marginal, res$marginal)
})

test_that("cac_rud() stops on unordered cut scores and non-positive standard errors", {
  expect_error(cac_rud(cutscore = c(1, -1), theta = c(0, 1), se = c(0.3, 0.3)), "ascending")
  expect_error(cac_rud(cutscore = 0, theta = c(0, 1), se = c(0.3, -0.3)), "positive")
  expect_error(cac_rud(cutscore = 0, theta = c(0, 1), se = c(0.3, 0)), "positive")
})

test_that("cac_rud() rescales weights that do not sum to one", {
  wts <- gen.weight(n = 21, dist = "norm", mu = 0, sigma = 1)
  se <- rep(0.3, 21)
  expect_no_warning(res1 <- cac_rud(cutscore = c(-0.5, 0.8), weights = wts, se = se))
  # frequencies instead of proportions give the same indices
  expect_warning(
    res2 <- cac_rud(cutscore = c(-0.5, 0.8), weights = cbind(wts[, 1], 100 * wts[, 2]), se = se),
    "rescaled"
  )
  expect_equal(res2$marginal, res1$marginal)
  expect_lte(res2$marginal$accuracy[4], 1)
})

test_that("cac_rud() stops on weights that are negative or have a non-positive sum", {
  wts <- gen.weight(n = 21, dist = "norm", mu = 0, sigma = 1)
  se <- rep(0.3, 21)
  # a negative weight with a sum of one
  wts_neg <- wts
  wts_neg[1, 2] <- -0.1
  wts_neg[2, 2] <- wts_neg[2, 2] + 1 - sum(wts_neg[, 2])
  expect_equal(sum(wts_neg[, 2]), 1)
  expect_error(cac_rud(cutscore = c(-0.5, 0.8), weights = wts_neg, se = se), "nonnegative")
  # a missing weight and a zero sum
  expect_error(cac_rud(cutscore = c(-0.5, 0.8), weights = cbind(wts[, 1], NA_real_), se = se),
               "weights")
  expect_error(cac_rud(cutscore = c(-0.5, 0.8), weights = cbind(wts[, 1], 0), se = se),
               "weights")
})
