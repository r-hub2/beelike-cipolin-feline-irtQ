# GRM thresholds must be in increasing order; otherwise grm() truncates the
# negative category probabilities and the probabilities no longer sum to one

# 10 3PLM + 5 five-category GRM items, D = 1.702, 300 examinees
make_mixed <- function(seed) {
  set.seed(1000 + seed)
  a3 <- rlnorm(10, 0, 0.3); b3 <- rnorm(10); g3 <- runif(10, .1, .25)
  ag <- rlnorm(5, 0, 0.3)
  bg <- t(sapply(1:5, function(i) sort(rnorm(4, 0, 1))))
  x <- shape_df(par.drm = list(a = a3, b = b3, g = g3),
                par.prm = list(a = ag, d = lapply(1:5, function(i) bg[i, ])),
                item.id = paste0("I", 1:15), cats = c(rep(2, 10), rep(5, 5)),
                model = c(rep("3PLM", 10), rep("GRM", 5)))
  set.seed(seed)
  list(x = x, data = simdat(x = x, theta = rnorm(300), D = 1.702))
}

grm_ordered <- function(pe) {
  idx <- which(pe$model == "GRM")
  all(sapply(idx, function(i) {
    b <- as.numeric(unlist(pe[i, 5:(3 + pe$cats[i])]))
    all(diff(b) > 0)
  }))
}

test_that("est_irt() returns ordered GRM thresholds when a category is rare", {
  skip_on_cran()
  # seed 5: item 14 has category counts 128/3/36/116/17
  md <- make_mixed(5)
  expect_equal(tabulate(md$data[, 14] + 1, 5), c(128, 3, 36, 116, 17))
  fit <- suppressWarnings(est_irt(x = md$x, data = md$data, D = 1.702, use.gprior = TRUE, se = FALSE, verbose = FALSE))
  expect_true(grm_ordered(fit$par.est))
  # the category probabilities of every GRM item sum to one
  for (i in which(fit$par.est$model == "GRM")) {
    P <- indep_prob(seq(-4, 4, by = 0.5), "GRM", as.numeric(unlist(fit$par.est[i, 4:8])), 5, 1.702)
    expect_true(all(P > 0))
    expect_equal(rowSums(P), rep(1, nrow(P)), tolerance = 1e-12)
  }
  # the reported log-likelihood is that of a valid model
  ref <- indep_marg(fit$par.est, md$data, fit$weights, D = 1.702)
  expect_equal(fit$loglikelihood, ref$loglik, tolerance = 1e-8)
})

test_that("est_irt() returns ordered GRM thresholds when a category is never observed", {
  skip_on_cran()
  # seed 3: item 15 has no response in category 1
  md <- make_mixed(3)
  expect_equal(tabulate(md$data[, 15] + 1, 5)[2], 0)
  fit <- suppressWarnings(est_irt(x = md$x, data = md$data, D = 1.702, use.gprior = TRUE, se = FALSE, verbose = FALSE))
  expect_true(grm_ordered(fit$par.est))
})

test_that("est_item() returns ordered GRM thresholds when a category is never observed", {
  set.seed(501)
  x <- shape_df(par.prm = list(a = rlnorm(5, 0, .3), d = lapply(1:5, function(i) sort(rnorm(4)))), cats = 5, model = "GRM")
  th <- rnorm(300)
  dat <- simdat(x = x, theta = th, D = 1.702)
  expect_equal(tabulate(dat[, 3] + 1, 5)[2], 0)
  fit <- suppressWarnings(est_item(x = x, data = dat, score = th, D = 1.702, verbose = FALSE))
  expect_true(grm_ordered(fit$par.est))
})

test_that("est_mg() returns ordered GRM thresholds when a category is rare", {
  skip_on_cran()
  # 8 3PLM and 4 five-category GRM items in two groups of 200 examinees
  set.seed(2001)
  x <- shape_df(
    par.drm = list(a = rlnorm(8, 0, 0.3), b = rnorm(8), g = runif(8, 0.1, 0.25)),
    par.prm = list(a = rlnorm(4, 0, 0.3), d = lapply(1:4, function(i) sort(rnorm(4)))),
    item.id = paste0("I", 1:12), cats = c(rep(2, 8), rep(5, 4)), model = c(rep("3PLM", 8), rep("GRM", 4)))
  set.seed(1)
  dat <- list(g1 = simdat(x = x, theta = rnorm(200), D = 1.702), g2 = simdat(x = x, theta = rnorm(200, 0.3), D = 1.702))
  fit <- suppressWarnings(est_mg(x = list(x, x), data = dat, group.name = c("G1", "G2"), free.group = "G2", D = 1.702,
                                 use.gprior = TRUE, se = FALSE, verbose = FALSE))
  expect_true(grm_ordered(fit$par.est$overall))
})
