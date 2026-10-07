# fixtures: the simMST 1-3-3 panel (dichotomous 3PLM items on the D = 1.702 scale)

x_mst     <- simMST$item_bank
mod_mst   <- simMST$module
map_mst   <- simMST$route_map
cut_mst   <- simMST$cut_score

# responses of 30 examinees to all 56 items of the bank
set.seed(2027)
theta_mst <- rnorm(30)
resp_mst  <- simdat(x = x_mst, theta = theta_mst, D = 1.702)

# fixtures: a small mixed-format 1-2-2 panel (three 3PLM items and one GRM item
# per module), so that every module has the same maximum sum score of 5
set.seed(2028)
n_mod_mix <- 5L
x_mix <- shape_df(
  par.drm = list(a = runif(3 * n_mod_mix, 0.8, 1.6),
                 b = rnorm(3 * n_mod_mix, 0, 0.8),
                 g = rep(0.15, 3 * n_mod_mix)),
  par.prm = list(a = runif(n_mod_mix, 0.8, 1.4),
                 d = replicate(n_mod_mix, sort(rnorm(2, 0, 0.7)), simplify = FALSE)),
  cats    = c(rep(2L, 3 * n_mod_mix), rep(3L, n_mod_mix)),
  model   = c(rep("3PLM", 3 * n_mod_mix), rep("GRM", n_mod_mix))
)
# items 1-3 of module 1, 4-6 of module 2, ...; the GRM items are 16-20
mod_mix <- matrix(0L, nrow = nrow(x_mix), ncol = n_mod_mix)
for (m in seq_len(n_mod_mix)) {
  mod_mix[c(((m - 1L) * 3L + 1L):(m * 3L), 15L + m), m] <- 1L
}
map_mix <- matrix(0L, n_mod_mix, n_mod_mix)
map_mix[1, 2:3] <- 1L
map_mix[2:3, 4:5] <- 1L
set.seed(2029)
theta_mix <- rnorm(30)
resp_mix  <- simdat(x = x_mix, theta = theta_mix, D = 1)

# helper: item rows of the modules administered up to stage s for examinee i
items_upto <- function(path_row, module, s) {
  unlist(lapply(path_row[seq_len(s)], function(m) which(module[, m] == 1)))
}

# helper: max abs difference between the routing estimates of run_mst() and
# est_score() applied to the cumulative responses of stages 1..s
max_diff_cum <- function(fit, x, module, resp, D, method, stages = 1:2) {
  diffs <- vapply(seq_len(nrow(fit$path)), function(i) {
    max(vapply(stages, function(s) {
      it  <- items_upto(fit$path[i, ], module, s)
      ref <- est_score(x = x[it, ], data = resp[i, it, drop = FALSE],
                       D = D, method = method)$est.theta
      abs(ref - fit$theta.route[i, s])
    }, numeric(1L)))
  }, numeric(1L))
  max(diffs)
}


# 1. structure of the result

test_that("run_mst() returns the documented structure for every routing method", {
  for (meth in c("ML", "WL", "MLF", "MAP", "EAP", "EAP.SUM", "INV.TCC")) {
    fit <- run_mst(
      x = x_mst, route_map = map_mst, module = mod_mst,
      theta = theta_mst[1:10], response = resp_mst[1:10, ], D = 1.702,
      route_method = "bmat",
      route_score = list(method = meth),
      final_score = list(method = "ML"), verbose = FALSE
    )
    expect_s3_class(fit, "run_mst")
    expect_equal(dim(fit$theta.route), c(10L, 3L))
    expect_equal(dim(fit$path), c(10L, 3L))
    # routing estimates of stages 1 and 2 are finite
    expect_true(all(is.finite(fit$theta.route[, 1:2])))
    # the last column of theta.route is the final estimate
    expect_equal(unname(fit$theta.route[, 3]), fit$est.theta)
  }
})


# 2. routing estimates use the responses to all modules administered so far

test_that("ML and EAP routing estimates equal est_score() on the cumulative responses", {
  for (meth in c("ML", "EAP")) {
    fit <- run_mst(
      x = x_mst, route_map = map_mst, module = mod_mst,
      theta = theta_mst, response = resp_mst, D = 1.702,
      route_method = "bmat",
      route_score = list(method = meth),
      final_score = list(method = "ML"), verbose = FALSE
    )
    expect_lt(max_diff_cum(fit, x_mst, mod_mst, resp_mst, 1.702, meth), 1e-6)
  }
})

test_that("WL, MAP, and MLF routing estimates equal est_score() on the cumulative responses", {
  for (meth in c("WL", "MAP", "MLF")) {
    fit <- run_mst(
      x = x_mst, route_map = map_mst, module = mod_mst,
      theta = theta_mst, response = resp_mst, D = 1.702,
      route_method = NULL, cut_score = cut_mst,
      route_score = list(method = meth),
      final_score = list(method = "ML"), verbose = FALSE
    )
    expect_lt(max_diff_cum(fit, x_mst, mod_mst, resp_mst, 1.702, meth), 1e-6)
  }
})

test_that("the stage 2 routing estimate differs from the estimate based on the stage 2 module alone", {
  fit <- run_mst(
    x = x_mst, route_map = map_mst, module = mod_mst,
    theta = theta_mst, response = resp_mst, D = 1.702,
    route_method = "bmat",
    route_score = list(method = "ML"),
    final_score = list(method = "ML"), verbose = FALSE
  )
  # estimate from the responses to the stage 2 module only
  stage_only <- vapply(seq_len(nrow(fit$path)), function(i) {
    it <- which(mod_mst[, fit$path[i, 2]] == 1)
    est_score(x = x_mst[it, ], data = resp_mst[i, it, drop = FALSE],
              D = 1.702, method = "ML")$est.theta
  }, numeric(1L))
  expect_gt(max(abs(stage_only - fit$theta.route[, 2])), 1e-3)
})

test_that("routing estimates are cumulative with mixed-format modules", {
  for (meth in c("ML", "EAP")) {
    fit <- run_mst(
      x = x_mix, route_map = map_mix, module = mod_mix,
      theta = theta_mix, response = resp_mix, D = 1,
      route_method = "bmat",
      route_score = list(method = meth),
      final_score = list(method = "ML"), verbose = FALSE
    )
    expect_lt(max_diff_cum(fit, x_mix, mod_mix, resp_mix, 1, meth, stages = 1:2), 1e-6)
  }
})

test_that("routing estimates use only the observed responses when some are missing", {
  resp_na <- resp_mst
  # set about 10 percent of the responses to missing at fixed positions
  set.seed(2030)
  na_pos <- sample(length(resp_na), size = round(0.10 * length(resp_na)))
  resp_na[na_pos] <- NA
  fit <- run_mst(
    x = x_mst, route_map = map_mst, module = mod_mst,
    theta = theta_mst, response = resp_na, D = 1.702,
    route_method = "bmat",
    route_score = list(method = "ML"),
    final_score = list(method = "ML"), verbose = FALSE
  )
  # est_score() drops the missing responses and scores the remaining items of
  # the administered modules
  expect_lt(max_diff_cum(fit, x_mst, mod_mst, resp_na, 1.702, "ML"), 1e-6)
})

test_that("INV.TCC routing estimates equal the inverse TCC lookup of the cumulative sum score", {
  fit <- run_mst(
    x = x_mst, route_map = map_mst, module = mod_mst,
    theta = theta_mst, response = resp_mst, D = 1.702,
    route_method = NULL, cut_score = cut_mst,
    route_score = list(method = "INV.TCC"),
    final_score = list(method = "INV.TCC"), verbose = FALSE
  )
  for (i in seq_len(nrow(fit$path))) {
    for (s in 1:2) {
      it  <- items_upto(fit$path[i, ], mod_mst, s)
      tbl <- est_score(x = x_mst[it, ], data = resp_mst[i, it, drop = FALSE],
                       D = 1.702, method = "INV.TCC")$score.table
      ref <- tbl$est.theta[tbl$sum.score == sum(resp_mst[i, it])]
      expect_equal(unname(fit$theta.route[i, s]), ref, tolerance = 1e-8)
    }
  }
})

test_that("cut-score routing assigns the module from the cumulative estimate", {
  fit <- run_mst(
    x = x_mst, route_map = map_mst, module = mod_mst,
    theta = theta_mst, response = resp_mst, D = 1.702,
    route_method = NULL, cut_score = cut_mst,
    route_score = list(method = "ML"),
    final_score = list(method = "ML"), verbose = FALSE
  )
  # stage 2 module follows the stage 1 estimate
  rank2 <- findInterval(fit$theta.route[, 1], cut_mst[[1]]) + 1L
  expect_equal(unname(fit$path[, 2]), c(2L, 3L, 4L)[rank2])

  # stage 3 module follows the estimate that pools the stage 1 and stage 2
  # responses; the rank is capped at the number of modules reachable from the
  # stage 2 module
  rank3 <- findInterval(fit$theta.route[, 2], cut_mst[[2]]) + 1L
  expected3 <- vapply(seq_len(nrow(fit$path)), function(i) {
    reach <- which(map_mst[fit$path[i, 2], ] == 1)
    reach[min(rank3[i], length(reach))]
  }, numeric(1L))
  expect_equal(unname(fit$path[, 3]), expected3)
})


# 3. agreement with the recursion-based evaluation

test_that("run_mst() with cut scores and inverse TCC scoring approaches reval_mst()", {
  skip_on_cran()

  # simulation with 1000 examinees at each of five ability levels
  grid <- seq(-2, 2, 1)
  n_rep <- 1000L
  set.seed(2031)
  theta_rep <- rep(grid, each = n_rep)
  mc <- run_mst(
    x = x_mst, route_map = map_mst, module = mod_mst,
    theta = theta_rep, D = 1.702,
    route_method = NULL, cut_score = cut_mst,
    route_score = list(method = "INV.TCC", range.tcc = c(-7, 7)),
    final_score = list(method = "INV.TCC", range.tcc = c(-7, 7)),
    verbose = FALSE
  )
  rv <- reval_mst(
    x = x_mst, D = 1.702, route_map = map_mst, module = mod_mst,
    cut_score = cut_mst, theta = grid, range.tcc = c(-7, 7)
  )$eval.tb

  grp       <- factor(theta_rep)
  bias_mc   <- as.numeric(tapply(mc$est.theta - theta_rep, grp, mean))
  csem_mc   <- as.numeric(tapply(mc$est.theta, grp, stats::sd))

  # bias differs from the analytical value within four Monte Carlo standard errors
  expect_true(all(abs(bias_mc - rv$bias) < 4 * rv$csem / sqrt(n_rep)))
  # CSEM differs from the analytical value by less than 15 percent
  expect_true(all(abs(csem_mc - rv$csem) / rv$csem < 0.15))
})


# 4. input validation

test_that("run_mst() stops on invalid input", {
  expect_error(
    run_mst(x = x_mst, route_map = map_mst, module = mod_mst, verbose = FALSE),
    "At least one of"
  )
  expect_error(
    run_mst(x = x_mst, route_map = map_mst, module = mod_mst,
            theta = theta_mst, route_method = NULL, verbose = FALSE),
    "cut_score"
  )
  expect_error(
    run_mst(x = x_mst, route_map = map_mst, module = mod_mst,
            theta = theta_mst, route_score = list(method = "XYZ"),
            verbose = FALSE),
    "route_score"
  )
})
