# fixtures: the simMST 1-3-3 panel (dichotomous 3PLM items on the D = 1.702 scale)

x_rv     <- simMST$item_bank
module_rv <- simMST$module
map_rv   <- simMST$route_map
cut_rv   <- simMST$cut_score


test_that("reval_mst() evaluates the default theta grid and returns the documented objects", {
  rv <- reval_mst(x = x_rv, D = 1.702, route_map = map_rv, module = module_rv,
                  cut_score = cut_rv)
  expect_named(rv, c("panel.info", "item.by.mod", "item.by.path", "eq.theta",
                     "cdist.by.mod", "jdist.by.path", "eval.tb"))
  # the default theta grid is seq(-5, 5, 1)
  expect_equal(rv$eval.tb$theta, seq(-5, 5, 1))
  expect_true(all(c("theta", "mu", "sigma2", "bias", "csem") %in% names(rv$eval.tb)))
  expect_true(all(is.finite(rv$eval.tb$csem)))
})

test_that("reval_mst() stops when modules in a stage differ in maximum sum score", {
  # drop one item from module 2 so that stage 2 has modules of 7, 8, and 8 items
  module_bad <- module_rv
  module_bad[which(module_bad[, 2] == 1)[1], 2] <- 0L
  expect_error(
    reval_mst(x = x_rv, D = 1.702, route_map = map_rv, module = module_bad,
              cut_score = cut_rv, theta = seq(-1, 1, 1)),
    "same maximum sum score"
  )
})


# fixtures: a mixed-format 1-2-3 panel and its exact enumeration over the stage
# sum scores (helper-mst-reference.R)
pnl <- ref_panel()

test_that("reval_mst() equals an exact enumeration on a mixed-format 1-2-3 panel", {
  skip_on_cran()
  th <- c(-2, 0, 1.5)
  for (D in c(1, 1.702)) {
    rv <- reval_mst(x = pnl$x, D = D, route_map = pnl$route_map, module = pnl$module,
                    cut_score = pnl$cut_score, theta = th)$eval.tb
    bf <- ref_brute_mst(pnl$x, pnl$module, pnl$route_map, pnl$cut_score, th,
                        cfg = list(1L, 2:3, 4:6), D = D)
    expect_equal(bf$ptot, rep(1, 3), tolerance = 1e-12)
    expect_equal(rv$mu, bf$mu, tolerance = 1e-10)
    expect_equal(rv$sigma2, bf$sigma2, tolerance = 1e-10)
  }
})

test_that("reval_mst() equals an exact enumeration on the simMST panel", {
  skip_on_cran()
  th <- c(-1, 0.3)
  rv <- reval_mst(x = simMST$item_bank, D = 1.702, route_map = simMST$route_map,
                  module = simMST$module, cut_score = simMST$cut_score, theta = th)$eval.tb
  bf <- ref_brute_mst(simMST$item_bank, simMST$module, simMST$route_map,
                      simMST$cut_score, th, cfg = list(1L, 2:4, 5:7), D = 1.702)
  expect_equal(rv$mu, bf$mu, tolerance = 1e-10)
  expect_equal(rv$sigma2, bf$sigma2, tolerance = 1e-10)
})

test_that("reval_mst() accepts 2PLM items whose guessing parameter is missing", {
  x2 <- shape_df(par.drm = list(a = rep(1, 9), b = c(-0.5, 0, 0.5, rep(-1, 3), rep(1, 3)),
                                g = rep(0, 9)),
                 cats = 2, model = "2PLM")
  md <- matrix(0, 9, 3)
  md[1:3, 1] <- 1
  md[4:6, 2] <- 1
  md[7:9, 3] <- 1
  rm <- matrix(0, 3, 3)
  rm[1, 2:3] <- 1
  rv <- reval_mst(x2, route_map = rm, module = md, cut_score = list(0), theta = c(-1, 0))$eval.tb
  bf <- ref_brute_mst(x2, md, rm, list(0), c(-1, 0), cfg = list(1L, 2:3))
  expect_equal(rv$mu, bf$mu, tolerance = 1e-10)
  expect_equal(rv$sigma2, bf$sigma2, tolerance = 1e-10)
})

test_that("reval_mst() stops on cut scores that do not match the panel", {
  args <- list(x = simMST$item_bank, D = 1.702, route_map = simMST$route_map,
               module = simMST$module, theta = 0)
  # an extra cut score in a stage
  expect_error(do.call(reval_mst, c(args, list(cut_score = list(c(-0.45, 0.46, 2), c(-0.42, 0.45))))),
               "cut_score\\[\\[1\\]\\]")
  # cut scores in descending order
  expect_error(do.call(reval_mst, c(args, list(cut_score = list(c(0.46, -0.45), c(-0.42, 0.45))))),
               "ascending")
  # one vector too many
  expect_error(do.call(reval_mst, c(args, list(cut_score = c(simMST$cut_score, list(1))))),
               "length 2")
})

test_that("reval_mst() stops with a clear message when an inverse TCC estimate is missing", {
  expect_error(
    reval_mst(x = simMST$item_bank, D = 1.702, route_map = simMST$route_map,
              module = simMST$module, cut_score = simMST$cut_score, theta = 0,
              intpol = FALSE),
    "Inverse TCC estimates"
  )
})

test_that("reval_mst() maps the cut scores to the modules in ascending order of module index", {
  # module 2 reaches modules 5 and 6, and module 3 reaches modules 4 and 5, so
  # the order of first appearance in the pathways (5, 6, 4) differs from the
  # order of the module indices (4, 5, 6)
  rm_x <- matrix(0L, 6, 6)
  rm_x[1, 2:3] <- 1L
  rm_x[2, 5:6] <- 1L
  rm_x[3, 4:5] <- 1L
  th <- c(-1, 0.5)
  rv <- reval_mst(x = pnl$x, D = 1, route_map = rm_x, module = pnl$module,
                  cut_score = pnl$cut_score, theta = th)$eval.tb
  bf <- ref_brute_mst(pnl$x, pnl$module, rm_x, pnl$cut_score, th,
                      cfg = list(1L, 2:3, 4:6))
  expect_equal(rv$mu, bf$mu, tolerance = 1e-10)
  expect_equal(rv$sigma2, bf$sigma2, tolerance = 1e-10)
})

test_that("reval_mst() stops when the module matrix does not match the route map", {
  expect_error(
    reval_mst(x = simMST$item_bank, D = 1.702, route_map = simMST$route_map,
              module = simMST$module[, 1:6], cut_score = simMST$cut_score, theta = 0),
    "one column per module"
  )
})

test_that("reval_mst() does not depend on the index of the routing module", {
  x <- simMST$item_bank
  md <- simMST$module
  rm <- simMST$route_map
  # renumber the modules so that the routing module is module 4
  new_of_old <- c(4, 1, 2, 3, 5, 6, 7)
  md2 <- md[, order(new_of_old)]
  rm2 <- matrix(0, 7, 7)
  rm2[new_of_old, new_of_old] <- rm
  rv1 <- reval_mst(x, D = 1.702, route_map = rm, module = md,
                   cut_score = simMST$cut_score, theta = c(-1, 1))$eval.tb
  rv2 <- reval_mst(x, D = 1.702, route_map = rm2, module = md2,
                   cut_score = simMST$cut_score, theta = c(-1, 1))$eval.tb
  expect_equal(rv2, rv1)
})

test_that("reval_mst() stops when the first stage has more than one module", {
  # a 2-2 panel: modules 1 and 2 form stage 1 and modules 3 and 4 form stage 2
  rm_x <- matrix(0L, 4, 4)
  rm_x[1:2, 3:4] <- 1L
  expect_error(
    reval_mst(x = simMST$item_bank, D = 1.702, route_map = rm_x,
              module = simMST$module[, 1:4], cut_score = list(0), theta = 0),
    "single routing module"
  )
})

test_that("reval_mst() treats a NULL cut score element as an empty vector", {
  # a 1-1-2 panel: the cut score element of the transition to the single module of stage 2 is empty
  keep <- which(rowSums(pnl$module[, 1:4]) > 0)
  md <- pnl$module[keep, 1:4]
  rm_x <- matrix(0L, 4, 4)
  rm_x[1, 2] <- 1L
  rm_x[2, 3:4] <- 1L
  rv_null <- reval_mst(x = pnl$x[keep, ], D = 1, route_map = rm_x, module = md,
                       cut_score = list(NULL, 0.2), theta = c(-1, 1))$eval.tb
  rv_empty <- reval_mst(x = pnl$x[keep, ], D = 1, route_map = rm_x, module = md,
                        cut_score = list(numeric(0), 0.2), theta = c(-1, 1))$eval.tb
  expect_equal(rv_null, rv_empty)
  # a non-empty element for a single-module stage stops with a message that names the stage
  expect_error(
    reval_mst(x = pnl$x[keep, ], D = 1, route_map = rm_x, module = md,
              cut_score = list(c(0.1, 0.2), 0.2), theta = 0),
    "must be empty"
  )
})
