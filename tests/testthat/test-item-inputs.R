# input checks of est_item() and the items it leaves out

set.seed(7)
xi <- shape_df(par.drm = list(a = c(1, 1, 1, rlnorm(5, 0, .3)), b = rnorm(8), g = c(rep(NA, 5), .15, .2, .1)),
               par.prm = list(a = rlnorm(2, 0, .3), d = list(sort(rnorm(2)), sort(rnorm(3)))),
               cats = c(rep(2, 8), 3, 4), model = c(rep("1PLM", 3), "2PLM", "2PLM", rep("3PLM", 3), "GRM", "GPCM"))
th <- rnorm(600)
di <- simdat(x = xi, theta = th, D = 1.702)
di[sample(length(di), 0.08 * length(di))] <- NA

test_that("score must give one finite value per examinee; matrix and data frame scores work", {
  expect_error(est_item(x = xi, data = di, score = replace(th, 3, NA), verbose = FALSE), "one finite ability value")
  expect_error(est_item(x = xi, data = di, score = th[-1], verbose = FALSE), "one finite ability value")
  f0 <- est_item(x = xi, data = di, score = th, D = 1.702, verbose = FALSE)
  f1 <- est_item(x = xi, data = di, score = data.frame(th = th), D = 1.702, verbose = FALSE)
  expect_identical(f1$par.est, f0$par.est)
  expect_equal(unname(f1$group.par), c(mean(th), sd(th)))
  f2 <- est_item(x = xi, data = di, score = matrix(th, ncol = 1), D = 1.702, verbose = FALSE)
  expect_identical(f2$par.est, f0$par.est)
})

test_that("responses must be scored 0, ..., cats - 1", {
  d2 <- di
  d2[1, 9] <- 3
  expect_error(est_item(x = xi, data = d2, score = th, verbose = FALSE), "scored as 0, 1")
  # responses stored as character strings are checked by their values
  dc <- as.data.frame(di)
  dc[, 4] <- as.character(dc[, 4])
  f0 <- est_item(x = xi, data = di, score = th, D = 1.702, verbose = FALSE)
  fc <- est_item(x = xi, data = dc, score = th, D = 1.702, verbose = FALSE)
  expect_equal(fc$par.est, f0$par.est)
  dc[1, 4] <- "a"
  expect_error(est_item(x = xi, data = dc, score = th, verbose = FALSE), "scored as 0, 1")
})

test_that("a single item left after removing the items with no responses is estimated", {
  d4 <- di
  d4[, 2:10] <- NA
  f <- suppressWarnings(est_item(x = xi, data = d4, score = th, D = 1.702, verbose = FALSE))
  expect_equal(f$deleted.item, 2:10)
  expect_true(all(!is.na(f$par.est[1, 4:5])))
  expect_true(all(is.na(f$par.est[-1, 4])))
})

test_that("the message about items that did not converge uses the positions in the supplied data", {
  d5 <- di
  d5[, 3] <- NA
  f <- suppressWarnings(est_item(x = xi, data = d5, score = th, D = 1.702, control = list(iter.max = 1, eval.max = 1),
                                 verbose = FALSE))
  listed <- as.integer(regmatches(f$convergence, gregexpr("[0-9]+", f$convergence))[[1]])
  expect_false(3 %in% listed)
  expect_true(10 %in% listed)
})
