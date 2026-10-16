# Numbers of free parameters used by AIC and BIC of est_irt() and est_mg()

xf <- bring.flexmirt(file = system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ"), "par")$Group1$full_df
x12 <- xf[c(1:10, 39:40), ]
set.seed(1)
d12 <- simdat(x = x12, theta = rnorm(500), D = 1)

# 10 3PLM items (3 parameters each) and 2 five-category GRM items (5 parameters each)
n.item <- 10 * 3 + 2 * 5
fit_irt <- function(...) suppressWarnings(est_irt(x = x12, data = d12, D = 1, MaxE = 5, se = FALSE, verbose = FALSE, ...))

# checks the definitions of npar.est, npar.detail, aic, and bic of a fitted object
check_npar <- function(fit, item, dist, weight, nobs) {
  expect_identical(names(fit$npar.detail), c("item", "dist", "weight"))
  expect_equal(unname(fit$npar.detail), c(item, dist, weight))
  expect_type(fit$npar.detail, "integer")
  expect_equal(sum(fit$npar.detail), fit$npar.est)
  expect_equal(fit$npar.est, item + dist + weight)
  ll <- if (is.list(fit$loglikelihood)) fit$loglikelihood$overall else fit$loglikelihood
  expect_equal(fit$aic, 2 * fit$npar.est - 2 * ll)
  expect_equal(fit$bic, fit$npar.est * log(nobs) - 2 * ll)
  # the counts are the last element of the object
  expect_identical(names(fit)[length(fit)], "npar.detail")
}

test_that("est_irt() counts the free parameters of the normal and empirical histogram priors", {
  skip_on_cran()
  check_npar(fit_irt(), n.item, 0, 0, nrow(d12))
  check_npar(fit_irt(EmpHist = TRUE), n.item, 0, 49 - 3, nrow(d12))
  w41 <- gen.weight(dist = "norm", theta = seq(-6, 6, length.out = 41))
  check_npar(fit_irt(weights = w41, EmpHist = TRUE), n.item, 0, 41 - 3, nrow(d12))
  check_npar(fit_irt(Quadrature = c(31, 5), EmpHist = TRUE), n.item, 0, 31 - 3, nrow(d12))
})

test_that("est_irt() counts the free parameters of FIPC without rescaling", {
  fix <- 1:6
  n.new <- n.item - (6 * 3)
  fipc_fit <- function(...) fit_irt(fipc = TRUE, ...)
  check_npar(fipc_fit(fix.loc = fix), n.new, 2, 0, nrow(d12))
  check_npar(fipc_fit(fix.loc = fix, EmpHist = TRUE), n.new, 0, 49 - 1, nrow(d12))
  check_npar(fipc_fit(fix.loc = fix, fipc.method = "OEM"), n.new, 2, 0, nrow(d12))
  check_npar(fipc_fit(fix.loc = seq_len(nrow(x12))), 0, 2, 0, nrow(d12))
  check_npar(fipc_fit(fix.loc = seq_len(nrow(x12)), EmpHist = TRUE), 0, 0, 49 - 1, nrow(d12))
})

test_that("est_mg() counts the free parameters of each group", {
  skip_on_cran()
  set.seed(2)
  dat <- list(g1 = simdat(x = x12, theta = rnorm(300), D = 1), g2 = simdat(x = x12, theta = rnorm(300, 0.5), D = 1),
              g3 = simdat(x = x12, theta = rnorm(300, 0.2), D = 1))
  mg_fit <- function(...) {
    suppressWarnings(est_mg(x = list(x12, x12, x12), data = dat, group.name = c("G1", "G2", "G3"), D = 1, MaxE = 5,
                            se = FALSE, verbose = FALSE, ...))
  }
  n.obs <- 900
  # one reference group and two free groups
  check_npar(mg_fit(free.group = c("G2", "G3")), n.item, 2 * 2, 0, n.obs)
  check_npar(mg_fit(free.group = c("G2", "G3"), EmpHist = TRUE), n.item, 0, 2 * (49 - 1) + 1 * (49 - 3), n.obs)
  # no free group: every group is a reference group
  check_npar(mg_fit(free.group = NULL, EmpHist = TRUE), n.item, 0, 3 * (49 - 3), n.obs)
  # FIPC with the first five items fixed in every group
  fipc_fit <- function(...) mg_fit(fipc = TRUE, fix.loc = list(1:5, 1:5, 1:5), ...)
  check_npar(fipc_fit(free.group = 1:3), n.item - 5 * 3, 3 * 2, 0, n.obs)
  check_npar(fipc_fit(free.group = 2:3, EmpHist = TRUE), n.item - 5 * 3, 0, 2 * (49 - 1) + (49 - 3), n.obs)
})

test_that("summary() shows the counts of the free parameters and skips them for older objects", {
  fit <- fit_irt(EmpHist = TRUE)
  out <- utils::capture.output(print(summary(fit)))
  i <- grep("Number of free parameters: 86", out)
  expect_length(i, 1)
  expect_identical(out[i + 1], paste0("   - Item parameters: ", n.item))
  expect_identical(out[i + 2], "   - Latent distribution parameters (mean and variance): 0")
  expect_identical(out[i + 3], "   - Latent distribution weights (empirical histogram): 46")

  # an object without npar.detail prints the total only
  old <- fit
  old$npar.detail <- NULL
  out.old <- utils::capture.output(print(summary(old)))
  expect_true(any(grepl("Number of free parameters: 86", out.old)))
  expect_false(any(grepl("Item parameters", out.old)))

  set.seed(2)
  dat <- list(g1 = simdat(x = x12, theta = rnorm(300), D = 1), g2 = simdat(x = x12, theta = rnorm(300, 0.5), D = 1))
  mg <- suppressWarnings(est_mg(x = list(x12, x12), data = dat, group.name = c("G1", "G2"), free.group = "G2", D = 1,
                                MaxE = 5, se = FALSE, verbose = FALSE))
  out.mg <- utils::capture.output(print(summary(mg)))
  j <- grep("Number of free parameters: 42", out.mg)
  expect_length(j, 1)
  expect_identical(out.mg[j + 1], paste0("   - Item parameters: ", n.item))
  expect_identical(out.mg[j + 2], "   - Latent distribution parameters (mean and variance): 2")
  expect_identical(out.mg[j + 3], "   - Latent distribution weights (empirical histogram): 0")
})
