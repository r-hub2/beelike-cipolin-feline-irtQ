# bring.mirt() reads mirt fits that mix dichotomous and polytomous items

# simulate responses to four dichotomous items and three four-category items
mixed_resp <- function() {
  set.seed(11)
  x <- shape_df(
    par.drm = list(a = c(1, 1.2, 0.8, 1.5), b = c(-1, 0, 0.5, 1), g = rep(0.2, 4)),
    par.prm = list(a = c(1, 1.3, 0.9), d = list(c(-1, 0, 1), c(-0.5, 0.3, 1.2), c(-1.2, -0.2, 0.8))),
    cats = c(2, 2, 2, 2, 4, 4, 4), model = c(rep("3PLM", 4), rep("GRM", 3)), item.id = paste0("I", 1:7)
  )
  resp <- as.data.frame(simdat(x = x, theta = rnorm(1000), D = 1))
  colnames(resp) <- paste0("I", 1:7)
  resp
}

# maximum difference between the category probabilities of irtQ and mirt
max_prob_diff <- function(fit, meta) {
  theta <- seq(-3, 3, 0.5)
  # the DRM label of the dichotomous items is read as the 3PLM, which warns
  trace <- suppressWarnings(traceline(meta, theta = theta, D = 1))
  max(vapply(seq_len(nrow(meta)), function(j) {
    prob_mirt <- mirt::probtrace(mirt::extract.item(fit, j), matrix(theta))
    max(abs(unname(prob_mirt) - unname(trace$prob.cats[[j]])))
  }, numeric(1)))
}

test_that("bring.mirt() reads mixed-format fits whatever the order of the item models", {
  skip_if_not_installed("mirt")
  resp <- mixed_resp()
  cases <- list(
    list(cols = c(5:7, 1:4), type = c(rep("graded", 3), rep("2PL", 4))),
    list(cols = c(5:7, 1:4), type = c(rep("gpcm", 3), rep("2PL", 4))),
    list(cols = c(5, 1, 6, 2, 7, 3, 4), type = c("graded", "3PL", "gpcm", "3PL", "graded", "3PL", "3PL")),
    list(cols = 1:7, type = rep("graded", 7)),
    list(cols = 1:7, type = rep("gpcm", 7))
  )
  for (cs in cases) {
    fit <- suppressMessages(mirt::mirt(resp[, cs$cols], 1, itemtype = cs$type, verbose = FALSE, TOL = 1e-4))
    meta <- bring.mirt(fit)$full_df
    # every parameter used by the item model is available
    npar <- ifelse(meta$cats == 2, 3, meta$cats)
    for (i in seq_len(nrow(meta))) {
      expect_false(anyNA(unlist(meta[i, 3 + seq_len(npar[i])])))
    }
    # two-category items are dichotomous items without guessing
    expect_true(all(meta$model[meta$cats == 2] == "DRM"))
    # the category probabilities agree with mirt
    expect_lt(max_prob_diff(fit, meta), 1e-10)
  }
})

test_that("bring.mirt() labels two-category items fitted with polytomous models as DRM", {
  skip_if_not_installed("mirt")
  resp <- mixed_resp()
  fit <- suppressMessages(mirt::mirt(resp, 1, itemtype = rep("graded", 7), verbose = FALSE, TOL = 1e-4))
  out <- bring.mirt(fit)
  expect_identical(out$full_df$model, c(rep("DRM", 4), rep("GRM", 3)))
  expect_identical(out$model$new, c(rep("DRM", 4), rep("GRM", 3)))
  expect_identical(out$model$old, rep("graded", 7))
  # the guessing parameter of a two-category item is zero
  expect_equal(out$full_df$par.3[1:4], rep(0, 4))
})

test_that("bring.mirt() moves the intercept of the gpcmIRT model into the thresholds", {
  skip_if_not_installed("mirt")
  resp <- mixed_resp()
  type <- c(rep("2PL", 4), rep("gpcmIRT", 3))
  start <- suppressMessages(mirt::mirt(resp, 1, itemtype = type, pars = "values"))
  start$value[start$name == "c"] <- 0.4
  start$est <- FALSE
  fit <- suppressMessages(mirt::mirt(resp, 1, itemtype = type, pars = start, TOL = NaN, verbose = FALSE))
  meta <- bring.mirt(fit)$full_df
  # thresholds b_v - c / a reproduce the mirt probabilities
  coefs <- mirt::coef(fit, simplify = TRUE)$items
  expect_equal(
    unname(unlist(meta[5, 5:7])),
    unname(unlist(coefs[5, c("b1", "b2", "b3")] - 0.4 / coefs[5, "a1"]))
  )
  expect_lt(max_prob_diff(fit, meta), 1e-10)
})

test_that("bring.mirt() stops for dichotomous items fitted with an unsupported model", {
  skip_if_not_installed("mirt")
  resp <- mixed_resp()
  type <- c(rep("4PL", 4), rep("graded", 3))
  # a loose tolerance may leave a convergence warning of mirt, which is not tested here
  fit <- suppressWarnings(suppressMessages(
    mirt::mirt(resp, 1, itemtype = type, verbose = FALSE, TOL = 1e-3)
  ))
  expect_error(bring.mirt(fit), "only Rasch, 2PL, 3PL, graded, gpcm, gpcmIRT")
})

test_that("bring.mirt() stops for a model fitted to multiple groups", {
  skip_if_not_installed("mirt")
  resp <- mixed_resp()
  group <- rep(c("A", "B"), length.out = nrow(resp))
  fit <- suppressMessages(mirt::multipleGroup(resp, 1, group = group, verbose = FALSE))
  expect_error(bring.mirt(fit), "does not support multiple-group mirt objects")
  # a single group extracted from the fit is read
  expect_s3_class(bring.mirt(mirt::extract.group(fit, "A"))$full_df, "data.frame")
})
