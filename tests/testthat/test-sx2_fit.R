# sx2_fit() computes the S-X2 item fit statistic of Orlando and Thissen (2000)

x_lsat <- shape_df(
  par.drm = list(a = c(0.8, 0.7, 0.9, 0.6, 0.7), b = c(-3.4, -1.3, -0.3, -1.8, -2.8), g = rep(0, 5)),
  item.id = paste0("V", 1:5), cats = 2, model = "2PLM"
)

test_that("sx2_fit() stops for data that do not match the items or the score categories", {
  # a column is missing
  expect_error(sx2_fit(x_lsat, data = LSAT6[, 1:4], D = 1), "number of columns")

  # a response above the highest category of an item
  resp <- LSAT6
  resp[, 1] <- resp[, 1] + 1
  expect_error(sx2_fit(x_lsat, data = resp, D = 1), "outside the score categories")

  # a non-integer response
  resp <- LSAT6
  resp[3, 2] <- 0.5
  expect_error(sx2_fit(x_lsat, data = resp, D = 1), "outside the score categories")
})

test_that("sx2_fit() reads character and factor responses like numeric responses", {
  resp_chr <- matrix(as.character(as.matrix(LSAT6)), nrow(LSAT6))
  df_chr <- as.data.frame(resp_chr, stringsAsFactors = FALSE)
  df_fac <- as.data.frame(lapply(df_chr, factor))
  ref <- sx2_fit(x_lsat, data = LSAT6, D = 1)
  for (dat in list(resp_chr, df_chr, df_fac)) {
    fit <- sx2_fit(x_lsat, data = dat, D = 1)
    expect_identical(fit$fit_stat, ref$fit_stat)
    expect_identical(unname(fit$exp_freq), unname(ref$exp_freq))
    expect_identical(unname(fit$obs_freq), unname(ref$obs_freq))
  }

  # a response that cannot be read as a number stops
  resp_chr[1, 1] <- "a"
  expect_error(sx2_fit(x_lsat, data = resp_chr, D = 1), "numeric scores")
})

test_that("sx2_fit() replaces missing responses with zeros and names the items without any response", {
  resp <- LSAT6
  resp[1:5, 2] <- NA
  resp0 <- resp
  resp0[is.na(resp0)] <- 0
  expect_warning(f1 <- sx2_fit(x_lsat, data = resp, D = 1), "Missing responses are replaced with 0")
  expect_identical(f1$fit_stat, sx2_fit(x_lsat, data = resp0, D = 1)$fit_stat)

  # an item without any response is named in the warning
  resp[, 4] <- NA
  expect_warning(sx2_fit(x_lsat, data = resp, D = 1), "item\\(s\\) V4")
})

prm_file <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
x_full <- bring.flexmirt(file = prm_file, "par")$Group1$full_df

test_that("sx2_fit() does not count a summed score group twice in a short test", {
  # two 5-category items: the groups pooled at the two ends overlap
  x_short <- x_full[53:54, ]
  set.seed(12)
  resp <- simdat(x = x_short, theta = rnorm(1000), D = 1)
  fit <- sx2_fit(x_short, data = resp, min.collapse = 0)
  n_mid <- sum(rowSums(resp) > 0 & rowSums(resp) < 8)
  for (i in 1:2) {
    expect_equal(sum(fit$obs_freq[[i]]), n_mid)
    expect_equal(sum(fit$exp_freq[[i]]), n_mid)
  }

  # two dichotomous items and one 5-category item
  x_mix3 <- x_full[c(1, 2, 53), ]
  set.seed(15)
  resp3 <- simdat(x = x_mix3, theta = rnorm(1000), D = 1)
  fit3 <- sx2_fit(x_mix3, data = resp3, min.collapse = 0)
  n_mid3 <- sum(rowSums(resp3) > 0 & rowSums(resp3) < 6)
  for (i in 1:3) {
    expect_equal(sum(fit3$obs_freq[[i]]), n_mid3)
    expect_equal(sum(fit3$exp_freq[[i]]), n_mid3)
  }
})

test_that("sx2_fit() keeps observed and expected totals of polytomous items equal", {
  x_poly <- x_full[c(1:10, 53:55), ]
  set.seed(11)
  resp <- simdat(x = x_poly, theta = rnorm(800), D = 1)
  fit <- sx2_fit(x_poly, data = resp)
  n_mid <- sum(rowSums(resp) > 0 & rowSums(resp) < sum(x_poly$cats - 1))
  for (i in 11:13) {
    e <- data.matrix(fit$exp_freq[[i]])
    o <- data.matrix(fit$obs_freq[[i]])
    expect_equal(rowSums(e, na.rm = TRUE), rowSums(o, na.rm = TRUE), ignore_attr = TRUE)
    expect_equal(sum(o, na.rm = TRUE), n_mid)
    # the df counts the remaining cells per score group minus the item parameters
    expect_equal(fit$fit_stat$df[i], sum(rowSums(!is.na(e)) - 1) - 5)
  }
})

test_that("sx2_fit() reports NA critical values and p-values when no degrees of freedom remain", {
  x3 <- x_lsat
  x3$model <- "3PLM"
  x3$par.2 <- c(-3.36, -1.0, 0.2, -1.3, -2.6)
  x3$par.3 <- 0.2
  expect_warning(fit <- sx2_fit(x3, data = LSAT6, D = 1), "No degrees of freedom.*V3")
  expect_equal(fit$fit_stat$df[3], 0)
  expect_true(is.na(fit$fit_stat$p[3]))
  expect_true(is.na(fit$fit_stat$crit.val[3]))

  # the other items keep their critical values and p-values
  ok <- fit$fit_stat$df > 0
  expect_false(anyNA(fit$fit_stat$p[ok]))
  expect_equal(fit$fit_stat$crit.val[ok], round(stats::qchisq(0.95, df = fit$fit_stat$df[ok]), 3))
})

test_that("sx2_fit() names the collapsed tables of a single polytomous item consistently", {
  x_one <- x_full[c(1:20, 55), ]
  set.seed(14)
  resp <- simdat(x = x_one, theta = rnorm(600), D = 1)
  fit <- sx2_fit(x_one, data = resp)
  expect_identical(names(fit$obs_freq[[21]]), paste0("score.", 0:4))
  expect_identical(names(fit$exp_freq[[21]]), paste0("score.", 0:4))

  # the observed and expected tables have the same shape and row names
  expect_identical(dim(fit$obs_freq[[21]]), dim(fit$exp_freq[[21]]))
  expect_identical(rownames(fit$obs_freq[[21]]), rownames(fit$exp_freq[[21]]))
})

test_that("sx2_fit() for an est_irt object uses the latent distribution stored in the object by default", {
  skip_on_cran()
  mod <- est_irt(data = LSAT6, D = 1.702, model = "2PLM", cats = 2, verbose = FALSE)
  base <- function(...) sx2_fit(mod$par.est, data = mod$data, D = mod$scale.D, ...)

  # all three distribution arguments omitted: the stored distribution
  expect_identical(sx2_fit(mod), base(weights = mod$weights))

  # norm.prior or nquad given: the normal distribution, as for the default method
  expect_identical(sx2_fit(mod, norm.prior = c(0, 1), nquad = 30), base())
  expect_identical(sx2_fit(mod, nquad = 21), base(nquad = 21))
  expect_identical(sx2_fit(mod, norm.prior = c(0.2, 1.1)), base(norm.prior = c(0.2, 1.1)))

  # weights given: those weights
  w <- gen.weight(n = 41, dist = "norm", mu = 0.3, sigma = 1.2)
  expect_identical(sx2_fit(mod, weights = w), base(weights = w))

  # the scaling constant comes from the object
  expect_identical(
    sx2_fit(mod)$fit_stat,
    sx2_fit(mod$par.est, data = mod$data, D = 1.702, weights = mod$weights)$fit_stat
  )
})

test_that("sx2_fit() follows the latent distribution estimated with EmpHist = TRUE", {
  skip_on_cran()
  set.seed(22)
  th <- c(rnorm(600, -1, 0.6), rnorm(600, 1, 0.6))
  x_emp <- x_full[1:12, ]
  dat <- simdat(x_emp, th, D = 1)
  mod <- suppressWarnings(
    est_irt(
      data = dat, D = 1, model = "3PLM", cats = 2, use.gprior = TRUE,
      EmpHist = TRUE, Etol = 0.001, verbose = FALSE
    )
  )
  fit_obj <- sx2_fit(mod)
  fit_norm <- sx2_fit(mod, norm.prior = c(0, 1), nquad = 30)
  expect_identical(fit_obj$fit_stat, sx2_fit(mod$par.est, data = dat, D = 1, weights = mod$weights)$fit_stat)
  expect_false(isTRUE(all.equal(fit_obj$fit_stat$chisq, fit_norm$fit_stat$chisq)))

  # an est_item object stores no latent distribution and keeps the normal distribution
  ei <- suppressWarnings(est_item(x = x_emp, data = dat, score = th, D = 1, use.gprior = TRUE, verbose = FALSE))
  expect_identical(sx2_fit(ei)$fit_stat, sx2_fit(ei$par.est, data = dat, D = 1)$fit_stat)
})

test_that("sx2_fit() matches an independent Lord-Wingersky computation without collapsing", {
  fit <- sx2_fit(x_lsat, data = LSAT6, D = 1, min.collapse = 0)

  # quadrature used by sx2_fit() by default
  wts <- gen.weight(n = 30, dist = "norm", mu = 0, sigma = 1)
  pmat <- lapply(seq_len(5), function(i) {
    p1 <- 1 / (1 + exp(-x_lsat$par.1[i] * (wts[, 1] - x_lsat$par.2[i])))
    cbind(1 - p1, p1)
  })

  # summed score distribution at each quadrature point
  lw <- function(p) {
    f <- matrix(1, nrow(wts), 1)
    for (q in p) f <- cbind(f * q[, 1], 0) + cbind(0, f * q[, 2])
    f
  }
  rs <- rowSums(LSAT6)
  n_s <- tabulate(rs + 1, 6)
  denom <- colSums(lw(pmat) * wts[, 2])
  for (i in 1:5) {
    rest <- lw(pmat[-i])
    s <- 1:4

    # expected number correct in each summed score group with the item removed
    e1 <- n_s[s + 1] * colSums(wts[, 2] * pmat[[i]][, 2] * rest[, s, drop = FALSE]) / denom[s + 1]
    o1 <- vapply(s, function(k) sum(rs == k & LSAT6[, i] == 1), numeric(1))
    x2 <- sum((o1 - e1)^2 / e1 + ((n_s[s + 1] - o1) - (n_s[s + 1] - e1))^2 / (n_s[s + 1] - e1))
    expect_equal(fit$fit_stat$chisq[i], round(x2, 3))
    expect_equal(fit$fit_stat$df[i], 4 - 2)
  }
  expect_equal(fit$fit_stat$chisq, c(0.481, 1.984, 3.985, 3.179, 0.320))
  expect_equal(fit$fit_stat$p, c(0.786, 0.371, 0.136, 0.204, 0.852))
})

test_that("sx2_fit() collapses sparse summed score groups and lowers the df", {
  fit <- sx2_fit(x_lsat, data = LSAT6, D = 1)
  expect_equal(fit$fit_stat$df, c(2, 2, 1, 2, 2))
  expect_equal(fit$fit_stat$chisq[3], 3.662)
  expect_equal(fit$fit_stat$p[3], 0.056)

  # every expected cell reaches the minimum after collapsing
  expect_true(all(vapply(fit$exp_freq, function(e) min(data.matrix(e), na.rm = TRUE), numeric(1)) >= 1))
})

test_that("sx2_fit() for est_irt and est_item objects uses the data and scaling constant of the object", {
  skip_on_cran()
  mod <- est_irt(data = LSAT6, D = 1.702, model = "2PLM", cats = 2, verbose = FALSE)
  expect_identical(
    sx2_fit(mod, norm.prior = c(0, 1), nquad = 30)$fit_stat,
    sx2_fit(mod$par.est, data = mod$data, D = 1.702)$fit_stat
  )

  # a different scaling constant gives different statistics
  expect_false(isTRUE(all.equal(
    sx2_fit(mod, norm.prior = c(0, 1), nquad = 30)$fit_stat$chisq,
    sx2_fit(mod$par.est, data = mod$data, D = 1)$fit_stat$chisq
  )))
  ei <- suppressWarnings(
    est_item(x = x_lsat, data = LSAT6, score = rnorm(nrow(LSAT6)), D = 1.702, verbose = FALSE)
  )
  expect_identical(sx2_fit(ei)$fit_stat, sx2_fit(ei$par.est, data = ei$data, D = 1.702)$fit_stat)
})

test_that("sx2_fit() lowers the df by one for each PCM item given in pcm.loc", {
  x_pcm <- x_full[c(1:10, 53:55), ]
  x_pcm[11:13, 3] <- "GPCM"
  x_pcm[11:13, 4] <- 1
  set.seed(13)
  resp <- simdat(x = x_pcm, theta = rnorm(1000), D = 1)
  f1 <- sx2_fit(x_pcm, data = resp)
  f2 <- sx2_fit(x_pcm, data = resp, pcm.loc = 11:13)
  expect_equal(f2$fit_stat$df - f1$fit_stat$df, c(rep(0, 10), rep(1, 3)))
  expect_equal(f2$fit_stat$chisq, f1$fit_stat$chisq)
})

test_that("sx2_fit() stops for a test with a single item", {
  expect_error(
    sx2_fit(x_lsat[1, ], data = LSAT6[, 1, drop = FALSE], D = 1),
    "at least two items"
  )
})

test_that("sx2_fit() reads logical responses like 0 and 1", {
  resp <- as.matrix(LSAT6)
  expect_identical(
    sx2_fit(x_lsat, data = resp == 1, D = 1)$fit_stat,
    sx2_fit(x_lsat, data = resp, D = 1)$fit_stat
  )
})

test_that("the error for responses outside the score categories names the items and columns", {
  resp <- LSAT6
  resp[, 3] <- resp[, 3] + 1
  expect_error(sx2_fit(x_lsat, data = resp, D = 1), "V3 [(]column 3[)]")
})
