# Independent implementation of CATSIB (Nandakumar & Roussos, 2004, Eqs. 7, 9-11): regression
# correction in each group, equal-width bins on the corrected scores with the number of bins
# reduced until at most max.del of each group is left out, and the weighted difference of the
# proportions correct with its standard error
cst_ref <- function(resp, score, se, group, focal, range = c(-5, 5), n.bin = c(80, 10),
                    min.binsize = 3, max.del = 0.075, weight.group = "comb", screen = FALSE) {
  isf <- !is.na(group) & group == focal
  isr <- !is.na(group) & group != focal
  cr <- rep(NA_real_, length(score))
  rho <- numeric(0)
  for (gi in list(isr, isf)) {
    s <- score[gi]
    e <- se[gi]
    ok <- !is.na(s) & s > range[1] & s < range[2]
    if (screen) {
      bad <- ok & !is.na(e) & e^2 >= var(s[ok])
      if (sum(bad) > 0 && sum(ok) - sum(bad) >= 2) ok <- ok & !bad
    }
    m <- mean(s[ok])
    rho2 <- max(0.05, min(1, 1 - mean(e[ok]^2, na.rm = TRUE) / var(s[ok])))
    rho <- c(rho, rho2)
    cr[gi] <- m + rho2 * (s - m)
  }
  lo <- min(cr, na.rm = TRUE)
  hi <- max(cr, na.rm = TRUE)
  out <- NULL
  for (j in seq_len(ncol(resp))) {
    y <- resp[, j]
    okr <- isr & !is.na(y) & !is.na(cr)
    okf <- isf & !is.na(y) & !is.na(cr)
    for (K in n.bin[1]:n.bin[2]) {
      bin <- findInterval(cr, seq(lo, hi, length.out = K + 1), left.open = TRUE)
      bin[bin == 0] <- 1
      nr <- tabulate(bin[okr], K)
      nf <- tabulate(bin[okf], K)
      keep <- nr >= min.binsize & nf >= min.binsize
      if (sum(nr[keep]) / sum(nr) >= 1 - max.del && sum(nf[keep]) / sum(nf) >= 1 - max.del) break
    }
    kk <- which(keep)
    w <- switch(weight.group, comb = (nr + nf)[kk], foc = nf[kk], ref = nr[kk])
    w <- w / sum(w)
    beta <- 0
    vb <- 0
    for (q in seq_along(kk)) {
      yr <- y[okr & bin == kk[q]]
      yf <- y[okf & bin == kk[q]]
      beta <- beta + w[q] * (mean(yr) - mean(yf))
      vb <- vb + w[q]^2 * (var(yr) / length(yr) + var(yf) / length(yf))
    }
    z <- beta / sqrt(vb)
    out <- rbind(out, data.frame(beta = beta, se = sqrt(vb), z.beta = z, p = 2 * pnorm(-abs(z)),
                                 n.ref = sum(nr[kk]), n.foc = sum(nf[kk])))
  }
  attr(out, "rho2") <- rho
  out
}

# Simulate dichotomous two-group data with uniform DIF on items 1 to 3 and impact
cst_sim <- function(seed = 1, J = 20, n = c(700, 600), impact = -0.3, miss = 0) {
  set.seed(seed)
  x <- shape_df(par.drm = list(a = runif(J, 0.8, 1.8), b = rnorm(J), g = 0.15), cats = 2, model = "3PLM")
  xf <- x
  xf$par.2[1:3] <- xf$par.2[1:3] + 0.7
  resp <- rbind(simdat(x, theta = rnorm(n[1]), D = 1), simdat(xf, theta = rnorm(n[2], impact), D = 1))
  if (miss > 0) resp[sample(length(resp), round(miss * length(resp)))] <- NA
  es <- suppressWarnings(est_score(x, resp, D = 1, range = c(-5, 5)))
  list(x = x, resp = resp, group = rep(0:1, n), score = es$est.theta, se = es$se.theta)
}

# Evaluate an expression and return its value together with all warning messages
cst_catch <- function(expr) {
  msgs <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    msgs <<- c(msgs, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = msgs)
}

test_that("catsib() statistics match an independent implementation", {
  sets <- list(
    list(seed = 1, miss = 0, wg = "comb", mb = 3, nb = c(80, 10), skip = NULL),
    list(seed = 2, miss = 0.3, wg = "foc", mb = 5, nb = c(40, 10), skip = c(4, 9)),
    list(seed = 3, miss = 0.5, wg = "ref", mb = 2, nb = c(60, 20), skip = NULL)
  )
  for (st in sets) {
    sim <- cst_sim(st$seed, miss = st$miss)
    for (scr in c(FALSE, TRUE)) {
      rst <- suppressWarnings(catsib(
        x = NULL, data = sim$resp, score = sim$score, se = sim$se, group = sim$group, focal.name = 1,
        weight.group = st$wg, min.binsize = st$mb, n.bin = st$nb, item.skip = st$skip, rho.screen = scr
      ))
      ref <- cst_ref(sim$resp, sim$score, sim$se, sim$group, 1, n.bin = st$nb, min.binsize = st$mb,
                     weight.group = st$wg, screen = scr)
      ref[st$skip, c("beta", "se", "z.beta", "p")] <- NA
      ds <- rst$no_purify$dif_stat
      expect_lt(max(abs(as.matrix(ds[, c("beta", "se", "z.beta", "p")]) -
                          as.matrix(ref[, c("beta", "se", "z.beta", "p")])), na.rm = TRUE), 5.1e-5)
      expect_identical(is.na(ds$p), is.na(ref$p))
      expect_equal(ds$n.ref, ref$n.ref)
      expect_equal(ds$n.foc, ref$n.foc)
      expect_equal(rst$no_purify$dif_item, as.numeric(which(ref$p <= 0.05)))
    }
  }
})

test_that("catsib() excludes examinees with a missing ability estimate with a warning", {
  sim <- cst_sim(5)
  s2 <- sim$score
  s2[c(4, 900)] <- NA
  out <- cst_catch(catsib(NULL, sim$resp, score = s2, se = sim$se, group = sim$group, focal.name = 1))
  r0 <- catsib(NULL, sim$resp[-c(4, 900), ], score = sim$score[-c(4, 900)], se = sim$se[-c(4, 900)],
               group = sim$group[-c(4, 900)], focal.name = 1)
  expect_equal(out$value$no_purify$dif_stat, r0$no_purify$dif_stat)
  expect_identical(out$warnings, "2 examinee(s) with a missing ability estimate were excluded.")
  expect_error(catsib(NULL, sim$resp, score = rep(NA_real_, length(s2)), se = sim$se, group = sim$group,
                      focal.name = 1), "No examinee has an ability estimate")
})

test_that("catsib() applies min.resp also when the scores are supplied", {
  skip_on_cran()
  sim <- cst_sim(6, miss = 0.6)
  few <- rowSums(!is.na(sim$resp)) < 6 & rowSums(!is.na(sim$resp)) > 0
  out <- cst_catch(catsib(sim$x, sim$resp, score = sim$score, se = sim$se, group = sim$group,
                          focal.name = 1, min.resp = 6))
  r0 <- suppressWarnings(catsib(sim$x, sim$resp[!few, ], score = sim$score[!few], se = sim$se[!few],
                                group = sim$group[!few], focal.name = 1))
  expect_equal(out$value$no_purify$dif_stat, r0$no_purify$dif_stat)
  expect_true(paste0(sum(few), " examinee(s) with fewer than 6 responses were excluded.") %in% out$warnings)

  # the scores are estimated without the excluded examinees when they are not supplied
  r1 <- suppressWarnings(catsib(sim$x, sim$resp, group = sim$group, focal.name = 1, min.resp = 6))
  r2 <- suppressWarnings(catsib(sim$x, sim$resp[!few, ], group = sim$group[!few], focal.name = 1))
  expect_equal(r1$no_purify$dif_stat, r2$no_purify$dif_stat)
})

test_that("catsib() purification works down to the last item", {
  skip_on_cran()
  sim <- cst_sim(9)
  # the abilities are estimated from all items, and two items are analyzed; after the first item is
  # removed, the abilities are estimated from the last item
  # with MAP, the corrected abilities of the groups from one item do not overlap, so no bin is
  # usable for the last item and the purification stops with a warning
  out1 <- cst_catch(catsib(sim$x[1:2, ], sim$resp[, 1:2], score = sim$score, se = sim$se,
                           group = sim$group, focal.name = 1, purify = TRUE, alpha = 0.999,
                           method = "MAP", verbose = FALSE))
  expect_true(any(grepl("no bin contained enough examinees", out1$warnings)))

  # the items of the iteration whose results are discarded are not reported as missing
  expect_false(any(grepl("^No bin contained", out1$warnings)))

  # the reliability of the iteration whose results are discarded is not reported
  expect_false(any(grepl("estimated reliability", out1$warnings)))
  expect_false(out1$value$with_purify$complete)
  expect_identical(out1$value$with_purify$n.iter, 0L)
  expect_equal(out1$value$with_purify$dif_item, out1$value$no_purify$dif_item)
  expect_false(anyNA(out1$value$with_purify$dif_stat$n.iter))

  # with ML, the ability estimates from one item are at the limits of range, so the purification
  # stops with a warning and keeps the results of the previous iteration
  out2 <- cst_catch(catsib(sim$x[1:2, ], sim$resp[, 1:2], score = sim$score, se = sim$se,
                           group = sim$group, focal.name = 1, purify = TRUE, alpha = 0.999,
                           verbose = FALSE))
  expect_true(any(grepl("fewer than two examinees of a group", out2$warnings)))
  expect_false(out2$value$with_purify$complete)
  expect_identical(out2$value$with_purify$n.iter, 0L)
  expect_equal(out2$value$with_purify$dif_item, out2$value$no_purify$dif_item)
  expect_equal(out2$value$with_purify$dif_stat[, 1:8], out2$value$no_purify$dif_stat)
})

test_that("catsib() reads character responses as numbers", {
  sim <- cst_sim(11)
  r0 <- catsib(NULL, sim$resp, score = sim$score, se = sim$se, group = sim$group, focal.name = 1)
  chr <- as.data.frame(apply(sim$resp, 2, as.character), stringsAsFactors = FALSE)
  r1 <- catsib(NULL, chr, score = sim$score, se = sim$se, group = sim$group, focal.name = 1)
  expect_identical(r1$no_purify$dif_stat, r0$no_purify$dif_stat)
})

test_that("catsib() ends the purification normally when only skipped items remain", {
  skip_on_cran()
  sim <- cst_sim(3, J = 8, n = c(800, 800))
  out <- cst_catch(catsib(sim$x, sim$resp, group = sim$group, focal.name = 1, item.skip = 2:8,
                          purify = TRUE, verbose = FALSE))
  expect_false(any(grepl("no bin contained enough examinees", out$warnings)))
  expect_true(out$value$with_purify$complete)
  expect_identical(out$value$with_purify$n.iter, 1L)
  expect_equal(out$value$with_purify$dif_item, 1)
})

test_that("catsib() stops with clear errors for invalid inputs", {
  sim <- cst_sim(13)
  n <- nrow(sim$resp)
  run <- function(...) {
    args <- utils::modifyList(
      list(x = sim$x, data = sim$resp, score = sim$score, se = sim$se, group = sim$group,
           focal.name = 1), list(...)
    )
    do.call(catsib, args)
  }
  resp_bad <- sim$resp
  resp_bad[1, 1] <- -1
  expect_error(run(data = resp_bad), "score categories")
  resp_bad[1, 1] <- 0.5
  expect_error(run(data = resp_bad), "score categories")
  resp_bad[1, 1] <- 2
  expect_error(run(data = resp_bad), "dichotomous")
  expect_error(run(data = sim$resp[, -1]), "number of columns")
  expect_error(run(x = NULL, score = NULL), "metadata")
  expect_error(run(score = sim$score[-1]), "'score'")
  expect_error(run(se = sim$se[-1]), "'se'")
  expect_error(run(group = sim$group[-1]), "'group'")
  expect_error(run(focal.name = 3), "focal.name")
  expect_error(run(focal.name = c(0, 1)), "focal.name")
  expect_error(run(group = replace(sim$group, 1:10, 2)), "exactly two groups")
  expect_error(run(item.skip = 30), "item.skip")
  expect_error(run(item.skip = "V1"), "item.skip")
  expect_error(run(alpha = 5), "alpha")
  expect_error(run(purify = TRUE, max.iter = 1.5), "max.iter")
  expect_error(run(min.resp = "a"), "min.resp")
  expect_error(run(n.bin = c(10, 80)), "n.bin")
  expect_error(run(min.binsize = 1), "min.binsize")
  expect_error(run(max.del = 1), "max.del")
  expect_error(run(score = NULL, method = "INV.TCC"), "not supported by catsib()", fixed = TRUE)

  # polytomous items and responses are refused with the name of the function
  x_poly <- sim$x
  x_poly$model[1] <- "GRM"
  x_poly$cats[1] <- 3L
  expect_error(run(x = x_poly), "catsib() supports dichotomous items only", fixed = TRUE)
})

test_that("catsib() leaves skipped items out of the tests and the purification", {
  skip_on_cran()
  sim <- cst_sim(14)
  r0 <- suppressWarnings(catsib(sim$x, sim$resp, score = sim$score, se = sim$se, group = sim$group,
                                focal.name = 1))
  r1 <- suppressWarnings(catsib(sim$x, sim$resp, score = sim$score, se = sim$se, group = sim$group,
                                focal.name = 1, item.skip = c(1, 5)))
  expect_true(all(is.na(r1$no_purify$dif_stat[c(1, 5), 2:5])))
  expect_equal(r1$no_purify$dif_stat[-c(1, 5), ], r0$no_purify$dif_stat[-c(1, 5), ])
  lg <- seq_len(nrow(sim$x)) %in% c(1, 5)
  p1 <- suppressWarnings(catsib(sim$x, sim$resp, score = sim$score, se = sim$se, group = sim$group,
                                focal.name = 1, item.skip = c(1, 5), purify = TRUE, verbose = FALSE))
  p2 <- suppressWarnings(catsib(sim$x, sim$resp, score = sim$score, se = sim$se, group = sim$group,
                                focal.name = 1, item.skip = lg, purify = TRUE, verbose = FALSE))
  expect_identical(p1$with_purify, p2$with_purify)
  expect_false(any(c(1, 5) %in% p1$with_purify$dif_item))
})

test_that("catsib() warns once about the number of cores during the purification", {
  skip_on_cran()
  sim <- cst_sim(15, miss = 0.2)
  out <- cst_catch(catsib(sim$x, sim$resp, group = sim$group, focal.name = 1, purify = TRUE, ncore = 2,
                          verbose = FALSE))
  expect_gte(out$value$with_purify$n.iter, 2)
  expect_equal(sum(grepl("^ncore > 1", out$warnings)), 1)
})

test_that("catsib() flags items with the unrounded p-values", {
  skip_on_cran()
  sim <- cst_sim(7)
  ref <- cst_ref(sim$resp, sim$score, sim$se, sim$group, 1, screen = TRUE)
  j <- which(abs(ref$p - round(ref$p, 4)) > 1e-5 & ref$p > 0.01)[1]
  alpha <- (ref$p[j] + round(ref$p[j], 4)) / 2
  rst <- suppressWarnings(catsib(NULL, sim$resp, score = sim$score, se = sim$se, group = sim$group,
                                 focal.name = 1, alpha = alpha))
  expect_identical(j %in% rst$no_purify$dif_item, ref$p[j] <= alpha)
})

test_that("catsib() purification removes the flagged item with the smallest p-value first", {
  skip_on_cran()
  sim <- cst_sim(8)
  # the reliability of the scores falls below the lower bound as many items are removed, and one
  # warning reports it
  out <- cst_catch(catsib(sim$x, sim$resp, score = sim$score, se = sim$se, group = sim$group,
                          focal.name = 1, purify = TRUE, max.iter = 20, verbose = FALSE))
  rst <- out$value
  expect_gte(rst$with_purify$n.iter, 2)
  expect_equal(sum(grepl("estimated reliability", out$warnings)), 1)
  first <- which(rst$with_purify$dif_stat$n.iter == 0)
  z0 <- rst$no_purify$dif_stat$z.beta
  expect_equal(first, which.max(abs(z0)))
  expect_true(all(rst$with_purify$dif_item %in% seq_len(nrow(sim$x))))
  expect_type(rst$with_purify$n.iter, "integer")

  # the statistics of the remaining items are those computed with the scores from these items
  keep <- setdiff(seq_len(nrow(sim$x)), rst$with_purify$dif_item)
  if (rst$with_purify$complete) {
    es <- est_score(sim$x[keep, ], sim$resp[, keep], D = 1, range = c(-5, 5))
    ref <- cst_ref(sim$resp[, keep], es$est.theta, es$se.theta, sim$group, 1, screen = TRUE)
    expect_lt(max(abs(rst$with_purify$dif_stat$beta[keep] - ref$beta)), 5.1e-5)
  }
})

test_that("catsib() reports an incomplete purification when max.iter is reached", {
  skip_on_cran()
  sim <- cst_sim(10)
  out <- cst_catch(catsib(sim$x, sim$resp, score = sim$score, se = sim$se, group = sim$group,
                          focal.name = 1, purify = TRUE, max.iter = 1, alpha = 0.3, verbose = FALSE))
  expect_true(any(grepl("maximum number of iterations", out$warnings)))
  expect_false(out$value$with_purify$complete)
  expect_true(length(out$value$with_purify$dif_item) > 1)
})

test_that("catsib() leaves the estimates with large standard errors out of rho^2 only when rho.screen = TRUE", {
  skip_on_cran()
  sim <- cst_sim(4)
  ref0 <- cst_ref(sim$resp, sim$score, sim$se, sim$group, 1)
  ref1 <- cst_ref(sim$resp, sim$score, sim$se, sim$group, 1, screen = TRUE)
  expect_true(any(attr(ref1, "rho2") > attr(ref0, "rho2")))
  r0 <- catsib(NULL, sim$resp, score = sim$score, se = sim$se, group = sim$group, focal.name = 1,
               rho.screen = FALSE)
  r1 <- catsib(NULL, sim$resp, score = sim$score, se = sim$se, group = sim$group, focal.name = 1)
  expect_lt(max(abs(r0$no_purify$dif_stat$beta - ref0$beta)), 5.1e-5)
  expect_lt(max(abs(r1$no_purify$dif_stat$beta - ref1$beta)), 5.1e-5)

  # the screened examinees are still matched and counted
  expect_equal(r1$no_purify$dif_stat$n.total, ref1$n.ref + ref1$n.foc)

  # the screening is not applied when fewer than two estimates would be left
  use <- c(TRUE, TRUE, TRUE, NA)
  expect_identical(rho_screen(score = c(-1, 0, 1, NA), se = c(0.1, 5, 5, 1), use = use), use)
  expect_identical(rho_screen(score = c(-1, 0, 1, 2), se = c(0.1, 0.2, 5, 1), use = rep(TRUE, 4)),
                   c(TRUE, TRUE, FALSE, TRUE))
  expect_error(catsib(NULL, sim$resp, score = sim$score, se = sim$se, group = sim$group, focal.name = 1,
                      rho.screen = NA), "rho.screen")
})

test_that("catsib() gives missing statistics when the corrected abilities of the groups do not overlap", {
  set.seed(21)
  n <- 300
  group <- rep(0:1, each = n)
  score <- c(rnorm(n, -2, 0.5), rnorm(n, 2, 0.5))
  se <- rep(0.5, 2 * n)
  resp <- matrix(rbinom(2 * n * 4, 1, 0.5), ncol = 4)
  out <- cst_catch(catsib(NULL, resp, score = score, se = se, group = group, focal.name = 1,
                          rho.screen = FALSE))
  ds <- out$value$no_purify$dif_stat
  expect_true(all(is.na(ds[, c("beta", "se", "z.beta", "p")])))
  expect_true(all(ds$n.total == 0))
  expect_null(out$value$no_purify$dif_item)
  expect_true(any(grepl("No bin contained enough examinees of both groups", out$warnings)))
})

test_that("print() of catsib results leaves the significance symbols of skipped items blank", {
  sim <- cst_sim(12)
  rst <- suppressWarnings(catsib(sim$x, sim$resp, score = sim$score, se = sim$se, group = sim$group,
                                 focal.name = 1, item.skip = 2))
  out <- utils::capture.output(print(rst))
  expect_false(any(grepl("?", out, fixed = TRUE)))
})

test_that("catsib() names in the warning only the items with missing statistics in the reported tables", {
  skip_on_cran()
  for (seed in c(43, 45, 46)) {
    sim <- cst_sim(seed, J = 14, n = c(200, 160), miss = 0.2)
    out <- cst_catch(catsib(sim$x, sim$resp, score = sim$score, se = sim$se, group = sim$group,
                            focal.name = 1, purify = TRUE, max.iter = 15, alpha = 0.3,
                            min.binsize = 25, n.bin = c(30, 5), verbose = FALSE))
    rst <- out$value
    na_item <- union(which(is.na(rst$no_purify$dif_stat$beta)), which(is.na(rst$with_purify$dif_stat$beta)))
    warned <- out$warnings[grepl("^No bin contained", out$warnings)]
    listed <- if (length(warned) == 0L) {
      character()
    } else {
      strsplit(sub("^.*item[(]s[)] (.*), so their statistics.*$", "\\1", warned), ", ")[[1]]
    }
    expect_identical(sort(listed), sort(rst$no_purify$dif_stat$id[na_item]))
  }

  # an item without a usable bin in an iteration of the purification only is not reported when its
  # statistics are computed in the reported tables
  sim <- cst_sim(45, J = 14, n = c(200, 160), miss = 0.2)
  out <- cst_catch(catsib(sim$x, sim$resp, score = sim$score, se = sim$se, group = sim$group,
                          focal.name = 1, purify = TRUE, max.iter = 15, alpha = 0.3,
                          min.binsize = 25, n.bin = c(30, 5), verbose = FALSE))
  expect_false(anyNA(out$value$with_purify$dif_stat$beta))
  expect_false(any(grepl("^No bin contained", out$warnings)))
})

test_that("catsib() does not name skipped items in the warning about missing statistics", {
  set.seed(21)
  n <- 300
  group <- rep(0:1, each = n)
  score <- c(rnorm(n, -2, 0.5), rnorm(n, 2, 0.5))
  se <- rep(0.5, 2 * n)
  resp <- matrix(rbinom(2 * n * 4, 1, 0.5), ncol = 4)
  out <- cst_catch(catsib(NULL, resp, score = score, se = se, group = group, focal.name = 1,
                          rho.screen = FALSE, item.skip = 2))
  warned <- out$warnings[grepl("^No bin contained", out$warnings)]
  expect_length(warned, 1L)
  expect_match(warned, "item(s) item.1, item.3, item.4,", fixed = TRUE)
  expect_false(grepl("item.2", warned, fixed = TRUE))
})
