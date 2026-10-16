# Independent implementation of the RDIF statistics (Lim, Choe, & Han, 2022,
# Eqs. 7-20), extended to polytomous items with the expected item score residual
rdt_probs <- function(x, theta, D) {
  lapply(seq_len(nrow(x)), function(j) {
    K <- x$cats[j]
    a <- x$par.1[j]
    if (K == 2) {
      g <- if (x$model[j] == "3PLM") x$par.3[j] else 0
      P <- g + (1 - g) / (1 + exp(-D * a * (theta - x$par.2[j])))
      return(cbind(1 - P, P))
    }
    d <- unlist(x[j, paste0("par.", 2:K)])
    if (x$model[j] == "GRM") {
      ps <- cbind(1, sapply(d, function(b) 1 / (1 + exp(-D * a * (theta - b)))), 0)
      return(ps[, 1:K] - ps[, 2:(K + 1)])
    }
    z <- cbind(0, t(apply(sapply(d, function(b) D * a * (theta - b)), 1, cumsum)))
    e <- exp(z - apply(z, 1, max))
    e / rowSums(e)
  })
}

rdt_ref <- function(x, resp, score, group, focal, D) {
  P <- rdt_probs(x, ifelse(is.na(score), 0, score), D)
  out <- NULL
  for (j in seq_len(nrow(x))) {
    k <- 0:(x$cats[j] - 1)
    ey <- drop(P[[j]] %*% k)
    dev <- outer(-ey, k, "+")
    m2 <- rowSums(dev^2 * P[[j]])
    m3 <- rowSums(dev^3 * P[[j]])
    m4 <- rowSums(dev^4 * P[[j]])
    ok <- !is.na(resp[, j]) & !is.na(score)
    f <- ok & group == focal
    r <- ok & group != focal
    e <- resp[, j] - ey
    nf <- sum(f)
    nr <- sum(r)
    vr <- mean(e[f]) - mean(e[r])
    vs <- mean(e[f]^2) - mean(e[r]^2)
    mus <- mean(m2[f]) - mean(m2[r])
    s2r <- sum(m2[f]) / nf^2 + sum(m2[r]) / nr^2
    s2s <- sum(m4[f] - m2[f]^2) / nf^2 + sum(m4[r] - m2[r]^2) / nr^2
    crs <- sum(m3[f]) / nf^2 + sum(m3[r]) / nr^2
    dv <- c(vr, vs - mus)
    chi <- drop(dv %*% solve(matrix(c(s2r, crs, crs, s2s), 2), dv))
    zr <- vr / sqrt(s2r)
    zs <- (vs - mus) / sqrt(s2s)
    out <- rbind(out, data.frame(
      rdifr = vr, z.rdifr = zr, rdifs = vs, z.rdifs = zs, rdifrs = chi,
      p.rdifr = 2 * pnorm(-abs(zr)), p.rdifs = 2 * pnorm(-abs(zs)),
      p.rdifrs = pchisq(chi, 2, lower.tail = FALSE), n.ref = nr, n.foc = nf,
      mu.rdifs = mus, sigma.rdifr = sqrt(s2r), sigma.rdifs = sqrt(s2s), covariance = crs
    ))
  }
  out
}

# Simulate a mixed-format two-group data set with DIF on items 1, 2, and 9
rdt_sim <- function(seed = 21, miss = 0.25) {
  set.seed(seed)
  x <- shape_df(
    par.drm = list(a = runif(8, 0.8, 2), b = rnorm(8), g = c(rep(NA, 4), runif(4, 0.05, 0.25))),
    par.prm = list(a = c(1.2, 0.9, 1.5), d = list(c(-1, 0, 1), c(-0.8, 0.2, 1.1, 1.6), c(-0.5, 0.7))),
    item.id = paste0("I", 1:11), cats = c(rep(2, 8), 4, 5, 3),
    model = c(rep("2PLM", 4), rep("3PLM", 4), "GRM", "GPCM", "GRM")
  )
  xf <- x
  xf$par.2[1:2] <- xf$par.2[1:2] + 0.6
  xf[9, c("par.2", "par.3", "par.4")] <- xf[9, c("par.2", "par.3", "par.4")] + 0.6
  resp <- rbind(simdat(x, theta = rnorm(600), D = 1), simdat(xf, theta = rnorm(450, -0.2), D = 1))
  resp[sample(length(resp), miss * length(resp))] <- NA
  keep <- rowSums(!is.na(resp)) > 0
  group <- c(rep("r", 600), rep("f", 450))
  list(x = x, resp = resp[keep, ], group = group[keep])
}

# Simulate a mixed-format data set in which the five-category item shows large DIF
rdt_sim_poly <- function(seed = 11) {
  set.seed(seed)
  x <- shape_df(
    par.drm = list(a = c(1, 1.3, 0.9, 1.6), b = c(-0.5, 0, 0.4, 0.8), g = rep(NA, 4)),
    par.prm = list(a = c(1.2, 1.4), d = list(c(-0.5, 0.7), c(-1.2, -0.3, 0.5, 1.4))),
    item.id = paste0("I", 1:6), cats = c(rep(2, 4), 3, 5), model = c(rep("2PLM", 4), "GRM", "GPCM")
  )
  xf <- x
  xf[6, c("par.2", "par.3", "par.4", "par.5")] <- xf[6, c("par.2", "par.3", "par.4", "par.5")] + 1.2
  resp <- rbind(simdat(x, theta = rnorm(500), D = 1), simdat(xf, theta = rnorm(500), D = 1))
  list(x = x, resp = resp, group = rep(c("r", "f"), each = 500))
}

# Simulate dichotomous two-group data in which only item 1 is answered by examinees with a common
# supplied ability, so that the covariance matrix of item 1 is singular and its RDIF_RS has one
# degree of freedom
rdt_sing_sim <- function(seed = 1) {
  set.seed(seed)
  nitem <- 6
  model <- sample(c("1PLM", "2PLM", "3PLM"), nitem, replace = TRUE)
  a <- ifelse(model == "1PLM", 1, runif(nitem, 0.7, 2))
  b <- rnorm(nitem)
  g <- ifelse(model == "3PLM", runif(nitem, 0.05, 0.25), 0)
  x <- data.frame(id = paste0("I", seq_len(nitem)), cats = 2, model = model, par.1 = a, par.2 = b,
                  par.3 = ifelse(model == "3PLM", g, NA), stringsAsFactors = FALSE)
  drift <- sort(sample(2:nitem, 2))
  bf <- b
  bf[drift] <- bf[drift] + sample(c(-1, 1), 2, TRUE) * 0.3
  af <- a
  af[drift[1]] <- af[drift[1]] * 0.5
  thr <- rnorm(300)
  thf <- rnorm(300, 0.2)
  sim <- function(th, a, b) {
    p <- sapply(seq_len(nitem), function(j) g[j] + (1 - g[j]) / (1 + exp(-a[j] * (th - b[j]))))
    matrix(rbinom(length(p), 1, p), nrow = length(th))
  }
  resp <- rbind(sim(thr, a, b), sim(thf, af, bf))
  resp[sample(length(resp), round(0.3 * length(resp)))] <- NA
  group <- rep(0:1, each = 300)
  score <- suppressWarnings(est_score(x, resp, D = 1)$est.theta)

  # only 200 examinees answer item 1, and their supplied abilities are all equal
  resp[, 1] <- NA
  idx <- sample(nrow(resp), 200)
  resp[idx, 1] <- rbinom(200, 1, ifelse(group[idx] == 1, 0.62, 0.5))
  score[idx] <- 0.7
  list(x = x, resp = resp, score = score, group = group)
}

# Evaluate an expression and return its value together with all warning messages
rdt_catch <- function(expr) {
  msgs <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    msgs <<- c(msgs, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = msgs)
}

test_that("rdif() methods for est_irt and est_item objects match the default method", {
  skip_on_cran()
  sim <- rdt_sim(miss = 0)
  x2 <- sim$x[1:8, ]
  resp <- sim$resp[, 1:8]
  fit <- est_irt(data = resp, D = 1, model = "2PLM", cats = 2, verbose = FALSE)
  r1 <- rdif(fit, group = sim$group, focal.name = "f", purify = TRUE, verbose = FALSE)
  r2 <- rdif(fit$par.est, data = fit$data, group = sim$group, focal.name = "f", D = fit$scale.D,
             purify = TRUE, verbose = FALSE)
  expect_identical(r1[c("no_purify", "with_purify")], r2[c("no_purify", "with_purify")])
  score <- est_score(x2, resp, D = 1)$est.theta
  fit2 <- est_item(x = x2, data = resp, score = score, D = 1, verbose = FALSE)
  r3 <- rdif(fit2, group = sim$group, focal.name = "f")
  r4 <- rdif(fit2$par.est, data = fit2$data, score = fit2$score, group = sim$group,
             focal.name = "f", D = fit2$scale.D)
  expect_identical(r3$no_purify, r4$no_purify)

  # the data and the scaling factor are taken from the object
  expect_error(rdif(fit, data = resp, group = sim$group, focal.name = "f"), "taken from the object")
  expect_error(rdif(fit2, D = 1.7, group = sim$group, focal.name = "f"), "taken from the object")

  # the ability estimates are also taken from an est_item object
  expect_error(rdif(fit2, score = score, group = sim$group, focal.name = "f"), "taken from the object")
})

test_that("rdif() warns once about min.resp when the scores are estimated", {
  skip_on_cran()
  sim <- rdt_sim()
  out <- rdt_catch(rdif(sim$x, sim$resp, group = sim$group, focal.name = "f", min.resp = 8))
  expect_length(out$warnings, 1)
  expect_match(out$warnings, "fewer than 8 responses")
})

test_that("rdif() warns once about the number of cores during the purification", {
  skip_on_cran()
  sim <- rdt_sim()
  out <- rdt_catch(rdif(sim$x, sim$resp, group = sim$group, focal.name = "f", purify = TRUE,
                        ncore = 2, verbose = FALSE))
  expect_gte(out$value$with_purify$n.iter, 2)
  expect_equal(sum(grepl("^ncore > 1", out$warnings)), 1)
})

test_that("rdif() leaves skipped items out of the tests and the purification", {
  skip_on_cran()
  sim <- rdt_sim()
  score <- suppressWarnings(est_score(sim$x, sim$resp, D = 1)$est.theta)
  r0 <- rdif(sim$x, sim$resp, score = score, group = sim$group, focal.name = "f")
  r1 <- rdif(sim$x, sim$resp, score = score, group = sim$group, focal.name = "f", item.skip = c(2, 5))
  expect_true(all(is.na(r1$no_purify$dif_stat[c(2, 5), 2:9])))
  expect_equal(r1$no_purify$dif_stat$n.total, r0$no_purify$dif_stat$n.total)
  expect_equal(r1$no_purify$dif_stat[-c(2, 5), ], r0$no_purify$dif_stat[-c(2, 5), ])
  expect_false(any(c(2, 5) %in% unlist(r1$no_purify$dif_item)))

  # a logical item.skip gives the same purification as the item positions
  lg <- seq_len(nrow(sim$x)) %in% c(2, 5)
  p1 <- rdif(sim$x, sim$resp, score = score, group = sim$group, focal.name = "f", item.skip = c(2, 5),
             purify = TRUE, verbose = FALSE)
  p2 <- rdif(sim$x, sim$resp, score = score, group = sim$group, focal.name = "f", item.skip = lg,
             purify = TRUE, verbose = FALSE)
  expect_identical(p1$with_purify, p2$with_purify)
  expect_true(all(is.na(p1$with_purify$dif_stat$rdifrs[c(2, 5)])))
})

test_that("rdif() excludes examinees with a missing ability estimate with a warning", {
  skip_on_cran()
  sim <- rdt_sim()
  score <- suppressWarnings(est_score(sim$x, sim$resp, D = 1)$est.theta)
  score[c(1:30, 701:720)] <- NA
  expect_warning(
    rst <- rdif(sim$x, sim$resp, score = score, group = sim$group, focal.name = "f"),
    "50 examinee\\(s\\) with a missing ability estimate"
  )
  ref <- rdt_ref(sim$x, sim$resp, score, sim$group, "f", 1)
  expect_lt(max(abs(rst$no_purify$dif_stat$rdifrs - ref$rdifrs)), 1e-3)
  expect_equal(rst$no_purify$dif_stat$n.ref, ref$n.ref)
  expect_equal(rst$no_purify$dif_stat$n.foc, ref$n.foc)
  expect_error(
    rdif(sim$x, sim$resp, score = rep(NA_real_, nrow(sim$resp)), group = sim$group, focal.name = "f"),
    "No examinee has an ability estimate"
  )

  # an examinee without any response is scored as NA by est_score()
  resp <- sim$resp
  resp[3, ] <- NA
  rst2 <- suppressWarnings(rdif(sim$x, resp, group = sim$group, focal.name = "f"))
  expect_true(is.na(rst2$no_purify$score[3]))
  expect_false(anyNA(rst2$no_purify$dif_stat$rdifrs))

  # the excluded examinees stay excluded during purification
  out <- rdt_catch(rdif(sim$x, sim$resp, score = score, group = sim$group, focal.name = "f",
                        purify = TRUE, verbose = FALSE))
  expect_true(any(grepl("50 examinee\\(s\\) with a missing ability estimate", out$warnings)))
  expect_false(anyNA(out$value$with_purify$dif_stat$rdifrs))
})

test_that("rdif() reports the examinees who lose their ability estimates in the purification once", {
  skip_on_cran()
  sim <- rdt_sim()
  score <- suppressWarnings(est_score(sim$x, sim$resp, D = 1)$est.theta)

  # five examinees respond only to the item that is flagged first and removed in the purification
  rst0 <- rdif(sim$x, sim$resp, score = score, group = sim$group, focal.name = "f")
  first <- which.max(rst0$no_purify$dif_stat$rdifrs)
  resp <- sim$resp
  resp[5:9, ] <- NA
  resp[5:9, first] <- 1
  out <- rdt_catch(rdif(sim$x, resp, score = score, group = sim$group, focal.name = "f",
                        purify = TRUE, verbose = FALSE))
  expect_true(first %in% out$value$with_purify$dif_item)
  msg <- grep("examinee\\(s\\) were excluded during the purification", out$warnings, value = TRUE)
  expect_length(msg, 1)
  expect_gte(as.numeric(sub(" .*", "", msg)), 5)
  expect_false(any(grepl("NA values are returned", out$warnings)))
})

test_that("rdif() statistics match an independent implementation", {
  sim <- rdt_sim()
  for (D in c(1, 1.702)) {
    score <- suppressWarnings(est_score(sim$x, sim$resp, D = D, method = "ML")$est.theta)
    rst <- rdif(sim$x, sim$resp, score = score, group = sim$group, focal.name = "f", D = D)
    ref <- rdt_ref(sim$x, sim$resp, score, sim$group, "f", D)
    stat <- rst$no_purify$dif_stat
    for (v in c("rdifr", "z.rdifr", "rdifs", "z.rdifs", "rdifrs", "p.rdifr", "p.rdifs", "p.rdifrs")) {
      expect_lt(max(abs(stat[[v]] - ref[[v]])), 1e-4)
    }
    expect_equal(stat$n.ref, ref$n.ref)
    expect_equal(stat$n.foc, ref$n.foc)
    mmt <- rst$no_purify$moments
    expect_equal(mmt$mu.rdifr, rep(0, nrow(sim$x)))
    expect_equal(mmt$mu.rdifs, ref$mu.rdifs, tolerance = 1e-8)
    expect_equal(mmt$sigma.rdifr, ref$sigma.rdifr, tolerance = 1e-8)
    expect_equal(mmt$sigma.rdifs, ref$sigma.rdifs, tolerance = 1e-8)
    expect_equal(mmt$covariance, ref$covariance, tolerance = 1e-8)
    expect_equal(rst$no_purify$dif_item$rdifr, which(ref$p.rdifr <= 0.05))
    expect_equal(rst$no_purify$dif_item$rdifrs, which(ref$p.rdifrs <= 0.05))
  }
})

test_that("rdif() estimates the scores with est_score() when score = NULL", {
  skip_on_cran()
  sim <- rdt_sim()
  score <- suppressWarnings(est_score(sim$x, sim$resp, D = 1, method = "MAP")$est.theta)
  r1 <- rdif(sim$x, sim$resp, score = score, group = sim$group, focal.name = "f")
  r2 <- rdif(sim$x, sim$resp, group = sim$group, focal.name = "f", method = "MAP")
  expect_identical(r1$no_purify$dif_stat, r2$no_purify$dif_stat)
  expect_identical(r2$no_purify$score, score)
})

test_that("rdif() purification removes the item with the smallest p-value first", {
  skip_on_cran()
  sim <- rdt_sim()
  rst <- suppressWarnings(rdif(sim$x, sim$resp, group = sim$group, focal.name = "f",
                               purify = TRUE, verbose = FALSE))
  ref <- rdt_ref(sim$x, sim$resp, rst$no_purify$score, sim$group, "f", 1)
  first <- which.min(ref$p.rdifrs)
  expect_true(first %in% rst$with_purify$dif_item)
  expect_equal(unlist(rst$with_purify$dif_stat[first, 1:12]), unlist(rst$no_purify$dif_stat[first, ]))
  expect_equal(rst$with_purify$dif_stat$n.iter[first], 0)
  expect_true(rst$with_purify$complete)
  expect_true(all(rst$with_purify$dif_stat$n.iter <= rst$with_purify$n.iter))

  # the final scores are the estimates without the flagged items
  keep <- setdiff(seq_len(nrow(sim$x)), rst$with_purify$dif_item)
  sc <- suppressWarnings(est_score(sim$x[keep, ], sim$resp[, keep], D = 1)$est.theta)
  expect_equal(rst$with_purify$score, sc)
})

test_that("rdif() reports an incomplete purification when max.iter is reached", {
  skip_on_cran()
  sim <- rdt_sim()
  expect_warning(
    rst <- rdif(sim$x, sim$resp, group = sim$group, focal.name = "f", purify = TRUE,
                max.iter = 1, alpha = 0.3, verbose = FALSE),
    "maximum number of iterations"
  )
  expect_false(rst$with_purify$complete)
  expect_equal(rst$with_purify$n.iter, 1)
  last <- which(rst$with_purify$dif_stat$n.iter == 1 & rst$with_purify$dif_stat$p.rdifrs <= 0.3)
  expect_true(all(last %in% rst$with_purify$dif_item))
})

test_that("rdif() applies min.resp when the scores are estimated", {
  skip_on_cran()
  sim <- rdt_sim()
  n_resp <- rowSums(!is.na(sim$resp))
  drop_row <- n_resp < 8
  expect_true(any(drop_row))
  resp <- sim$resp
  resp[drop_row, ] <- NA
  out <- rdt_catch(rdif(sim$x, sim$resp, group = sim$group, focal.name = "f", min.resp = 8))
  r1 <- out$value
  expect_true(any(grepl("fewer than 8 responses", out$warnings)))
  r2 <- suppressWarnings(rdif(sim$x, resp, group = sim$group, focal.name = "f"))
  expect_identical(r1$no_purify$dif_stat, r2$no_purify$dif_stat)
  expect_equal(r1$no_purify$dif_stat$n.total, unname(colSums(!is.na(resp))))
})

test_that("rdif() applies min.resp also when the scores are supplied", {
  skip_on_cran()
  sim <- rdt_sim()
  n_resp <- rowSums(!is.na(sim$resp))
  drop_row <- n_resp < 8
  resp <- sim$resp
  resp[drop_row, ] <- NA
  score <- suppressWarnings(est_score(sim$x, sim$resp, D = 1)$est.theta)
  score_na <- score
  score_na[drop_row] <- NA
  out <- rdt_catch(rdif(sim$x, sim$resp, score = score, group = sim$group, focal.name = "f", min.resp = 8))
  r3 <- out$value
  expect_true(any(grepl("fewer than 8 responses", out$warnings)))
  r4 <- rdif(sim$x, resp, score = score_na, group = sim$group, focal.name = "f")
  expect_identical(r3$no_purify$dif_stat, r4$no_purify$dif_stat)
  expect_equal(r3$no_purify$dif_stat$n.total, unname(colSums(!is.na(resp))))

  # the examinees stay excluded during purification
  keep4 <- n_resp >= 4 | n_resp == 0
  out2 <- rdt_catch(rdif(sim$x, sim$resp, score = score, group = sim$group, focal.name = "f",
                         min.resp = 4, purify = TRUE, verbose = FALSE))
  expect_false(anyNA(out2$value$with_purify$dif_stat$rdifrs))
  expect_true(all(out2$value$with_purify$dif_stat$n.total <= sum(keep4 & n_resp > 0)))
})

test_that("rdif() purification works when items are removed down to the last one", {
  skip_on_cran()
  x1 <- shape_df(par.drm = list(a = c(1, 1.2, 0.9), b = c(0, 0.5, -0.5), g = 0), cats = 2, model = "2PLM")
  x2 <- x1
  x2$par.2 <- x2$par.2 + c(1.5, 1.2, 0)
  set.seed(5)
  resp <- rbind(simdat(x1, theta = rnorm(500), D = 1), simdat(x2, theta = rnorm(500), D = 1))
  group <- rep(0:1, each = 500)
  rst <- suppressWarnings(rdif(x1, resp, group = group, focal.name = 1, purify = TRUE, verbose = FALSE))
  expect_true(length(rst$with_purify$dif_item) >= 2)
  expect_false(anyNA(rst$with_purify$dif_stat$rdifrs))

  # the purification stops with a warning when every item is flagged
  out <- rdt_catch(rdif(x1[1:2, ], resp[, 1:2], group = group, focal.name = 1, purify = TRUE,
                        alpha = 0.999, verbose = FALSE))
  expect_true(any(grepl("All items were flagged", out$warnings)))
  expect_false(out$value$with_purify$complete)
  expect_equal(out$value$with_purify$dif_item, 1:2)
  expect_equal(out$value$with_purify$n.iter, 1)
  expect_false(anyNA(out$value$with_purify$dif_stat$rdifrs))
  rst0 <- rdif(x1[1:2, ], resp[, 1:2], group = group, focal.name = 1, alpha = 0.999)
  expect_identical(out$value$no_purify, rst0$no_purify)
})

test_that("rdif() purification stops with a warning when a group has no ability estimate left", {
  skip_on_cran()
  sim <- rdt_sim()
  score <- suppressWarnings(est_score(sim$x, sim$resp, D = 1)$est.theta)

  # min.resp counts the remaining items, so the examinees drop out as the items are removed
  out <- rdt_catch(rdif(sim$x, sim$resp, score = score, group = sim$group, focal.name = "f",
                        min.resp = 10, purify = TRUE, alpha = 0.5, verbose = FALSE))
  expect_true(any(grepl("no examinee of a group", out$warnings)))
  expect_false(out$value$with_purify$complete)

  # the items that remain keep the results of the previous iteration, and the flagged ones are reported
  expect_false(anyNA(out$value$with_purify$dif_stat$rdifrs))
  expect_false(anyNA(out$value$with_purify$dif_stat$n.total))
  expect_false(anyNA(out$value$with_purify$moments$mu.rdifr))
  expect_true(all(out$value$with_purify$dif_item %in% seq_len(nrow(sim$x))))
  expect_lte(max(out$value$with_purify$dif_stat$n.iter), out$value$with_purify$n.iter)
  rst0 <- suppressWarnings(rdif(sim$x, sim$resp, score = score, group = sim$group, focal.name = "f",
                                min.resp = 10, alpha = 0.5))
  expect_identical(out$value$no_purify, rst0$no_purify)

  # the returned scores and the examinees reported as excluded are those of the last iteration counted
  out2 <- rdt_catch(rdif(sim$x, sim$resp, score = score, group = sim$group, focal.name = "f",
                         min.resp = 10, purify = TRUE, alpha = 0.5, verbose = FALSE,
                         max.iter = out$value$with_purify$n.iter))
  expect_identical(out$value$with_purify$dif_stat, out2$value$with_purify$dif_stat)
  expect_identical(out$value$with_purify$score, out2$value$with_purify$score)
  expect_false(all(is.na(out$value$with_purify$score[sim$group == "f"])))
  expect_identical(grep("excluded during the purification", out$warnings, value = TRUE),
                   grep("excluded during the purification", out2$warnings, value = TRUE))
})

test_that("rdif() stops for a scoring method without pattern scores when it estimates abilities", {
  skip_on_cran()
  sim <- rdt_sim()
  score <- suppressWarnings(est_score(sim$x, sim$resp, D = 1)$est.theta)
  for (m in c("EAP.SUM", "INV.TCC")) {
    expect_error(rdif(sim$x, sim$resp, group = sim$group, focal.name = "f", method = m),
                 "not supported by rdif()", fixed = TRUE)
    expect_error(rdif(sim$x, sim$resp, score = score, group = sim$group, focal.name = "f",
                      method = m, purify = TRUE, verbose = FALSE),
                 "not supported by rdif()", fixed = TRUE)

    # the method is not used when the scores are supplied and no purification is carried out
    expect_no_error(rdif(sim$x, sim$resp, score = score, group = sim$group, focal.name = "f",
                         method = m))
  }
})

test_that("rdif() purification works when the item with the most categories is removed", {
  skip_on_cran()
  sim <- rdt_sim_poly()
  rst <- rdif(sim$x, sim$resp, group = sim$group, focal.name = "f", purify = TRUE, verbose = FALSE)
  expect_true(6 %in% rst$with_purify$dif_item)
  expect_equal(rst$with_purify$dif_stat$n.iter[6], 0)
  expect_false(anyNA(rst$with_purify$dif_stat$rdifrs))
})

test_that("rdif() reads character responses as numbers", {
  skip_on_cran()
  sim <- rdt_sim()
  score <- suppressWarnings(est_score(sim$x, sim$resp, D = 1)$est.theta)
  resp_chr <- sim$resp
  storage.mode(resp_chr) <- "character"
  r1 <- rdif(sim$x, sim$resp, score = score, group = sim$group, focal.name = "f")
  r2 <- rdif(sim$x, resp_chr, score = score, group = sim$group, focal.name = "f")
  r3 <- rdif(sim$x, as.data.frame(resp_chr), score = score, group = sim$group, focal.name = "f")
  expect_identical(r1$no_purify$dif_stat, r2$no_purify$dif_stat)
  expect_identical(r1$no_purify$dif_stat, r3$no_purify$dif_stat)
})

test_that("rdif() flags items with the unrounded p-values", {
  skip_on_cran()
  sim <- rdt_sim()
  score <- suppressWarnings(est_score(sim$x, sim$resp, D = 1)$est.theta)
  ref <- rdt_ref(sim$x, sim$resp, score, sim$group, "f", 1)

  # an item whose rounded p-value exceeds its unrounded p-value
  j <- which(round(ref$p.rdifr, 4) > ref$p.rdifr & ref$p.rdifr < 0.9)[1]
  expect_false(is.na(j))
  a <- ref$p.rdifr[j] * (1 + 1e-9)
  rst <- rdif(sim$x, sim$resp, score = score, group = sim$group, focal.name = "f", D = 1, alpha = a)
  expect_true(j %in% rst$no_purify$dif_item$rdifr)
  expect_equal(rst$no_purify$dif_stat$p.rdifr[j], round(ref$p.rdifr[j], 4))

  # the same holds for RDIF_S and RDIF_RS
  j2 <- which(round(ref$p.rdifrs, 4) > ref$p.rdifrs & ref$p.rdifrs < 0.9)[1]
  expect_false(is.na(j2))
  rst2 <- rdif(sim$x, sim$resp, score = score, group = sim$group, focal.name = "f", D = 1,
               alpha = ref$p.rdifrs[j2] * (1 + 1e-9))
  expect_true(j2 %in% rst2$no_purify$dif_item$rdifrs)
})

test_that("rdif() purification removes only flagged items when an item has a singular covariance", {
  sim <- rdt_sing_sim()
  st <- suppressWarnings(rdif(sim$x, sim$resp, score = sim$score, group = sim$group,
                           focal.name = 1))$no_purify$dif_stat
  p <- st$p.rdifrs

  # item 1 has the smallest p-value, but another item has the largest statistic
  j <- which.max(st$rdifrs)
  expect_false(j == 1)
  expect_lt(p[1], p[j])

  # only item 1 is flagged, so it is the only item that purification may remove first
  alpha <- (p[1] + p[j]) / 2
  rst <- suppressWarnings(rdif(sim$x, sim$resp, score = sim$score, group = sim$group,
                           focal.name = 1, alpha = alpha, purify = TRUE, verbose = FALSE))
  expect_identical(rst$no_purify$dif_item$rdifrs, 1L)
  expect_true(1 %in% rst$with_purify$dif_item)
  expect_false(j %in% rst$with_purify$dif_item)
})

test_that("rdif() uses a generalized inverse when the covariance matrix is singular", {
  sim <- rdt_sim(miss = 0)
  x2 <- sim$x[1:8, ]
  expect_warning(
    rst <- rdif(x2, sim$resp[, 1:8], score = rep(1, nrow(sim$resp)), group = sim$group, focal.name = "f"),
    "generalized inverse"
  )
  stat <- rst$no_purify$dif_stat

  # with a common ability the squared residual is a linear function of the raw residual
  expect_equal(stat$rdifrs, stat$z.rdifr^2, tolerance = 1e-3)
  expect_equal(stat$p.rdifrs, stat$p.rdifr, tolerance = 1e-3)
  expect_true(all(stat$rdifrs >= 0))
})

test_that("pick_flagged_item() chooses the flagged item with the smallest log p-value", {
  log_p <- c(-1, -50, -3, -4, -10)

  # an unflagged item with a smaller p-value is not chosen
  expect_equal(irtQ:::pick_flagged_item(c(3, 5), log_p), 5)
  expect_equal(irtQ:::pick_flagged_item(c(2, 5), log_p), 2)

  # the first item is chosen when the log p-values are tied
  expect_equal(irtQ:::pick_flagged_item(c(1, 3, 4), c(-2, 0, -7, -7)), 3)
})

test_that("rdif() keeps the log p-values used to choose the item to be removed", {
  skip_on_cran()
  sim <- rdt_sim()
  score <- suppressWarnings(est_score(sim$x, sim$resp, D = 1)$est.theta)
  out <- irtQ:::rdif_one(sim$x, sim$resp, score, sim$group, "f", D = 1)
  expect_equal(exp(out$log_p$rdifr), out$dif_stat$p.rdifr, tolerance = 1e-3, ignore_attr = TRUE)
  expect_equal(exp(out$log_p$rdifs), out$dif_stat$p.rdifs, tolerance = 1e-3, ignore_attr = TRUE)
  expect_equal(exp(out$log_p$rdifrs), out$dif_stat$p.rdifrs, tolerance = 1e-3, ignore_attr = TRUE)
})

test_that("quad_form() uses the inverse for a regular matrix and the rank for a singular one", {
  s <- matrix(c(2, 0.5, 0.5, 1), 2)
  d <- c(0.3, -0.2)
  qf <- irtQ:::quad_form(s, d, df = 2)
  expect_equal(qf$stat, drop(t(d) %*% solve(s) %*% d))
  expect_equal(qf$df, 2)
  expect_false(qf$reduced)

  # a singular matrix of rank one gives the generalized inverse and one degree of freedom
  s1 <- matrix(c(1, 2, 2, 4), 2)
  d1 <- c(1, 2)
  qf1 <- irtQ:::quad_form(s1, d1, df = 2)
  expect_equal(qf1$stat, drop(t(d1) %*% (s1 / 25) %*% d1))
  expect_equal(qf1$df, 1)
  expect_true(qf1$reduced)

  # a matrix with non-finite values gives a missing statistic
  expect_true(is.na(irtQ:::quad_form(matrix(NA_real_, 2, 2), d, df = 2)$stat))
})

test_that("print() of rdif results leaves the significance symbols of skipped items blank", {
  skip_on_cran()
  sim <- rdt_sim()
  score <- suppressWarnings(est_score(sim$x, sim$resp, D = 1)$est.theta)
  rst <- rdif(sim$x, sim$resp, score = score, group = sim$group, focal.name = "f", item.skip = 2)
  out <- utils::capture.output(print(rst))
  expect_false(any(grepl("?", out, fixed = TRUE)))
})

test_that("rdif() returns the flagged items and the number of iterations as integers", {
  skip_on_cran()
  sim <- rdt_sim()
  score <- suppressWarnings(est_score(sim$x, sim$resp, D = 1)$est.theta)
  rst <- rdif(sim$x, sim$resp, score = score, group = sim$group, focal.name = "f", purify = TRUE,
              verbose = FALSE)
  expect_false(is.null(rst$no_purify$dif_item$rdifr))
  for (v in rst$no_purify$dif_item) {
    if (!is.null(v)) expect_type(v, "integer")
  }
  expect_type(rst$with_purify$dif_item, "integer")
  expect_type(rst$with_purify$n.iter, "integer")

  # no item is flagged, so the purification is not carried out
  none <- rdif(sim$x, sim$resp, score = score, group = sim$group, focal.name = "f", alpha = 1e-12,
               purify = TRUE, verbose = FALSE)
  expect_null(none$with_purify$dif_item)
  expect_identical(none$with_purify$n.iter, 0L)
})

test_that("rdif() stops with clear errors for invalid inputs", {
  skip_on_cran()
  sim <- rdt_sim()
  score <- suppressWarnings(est_score(sim$x, sim$resp, D = 1)$est.theta)
  n <- nrow(sim$resp)
  run <- function(...) {
    args <- utils::modifyList(
      list(x = sim$x, data = sim$resp, score = score, group = sim$group, focal.name = "f"), list(...)
    )
    do.call(rdif, args)
  }
  resp_bad <- sim$resp
  resp_bad[1, 1] <- -1
  expect_error(run(data = resp_bad), "score categories")
  resp_bad[1, 1] <- 0.5
  expect_error(run(data = resp_bad), "score categories")
  resp_bad <- sim$resp
  resp_bad[is.na(resp_bad)] <- -9
  expect_error(run(data = resp_bad), "score categories")
  expect_error(run(data = sim$resp[, -1]), "number of columns")
  expect_error(run(score = score[-1]), "'score'")
  expect_error(run(group = sim$group[-1]), "'group'")
  expect_error(run(focal.name = "F"), "focal.name")
  expect_error(run(focal.name = c("f", "r")), "focal.name")
  expect_error(run(group = rep("f", n)), "reference group")
  expect_error(run(group = replace(sim$group, 1:10, "z")), "exactly two groups")
  expect_error(run(item.skip = 20), "item.skip")
  expect_error(run(item.skip = "I2"), "item.skip")
  expect_error(run(item.skip = 1.5), "item.skip")
  expect_error(run(alpha = 5), "alpha")
  expect_error(run(alpha = NA), "alpha")
  expect_error(run(purify = TRUE, max.iter = 0), "max.iter")
  expect_error(run(purify = TRUE, max.iter = 2.5), "max.iter")
  expect_error(run(min.resp = "a"), "min.resp")
})
