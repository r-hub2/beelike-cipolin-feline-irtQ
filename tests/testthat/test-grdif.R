# Category probabilities of the items computed directly from the item parameters
grt_probs <- function(x, theta, D) {
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

# Independent implementation of the GRDIF statistics (Lim, Zhu, Choe, & Han, 2024, Eqs. 1-10),
# with the expected item score residual for polytomous items
grt_ref <- function(x, resp, score, group, focal, D) {
  ref <- setdiff(unique(group[!is.na(group)]), focal)
  lev <- c(ref, focal)
  G <- length(lev)
  P <- grt_probs(x, ifelse(is.na(score), 0, score), D)
  C <- cbind(-1, diag(G - 1))
  Crs <- kronecker(C, diag(2))
  out <- NULL
  for (j in seq_len(nrow(x))) {
    k <- 0:(x$cats[j] - 1)
    ey <- drop(P[[j]] %*% k)
    dev <- outer(-ey, k, "+")
    m2 <- rowSums(dev^2 * P[[j]])
    m3 <- rowSums(dev^3 * P[[j]])
    m4 <- rowSums(dev^4 * P[[j]])
    mrr <- msr <- mus <- vr <- vs <- cv <- nn <- numeric(G)
    for (g in seq_len(G)) {
      ok <- which(!is.na(group) & group == lev[g] & !is.na(resp[, j]) & !is.na(score))
      n <- length(ok)
      e <- resp[ok, j] - ey[ok]
      mrr[g] <- mean(e)
      msr[g] <- mean(e^2)
      mus[g] <- mean(m2[ok])
      vr[g] <- sum(m2[ok]) / n^2
      vs[g] <- sum(m4[ok] - m2[ok]^2) / n^2
      cv[g] <- sum(m3[ok]) / n^2
      nn[g] <- n
    }
    Srs <- matrix(0, 2 * G, 2 * G)
    for (g in seq_len(G)) {
      Srs[2 * g - 1:0, 2 * g - 1:0] <- matrix(c(vr[g], cv[g], cv[g], vs[g]), 2)
    }
    dr <- C %*% mrr
    ds <- C %*% (msr - mus)
    drs <- Crs %*% (as.vector(rbind(mrr, msr)) - as.vector(rbind(0, mus)))
    chi_r <- drop(t(dr) %*% solve(C %*% diag(vr, G) %*% t(C), dr))
    chi_s <- drop(t(ds) %*% solve(C %*% diag(vs, G) %*% t(C), ds))
    chi_rs <- drop(t(drs) %*% solve(Crs %*% Srs %*% t(Crs), drs))
    out <- rbind(out, data.frame(
      grdifr = chi_r, grdifs = chi_s, grdifrs = chi_rs,
      p.grdifr = pchisq(chi_r, G - 1, lower.tail = FALSE),
      p.grdifs = pchisq(chi_s, G - 1, lower.tail = FALSE),
      p.grdifrs = pchisq(chi_rs, 2 * (G - 1), lower.tail = FALSE),
      n.ref = nn[1], n.total = sum(nn), msr.ref = mus[1], var.mrr = vr[1], var.msr = vs[1], cov = cv[1]
    ))
  }
  out
}

# Simulate dichotomous data for a reference group and three focal groups with DIF on items 1, 2, 3,
# and 11
grt_sim <- function(seed = 1, n = c(400, 300, 300, 250), miss = 0, dif10 = FALSE) {
  set.seed(seed)
  x <- shape_df(
    par.drm = list(a = runif(11, 0.8, 1.8), b = rnorm(11), g = c(rep(NA, 4), runif(7, 0.05, 0.2))),
    item.id = paste0("I", 1:11), cats = 2, model = c(rep("2PLM", 4), rep("3PLM", 7))
  )
  lab <- c("r", "f1", "f2", "f3")[seq_along(n)]
  resp <- NULL
  for (k in seq_along(n)) {
    xk <- x
    if (k > 1) {
      xk$par.2[1:2] <- xk$par.2[1:2] + c(0.6, -0.5, 0.4)[k - 1]
      if (k == 2) xk$par.1[3] <- xk$par.1[3] * 0.4
      if (k == 3) xk$par.2[11] <- xk$par.2[11] + 0.8
      if (dif10 && k == 4) xk$par.2[10] <- xk$par.2[10] + 1
    }
    resp <- rbind(resp, simdat(xk, theta = rnorm(n[k], c(0, -0.4, 0.3, 0.1)[k]), D = 1))
  }
  if (miss > 0) resp[sample(length(resp), round(miss * length(resp)))] <- NA
  list(x = x, resp = resp, group = rep(lab, n), focal = lab[-1])
}

# Simulate dichotomous data of three items for three groups, in which the first focal group shows
# large DIF on the first two items
grt_sim3 <- function(seed = 5) {
  x <- shape_df(par.drm = list(a = c(1, 1.2, 0.9), b = c(0, 0.5, -0.5), g = 0), cats = 2, model = "2PLM")
  x2 <- x
  x2$par.2 <- x2$par.2 + c(1.5, 1.2, 0)
  set.seed(seed)
  resp <- rbind(simdat(x, theta = rnorm(500), D = 1), simdat(x2, theta = rnorm(500), D = 1),
                simdat(x, theta = rnorm(500), D = 1))
  list(x = x, resp = resp, group = rep(c(0, 1, 2), each = 500))
}

# Evaluate an expression and return its value together with all warning messages
grt_catch <- function(expr) {
  msgs <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    msgs <<- c(msgs, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = msgs)
}

test_that("grdif() statistics match an independent implementation", {
  skip_on_cran()
  sim <- grt_sim(miss = 0.25)
  for (D in c(1, 1.702)) {
    score <- suppressWarnings(est_score(sim$x, sim$resp, D = D)$est.theta)
    score[c(5, 450, 900)] <- NA
    out <- grt_catch(grdif(sim$x, sim$resp, score = score, group = sim$group, focal.name = sim$focal,
                           D = D, post.hoc = FALSE))
    rst <- out$value
    ref <- grt_ref(sim$x, sim$resp, score, sim$group, sim$focal, D)
    ds <- rst$no_purify$dif_stat
    expect_equal(ds$grdifr, ref$grdifr, tolerance = 1e-4)
    expect_equal(ds$grdifs, ref$grdifs, tolerance = 1e-4)
    expect_equal(ds$grdifrs, ref$grdifrs, tolerance = 1e-4)
    expect_lt(max(abs(ds[, c("p.grdifr", "p.grdifs", "p.grdifrs")] -
                        ref[, c("p.grdifr", "p.grdifs", "p.grdifrs")])), 5.1e-5)
    expect_equal(ds$n.ref, ref$n.ref)
    expect_equal(ds$n.total, ref$n.total)
    expect_equal(rst$no_purify$moments$mu$msr.ref, ref$msr.ref, tolerance = 1e-8)
    expect_equal(rst$no_purify$moments$var$mrr.ref, ref$var.mrr, tolerance = 1e-8)
    expect_equal(rst$no_purify$moments$var$msr.ref, ref$var.msr, tolerance = 1e-8)
    expect_equal(rst$no_purify$moments$covariance$ref, ref$cov, tolerance = 1e-8)
    expect_identical(as.integer(rst$no_purify$dif_item$grdifrs), which(ref$p.grdifrs <= 0.05))
    expect_identical(out$warnings, "3 examinee(s) with a missing ability estimate were excluded.")
  }
})

test_that("grdif() methods for est_irt and est_item objects match the default method", {
  skip_on_cran()
  sim <- grt_sim3()
  fit <- suppressWarnings(est_irt(data = sim$resp, D = 1, model = "2PLM", cats = 2, verbose = FALSE))
  r1 <- grdif(fit, group = sim$group, focal.name = c(1, 2), purify = TRUE, verbose = FALSE)
  r0 <- grdif(fit$par.est, sim$resp, group = sim$group, focal.name = c(1, 2), D = 1, purify = TRUE,
              verbose = FALSE)
  r1$call <- r0$call <- NULL
  expect_identical(r1, r0)
  score <- est_score(fit$par.est, sim$resp, D = 1)$est.theta
  fit2 <- suppressWarnings(est_item(x = fit$par.est, data = sim$resp, score = score, D = 1, verbose = FALSE))
  r2 <- grdif(fit2, group = sim$group, focal.name = c(1, 2))
  r3 <- grdif(fit2$par.est, sim$resp, score = score, group = sim$group, focal.name = c(1, 2), D = 1)
  r2$call <- r3$call <- NULL
  expect_identical(r2, r3)
  expect_error(grdif(fit, data = sim$resp, group = sim$group, focal.name = c(1, 2)), "taken from the object")
  expect_error(grdif(fit2, score = score, group = sim$group, focal.name = c(1, 2)), "taken from the object")
})

test_that("grdif() excludes examinees with a missing ability estimate or group with a warning", {
  skip_on_cran()
  sim <- grt_sim(seed = 2)
  score <- suppressWarnings(est_score(sim$x, sim$resp, D = 1)$est.theta)
  r0 <- grdif(sim$x, sim$resp, score = score, group = sim$group, focal.name = sim$focal, post.hoc = FALSE)

  # the examinees with a missing group are not part of the reference group
  group <- sim$group
  group[401:450] <- NA
  r1 <- grdif(sim$x, sim$resp, score = score, group = group, focal.name = sim$focal, post.hoc = FALSE)
  expect_equal(r1$no_purify$dif_stat$n.ref, r0$no_purify$dif_stat$n.ref)
  expect_equal(r1$no_purify$dif_stat$n.foc1, r0$no_purify$dif_stat$n.foc1 - 50)

  # a missing ability estimate gives the same results as removing the examinee
  s2 <- score
  s2[c(3, 600)] <- NA
  out <- grt_catch(grdif(sim$x, sim$resp, score = s2, group = sim$group, focal.name = sim$focal))
  r3 <- grdif(sim$x, sim$resp[-c(3, 600), ], score = score[-c(3, 600)], group = sim$group[-c(3, 600)],
              focal.name = sim$focal)
  expect_equal(out$value$no_purify$dif_stat, r3$no_purify$dif_stat)
  expect_equal(out$value$no_purify$post.hoc, r3$no_purify$post.hoc)
  expect_identical(out$warnings, "2 examinee(s) with a missing ability estimate were excluded.")

  # an examinee without any response is excluded without a warning when the scores are estimated
  resp <- sim$resp
  resp[7, ] <- NA
  expect_no_warning(grdif(sim$x, resp, group = sim$group, focal.name = sim$focal, post.hoc = FALSE))
  expect_error(grdif(sim$x, sim$resp, score = rep(NA_real_, length(score)), group = sim$group,
                     focal.name = sim$focal), "No examinee has an ability estimate")
})

test_that("grdif() applies min.resp also when the scores are supplied", {
  skip_on_cran()
  sim <- grt_sim(seed = 3, miss = 0.4)
  score <- suppressWarnings(est_score(sim$x, sim$resp, D = 1)$est.theta)
  few <- rowSums(!is.na(sim$resp)) < 5 & rowSums(!is.na(sim$resp)) > 0
  out <- grt_catch(grdif(sim$x, sim$resp, score = score, group = sim$group, focal.name = sim$focal,
                         min.resp = 5, post.hoc = FALSE))
  r0 <- grdif(sim$x, sim$resp[!few, ], score = score[!few], group = sim$group[!few],
              focal.name = sim$focal, post.hoc = FALSE)
  expect_equal(out$value$no_purify$dif_stat, r0$no_purify$dif_stat)
  expect_true(all(is.na(out$value$no_purify$score[few])))
  expect_identical(out$warnings, paste0(sum(few), " examinee(s) with fewer than 5 responses were excluded."))
})

test_that("grdif() gives missing statistics for an item that a group did not answer", {
  skip_on_cran()
  sim <- grt_sim(seed = 4)
  score <- suppressWarnings(est_score(sim$x, sim$resp, D = 1)$est.theta)
  resp <- sim$resp
  resp[sim$group == "f2", 4] <- NA
  rst <- grdif(sim$x, resp, score = score, group = sim$group, focal.name = sim$focal)
  ds <- rst$no_purify$dif_stat
  expect_true(all(is.nan(unlist(ds[4, c("grdifr", "grdifs", "grdifrs")]))))
  expect_equal(ds$n.foc2[4], 0)
  expect_false(4 %in% unlist(rst$no_purify$dif_item))
  expect_false(anyNA(ds$grdifrs[-4]))
})

test_that("grdif() stops for a scoring method without pattern scores when it estimates abilities", {
  skip_on_cran()
  sim <- grt_sim3()
  score <- est_score(sim$x, sim$resp, D = 1)$est.theta
  for (m in c("EAP.SUM", "INV.TCC")) {
    expect_error(grdif(sim$x, sim$resp, group = sim$group, focal.name = c(1, 2), method = m),
                 "not supported by grdif()", fixed = TRUE)
    expect_error(grdif(sim$x, sim$resp, score = score, group = sim$group, focal.name = c(1, 2),
                       method = m, purify = TRUE, verbose = FALSE),
                 "not supported by grdif()", fixed = TRUE)
    expect_no_error(grdif(sim$x, sim$resp, score = score, group = sim$group, focal.name = c(1, 2),
                          method = m))
  }
})

test_that("grdif() reads character responses as numbers", {
  sim <- grt_sim(seed = 7)
  score <- suppressWarnings(est_score(sim$x, sim$resp, D = 1)$est.theta)
  r0 <- grdif(sim$x, sim$resp, score = score, group = sim$group, focal.name = sim$focal, post.hoc = FALSE)
  chr <- as.data.frame(apply(sim$resp, 2, as.character), stringsAsFactors = FALSE)
  r1 <- grdif(sim$x, chr, score = score, group = sim$group, focal.name = sim$focal, post.hoc = FALSE)
  expect_identical(r1$no_purify$dif_stat, r0$no_purify$dif_stat)
})


test_that("grdif() stops with clear errors for invalid inputs", {
  sim <- grt_sim(seed = 8)
  score <- suppressWarnings(est_score(sim$x, sim$resp, D = 1)$est.theta)
  run <- function(...) {
    args <- utils::modifyList(
      list(x = sim$x, data = sim$resp, score = score, group = sim$group, focal.name = sim$focal,
           post.hoc = FALSE), list(...)
    )
    do.call(grdif, args)
  }
  resp_bad <- sim$resp
  resp_bad[1, 1] <- -1
  expect_error(run(data = resp_bad), "score categories")
  resp_bad[1, 1] <- 2
  expect_error(run(data = resp_bad), "score categories")
  expect_error(run(data = sim$resp[, -1]), "number of columns")
  expect_error(run(score = score[-1]), "'score'")
  expect_error(run(group = sim$group[-1]), "'group'")
  expect_error(run(focal.name = c("f1", "f9")), "focal.name")
  expect_error(run(focal.name = c("f1", "f1")), "focal.name")
  expect_error(run(focal.name = c("r", sim$focal)), "reference group")
  expect_error(run(group = replace(sim$group, sim$group == "f3", "z"), focal.name = c("f1", "f2")),
               "exactly one reference group")
  expect_error(run(alpha = 5), "alpha")
  expect_error(run(alpha = NA), "alpha")
  expect_error(run(purify = TRUE, max.iter = 0), "max.iter")
  expect_error(run(purify = TRUE, max.iter = 1.5), "max.iter")
  expect_error(run(min.resp = "a"), "min.resp")
})

test_that("grdif() stops when an item is polytomous", {
  skip_on_cran()
  sim <- grt_sim(seed = 9)
  x <- sim$x
  x$cats[10] <- 3L
  x$model[10] <- "GRM"
  x$par.3[10] <- 1
  expect_error(
    grdif(x, sim$resp, group = sim$group, focal.name = sim$focal),
    "supports dichotomous items only"
  )
  fit <- suppressWarnings(est_irt(
    data = sim$resp[, 1:3], D = 1, model = "2PLM", cats = 2, verbose = FALSE
  ))
  expect_no_error(grdif(fit, group = sim$group, focal.name = sim$focal, post.hoc = FALSE))

  # an item of the "GRM" or "GPCM" model is refused also when it has two categories
  x2 <- sim$x
  x2$model[1] <- "GRM"
  expect_error(
    grdif(x2, sim$resp, group = sim$group, focal.name = sim$focal),
    "supports dichotomous items only"
  )
  x2$model[1] <- "GPCM"
  expect_error(
    grdif(x2, sim$resp, group = sim$group, focal.name = sim$focal),
    "supports dichotomous items only"
  )

  # the objects with polytomous items are refused
  set.seed(3)
  xp <- shape_df(
    par.prm = list(a = c(1, 1.2, 0.9, 1.1), d = list(c(-1, 0.5), c(-0.5, 1), c(-0.8, 0.8), c(-1.2, 0.2))),
    cats = 3, model = "GRM"
  )
  resp_p <- simdat(xp, theta = rnorm(300), D = 1)
  group_p <- rep(0:1, each = 150)
  fit_p <- suppressWarnings(est_irt(data = resp_p, D = 1, model = "GRM", cats = 3, verbose = FALSE))
  expect_error(grdif(fit_p, group = group_p, focal.name = 1), "supports dichotomous items only")
  score_p <- est_score(fit_p$par.est, resp_p, D = 1)$est.theta
  fit_p2 <- suppressWarnings(est_item(
    x = fit_p$par.est, data = resp_p, score = score_p, D = 1, verbose = FALSE
  ))
  expect_error(grdif(fit_p2, group = group_p, focal.name = 1), "supports dichotomous items only")
})

test_that("grdif() warns once about the number of cores during the purification", {
  skip_on_cran()
  sim <- grt_sim(miss = 0.2)
  out <- grt_catch(grdif(sim$x, sim$resp, group = sim$group, focal.name = sim$focal, purify = TRUE,
                         post.hoc = FALSE, ncore = 2, verbose = FALSE))
  expect_gte(out$value$with_purify$n.iter, 2)
  expect_equal(sum(grepl("^ncore > 1", out$warnings)), 1)
})

test_that("grdif() flags items with the unrounded p-values", {
  sim <- grt_sim(seed = 6)
  score <- suppressWarnings(est_score(sim$x, sim$resp, D = 1)$est.theta)
  ref <- grt_ref(sim$x, sim$resp, score, sim$group, sim$focal, 1)
  j <- which(abs(ref$p.grdifrs - round(ref$p.grdifrs, 4)) > 1e-5 & ref$p.grdifrs > 0.01)[1]
  alpha <- (ref$p.grdifrs[j] + round(ref$p.grdifrs[j], 4)) / 2
  rst <- grdif(sim$x, sim$resp, score = score, group = sim$group, focal.name = sim$focal,
               alpha = alpha, post.hoc = FALSE)
  expect_identical(j %in% rst$no_purify$dif_item$grdifrs, ref$p.grdifrs[j] <= alpha)
})

test_that("grdif() uses a generalized inverse when a covariance matrix is singular", {
  sim <- grt_sim(seed = 5)

  # every examinee has the same ability estimate, so the MSR is a linear function of the MRR
  out <- grt_catch(grdif(sim$x, sim$resp, score = rep(0.5, nrow(sim$resp)), group = sim$group,
                         focal.name = sim$focal, post.hoc = FALSE))
  ds <- out$value$no_purify$dif_stat
  expect_true(any(grepl("singular or nearly singular", out$warnings)))
  expect_equal(ds$grdifrs, ds$grdifr, tolerance = 1e-3)
  expect_equal(ds$p.grdifrs, ds$p.grdifr, tolerance = 1e-3)
  expect_false(anyNA(ds$p.grdifrs))
})

test_that("grdif() purification removes the flagged item with the smallest p-value first", {
  skip_on_cran()
  sim <- grt_sim()
  score <- suppressWarnings(est_score(sim$x, sim$resp, D = 1)$est.theta)
  rst <- grdif(sim$x, sim$resp, score = score, group = sim$group, focal.name = sim$focal,
               purify = TRUE, purify.by = "grdifr", post.hoc = FALSE, verbose = FALSE)
  first <- rst$with_purify$dif_stat$n.iter == 0
  p0 <- rst$no_purify$dif_stat$grdifr
  expect_equal(which(first), which.max(p0))
  expect_true(rst$with_purify$complete)

  # the final scores are those of the items that were not flagged
  keep <- setdiff(seq_len(nrow(sim$x)), rst$with_purify$dif_item)
  sc <- est_score(trim_par_cols(sim$x[keep, ]), sim$resp[, keep], D = 1, range = c(-5, 5))$est.theta
  expect_equal(rst$with_purify$score, sc)
  expect_type(rst$with_purify$n.iter, "integer")
})

test_that("grdif() purification stops with a warning when all items are flagged", {
  skip_on_cran()
  sim <- grt_sim3()
  out <- grt_catch(grdif(sim$x[1:2, ], sim$resp[, 1:2], group = sim$group, focal.name = c(1, 2),
                         purify = TRUE, alpha = 0.999, verbose = FALSE))
  expect_true(any(grepl("All items were flagged", out$warnings)))
  expect_false(out$value$with_purify$complete)
  expect_equal(out$value$with_purify$dif_item, 1:2)
  expect_identical(out$value$with_purify$n.iter, 1L)
  expect_false(anyNA(out$value$with_purify$dif_stat$grdifrs))
  expect_setequal(unique(out$value$with_purify$post.hoc$id), sim$x$id[1:2])
  r0 <- grdif(sim$x[1:2, ], sim$resp[, 1:2], group = sim$group, focal.name = c(1, 2), alpha = 0.999)
  expect_identical(out$value$no_purify, r0$no_purify)
})

test_that("grdif() purification stops with a warning when a group has no ability estimate left", {
  skip_on_cran()
  sim <- grt_sim(seed = 3, miss = 0.3)
  score <- suppressWarnings(est_score(sim$x, sim$resp, D = 1)$est.theta)
  args <- list(x = sim$x, data = sim$resp, score = score, group = sim$group, focal.name = sim$focal,
               min.resp = 8, purify = TRUE, alpha = 0.9, post.hoc = FALSE, verbose = FALSE)
  out <- grt_catch(do.call(grdif, args))
  expect_true(any(grepl("no examinee of a group", out$warnings)))
  expect_false(out$value$with_purify$complete)
  expect_false(anyNA(out$value$with_purify$dif_stat$grdifrs))
  expect_false(anyNA(out$value$with_purify$moments$mu$msr.ref))

  # the returned results are those of the last iteration counted
  out2 <- grt_catch(do.call(grdif, utils::modifyList(args, list(max.iter = out$value$with_purify$n.iter))))
  keep <- out$value$with_purify$dif_stat$n.iter == out$value$with_purify$n.iter
  expect_identical(out$value$with_purify$dif_stat[keep, ], out2$value$with_purify$dif_stat[keep, ])
  expect_identical(out$value$with_purify$score, out2$value$with_purify$score)
  for (g in unique(sim$group)) {
    expect_false(all(is.na(out$value$with_purify$score[sim$group == g])))
  }
})

test_that("grdif() reports an incomplete purification when max.iter is reached", {
  skip_on_cran()
  sim <- grt_sim()
  score <- suppressWarnings(est_score(sim$x, sim$resp, D = 1)$est.theta)
  out <- grt_catch(grdif(sim$x, sim$resp, score = score, group = sim$group, focal.name = sim$focal,
                         purify = TRUE, max.iter = 1, post.hoc = FALSE, verbose = FALSE))
  expect_true(any(grepl("maximum number of iterations", out$warnings)))
  expect_false(out$value$with_purify$complete)
  expect_true(length(out$value$with_purify$dif_item) > 1)
})


test_that("grdif() post-hoc analyses match rdif() for each pair, with the reference group as the reference of its pairs", {
  skip_on_cran()
  sim <- grt_sim()
  score <- suppressWarnings(est_score(sim$x, sim$resp, D = 1)$est.theta)

  # the reference group is sorted after the focal groups, and the focal groups are not in sorted order
  group <- unname(c(r = "z", f1 = "c", f2 = "a", f3 = "b")[sim$group])
  focal <- c("c", "a", "b")
  rst <- grdif(sim$x, sim$resp, score = score, group = group, focal.name = focal)
  flag <- rst$no_purify$dif_item$grdifrs
  expect_gt(length(flag), 0)
  ph <- rst$no_purify$post.hoc$by.grdifrs
  pairs <- utils::combn(sort(unique(group)), 2, simplify = FALSE)
  expect_equal(nrow(ph), length(flag) * length(pairs))
  for (p in pairs) {
    p <- p[order(match(p, c("z", focal)))]
    keep <- group %in% p
    two <- rdif(sim$x[flag, ], sim$resp[keep, flag, drop = FALSE], score = score[keep],
                group = group[keep], focal.name = p[2])
    got <- ph[ph$group.pair == paste(p, collapse = " & "), ]
    rownames(got) <- NULL
    expect_equal(got$id, sim$x$id[flag])
    expect_equal(got[, setdiff(names(got), "group.pair")], two$no_purify$dif_stat)
  }
  expect_true(all(c("z & c", "z & a", "z & b", "c & a", "c & b", "a & b") %in% ph$group.pair))

  # a factor group gives the same pairs and signs as the character group
  rst_f <- grdif(sim$x, sim$resp, score = score, group = factor(group), focal.name = focal)
  expect_equal(rst_f$no_purify$post.hoc, rst$no_purify$post.hoc)
})

test_that("print() of grdif results leaves the significance symbols of missing statistics blank", {
  sim <- grt_sim(seed = 4)
  score <- suppressWarnings(est_score(sim$x, sim$resp, D = 1)$est.theta)
  resp <- sim$resp
  resp[sim$group == "f2", 4] <- NA
  rst <- grdif(sim$x, resp, score = score, group = sim$group, focal.name = sim$focal)
  expect_false(any(grepl("?", utils::capture.output(print(rst)), fixed = TRUE)))
})
