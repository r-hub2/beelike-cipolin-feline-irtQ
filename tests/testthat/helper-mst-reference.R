# independent reference computations for the classification and MST tests:
# item category probabilities, the Lord-Wingersky recursion, Lee's (2010)
# indices, Rudner's (2001, 2005) indices, and an exact enumeration of an MST
# panel over stage sum scores

# category probabilities of each item (a list of length(theta) x cats matrices)
ref_prob <- function(x, theta, D = 1) {
  lapply(seq_len(nrow(x)), function(j) {
    a <- x$par.1[j]
    K <- x$cats[j]
    if (K == 2) {
      g <- if (x$model[j] == "3PLM") x$par.3[j] else 0
      p <- g + (1 - g) / (1 + exp(-D * a * (theta - x$par.2[j])))
      cbind(1 - p, p)
    } else {
      b <- unlist(x[j, paste0("par.", 2:K)])
      if (x$model[j] == "GRM") {
        ps <- matrix(sapply(b, function(bk) 1 / (1 + exp(-D * a * (theta - bk)))),
                     nrow = length(theta))
        cs <- cbind(1, ps, 0)
        cs[, 1:K, drop = FALSE] - cs[, 2:(K + 1), drop = FALSE]
      } else {
        z <- sapply(seq_len(K), function(k) {
          if (k == 1) rep(0, length(theta)) else
            rowSums(matrix(sapply(b[1:(k - 1)], function(bk) D * a * (theta - bk)),
                           nrow = length(theta)))
        })
        e <- exp(matrix(z, nrow = length(theta)))
        e / rowSums(e)
      }
    }
  })
}

# summed-score distributions (rows = scores 0..max, columns = theta values)
ref_lw <- function(x, theta, D = 1) {
  pl <- ref_prob(x, theta, D)
  sapply(seq_along(theta), function(t) {
    d <- 1
    for (pj in pl) {
      p <- pj[t, ]
      nd <- numeric(length(d) + length(p) - 1)
      for (k in seq_along(p)) {
        idx <- k:(k + length(d) - 1)
        nd[idx] <- nd[idx] + d * p[k]
      }
      d <- nd
    }
    d
  })
}

# Lee (2010): conditional and marginal accuracy and consistency
ref_lee <- function(x, cut, theta, w, D = 1) {
  lk <- ref_lw(x, theta, D)
  sc <- 0:(nrow(lk) - 1)
  tau <- colSums(lk * sc)
  K <- length(cut) + 1
  P <- matrix(sapply(1:K, function(h) colSums(lk[findInterval(sc, cut) + 1 == h, , drop = FALSE])),
              ncol = K)
  lev <- findInterval(tau, cut) + 1
  acc <- P[cbind(seq_along(theta), lev)]
  con <- rowSums(P^2)
  list(acc = sum(w * acc), con = sum(w * con), cond_acc = acc, cond_con = con, P = P)
}

# Rudner (2001, 2005): conditional and marginal accuracy and consistency
ref_rud <- function(cut, theta, se, w) {
  br <- c(-Inf, cut, Inf)
  P <- t(sapply(seq_along(theta), function(i) diff(pnorm(br, theta[i], se[i]))))
  lev <- findInterval(theta, cut) + 1
  acc <- P[cbind(seq_along(theta), lev)]
  con <- rowSums(P^2)
  list(acc = sum(w * acc), con = sum(w * con), cond_acc = acc, cond_con = con)
}

# exact evaluation of an MST panel by enumerating every sequence of stage
# sum scores; routing uses the inverse TCC estimate of the cumulative sum
# score and the cut scores between the modules of the next stage in module
# index order (a tie goes to the higher module); cfg lists the modules of each
# stage so that the reference does not depend on panel_info()
ref_brute_mst <- function(x, module, route_map, cut_score, theta, cfg, D = 1,
                          range.tcc = c(-7, 7)) {
  n.stg <- length(cfg)
  items <- lapply(seq_len(ncol(module)), function(m) which(module[, m] == 1))
  cache <- new.env()
  eqth <- function(mods) {
    key <- paste(mods, collapse = "_")
    if (is.null(cache[[key]])) {
      cache[[key]] <- irtQ:::inv_tcc_nr(confirm_df(x[unlist(items[mods]), ]), D = D,
                                        range.tcc = range.tcc)$est.theta
    }
    cache[[key]]
  }
  out <- t(sapply(theta, function(th) {
    pm <- lapply(items, function(it) ref_lw(x[it, ], th, D)[, 1])
    m1 <- cfg[[1]][1]
    states <- lapply(seq_along(pm[[m1]]), function(k) list(mods = m1, sc = k - 1, p = pm[[m1]][k]))
    for (s in 2:n.stg) {
      new <- list()
      for (st in states) {
        nxt <- which(route_map[utils::tail(st$mods, 1), ] == 1)
        est <- eqth(st$mods)[st$sc + 1]
        pos <- match(nxt, sort(cfg[[s]]))
        cu <- cut_score[[s - 1]][pos[-length(pos)]]
        m <- nxt[sum(est >= cu) + 1]
        for (k in seq_along(pm[[m]])) {
          new[[length(new) + 1]] <- list(mods = c(st$mods, m), sc = st$sc + k - 1,
                                         p = st$p * pm[[m]][k])
        }
      }
      states <- new
    }
    fin <- vapply(states, function(st) eqth(st$mods)[st$sc + 1], numeric(1))
    p <- vapply(states, function(st) st$p, numeric(1))
    mu <- sum(p * fin)
    c(mu = mu, sigma2 = sum(p * fin^2) - mu^2, ptot = sum(p))
  }))
  data.frame(theta = theta, out)
}

# a small mixed-format 1-2-3 panel: three 3PLM items and one GRM or GPCM item
# per module, so that every module has the same maximum sum score of 5
ref_panel <- function() {
  set.seed(2028)
  n_mod <- 6L
  x <- shape_df(
    par.drm = list(a = runif(3 * n_mod, 0.8, 1.6), b = rnorm(3 * n_mod, 0, 0.8),
                   g = rep(0.15, 3 * n_mod)),
    par.prm = list(a = runif(n_mod, 0.8, 1.4),
                   d = replicate(n_mod, sort(rnorm(2, 0, 0.7)), simplify = FALSE)),
    cats = c(rep(2L, 3 * n_mod), rep(3L, n_mod)),
    model = c(rep("3PLM", 3 * n_mod), rep(c("GRM", "GPCM"), length.out = n_mod))
  )
  module <- matrix(0L, nrow(x), n_mod)
  for (m in seq_len(n_mod)) module[c(((m - 1) * 3 + 1):(m * 3), 3 * n_mod + m), m] <- 1L
  route_map <- matrix(0L, n_mod, n_mod)
  route_map[1, 2:3] <- 1L
  route_map[2, 4:5] <- 1L
  route_map[3, 5:6] <- 1L
  list(x = x, module = module, route_map = route_map, cut_score = list(0.1, c(-0.6, 0.7)))
}
