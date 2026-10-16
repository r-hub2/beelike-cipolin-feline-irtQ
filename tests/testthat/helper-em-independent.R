# Independent IRT formulas for the EM tests (no irtQ probability functions)

# category probabilities (ntheta x cats) of one item
indep_prob <- function(theta, model, pars, cats, D) {
  a <- pars[1]
  if (model %in% c("1PLM", "2PLM", "3PLM")) {
    g <- if (model == "3PLM") pars[3] else 0
    p <- g + (1 - g) / (1 + exp(-D * a * (theta - pars[2])))
    return(cbind(1 - p, p))
  }
  b <- pars[2:cats]
  if (model == "GRM") {
    ps <- matrix(sapply(b, function(bk) 1 / (1 + exp(-D * a * (theta - bk)))), nrow = length(theta))
    return(cbind(1, ps) - cbind(ps, 0))
  }
  z <- matrix(sapply(c(0, b), function(bk) D * a * (theta - bk)), nrow = length(theta))
  cz <- t(apply(z, 1, cumsum))
  e <- exp(cz - apply(cz, 1, max))
  e / rowSums(e)
}

# examinee x quadrature log-likelihood matrix (missing responses skipped)
indep_llmat <- function(x, data, theta, D) {
  data <- as.matrix(data)
  ll <- matrix(0, nrow(data), length(theta))
  for (j in seq_len(nrow(x))) {
    pars <- as.numeric(unlist(x[j, -(1:3)]))
    lp <- log(indep_prob(theta, x$model[j], pars, x$cats[j], D))
    obs <- which(!is.na(data[, j]))
    ll[obs, ] <- ll[obs, ] + t(lp[, data[obs, j] + 1, drop = FALSE])
  }
  ll
}

# log-sum-exp by row
indep_lse <- function(m) {
  mx <- apply(m, 1, max)
  mx + log(rowSums(exp(m - mx)))
}

# marginal log-likelihood per examinee and posterior over the quadrature grid
indep_marg <- function(x, data, weights, D) {
  lj <- sweep(indep_llmat(x, data, weights[, 1], D), 2, log(weights[, 2]), "+")
  li <- indep_lse(lj)
  list(loglik = sum(li), loglik_i = li, post = exp(lj - li))
}

# evaluates an expression and returns its value together with the warning messages it issued
with_warnings <- function(expr) {
  msg <- character()
  val <- withCallingHandlers(expr, warning = function(cond) {
    msg <<- c(msg, conditionMessage(cond))
    invokeRestart("muffleWarning")
  })
  list(value = val, warnings = msg)
}
