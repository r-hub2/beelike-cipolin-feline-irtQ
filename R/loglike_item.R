# Negative log-likelihood (plus negative log priors) of a GRM or GPCM item
# @description This function computes the negative loglikelihood of an item with the
# polytomous IRT model
# @param item_par A vector of item parameters. Unless fix.a = TRUE, the first element is the slope parameter and
# the remaining elements are the threshold parameters; if fix.a = TRUE, all elements are threshold parameters and
# the slope is a.val.
# @param r_i A matrix of (expected) frequencies of the score categories, with one row per theta value and one
# column per category
# @param theta A vector of theta for an item.
# @param pr.mod A character string, "GRM" for the graded response model or "GPCM" for the (generalized) partial
# credit model.
# @param D A scaling constant. A value of 1.702 makes the logistic function approximate the normal ogive
# function. Default is 1.
#
# @return The negative log-likelihood summed over theta values, plus the negative log prior densities of the
# parameters when priors are used
loglike_prm <- function(item_par, r_i, theta, pr.mod = c("GRM", "GPCM"), D = 1,
                        fix.a = FALSE, a.val = 1,
                        aprior = list(dist = "lnorm", params = c(0, 0.5)),
                        bprior = list(dist = "norm", params = c(0.0, 1.0)),
                        use.aprior = FALSE, use.bprior = FALSE,
                        prob_cache = NULL) {
  # `prob_cache`, when non-NULL, is a list produced by make_prm_optim_fns()
  # carrying pre-computed probability objects for the current (item_par, theta):
  #   GRM  -> $P (clipped category probabilities)
  #   GPCM -> $P (= numer / denom, equal to the result of prm() up to rounding)
  # loglike_prm only needs $P, so the same cache structure works for both.
  ## -------------------------------------------------------------------------
  if (!fix.a) {
    # compute category probabilities for all thetas;
    # use the cached P matrix when supplied
    ps <- if (is.null(prob_cache)) {
      prm(theta, a = item_par[1], d = item_par[-1], D = D, pr.model = pr.mod)
    } else {
      prob_cache$P
    }

    # compute log-likelihood
    log_ps <- log(ps)

    # log-likelihood
    llike <- sum(r_i * log_ps)

    # when the slope parameter prior is used
    if (use.aprior) {
      ln.aprior <- logprior(
        val = item_par[1], is.aprior = TRUE, D = D, dist = aprior$dist,
        par.1 = aprior$params[1], par.2 = aprior$params[2]
      )
      llike <- llike + ln.aprior
    }

    # when the difficulty parameter prior is used
    if (use.bprior) {
      ln.bprior <- logprior(
        val = item_par[-1], is.aprior = FALSE, D = NULL, dist = bprior$dist,
        par.1 = bprior$params[1], par.2 = bprior$params[2]
      )
      llike <- llike + sum(ln.bprior)
    }
  } else {
    # compute category probabilities for all thetas (fix.a / PCM path)
    # use the cached P matrix when supplied
    ps <- if (is.null(prob_cache)) {
      prm(theta, a = a.val, d = item_par, D = D, pr.model = pr.mod)
    } else {
      prob_cache$P
    }

    # compute loglikelihood
    log_ps <- log(ps)

    # log-likelihood
    llike <- sum(r_i * log_ps)

    # when the difficulty parameter prior is used
    if (use.bprior) {
      ln.bprior <- logprior(
        val = item_par, is.aprior = FALSE, D = NULL, dist = bprior$dist,
        par.1 = bprior$params[1], par.2 = bprior$params[2]
      )
      llike <- llike + sum(ln.bprior)
    }
  }

  # return a negative log-likelihood
  -llike
}

# Negative log-likelihood (plus negative log priors) of dichotomous items
# @description This function computes the negative loglikelihood of an item with the
# dichotomous IRT model
# @param item_par A vector of item parameters. The first element is the item discrimination (or slope)
# parameter, the second element is the item difficulty parameter, and the third element is the item guessing
# parameter. The third element is necessary only when the 3PLM is used.
# @param r_i A vector of the (expected) numbers of correct responses at each theta value
# @param s_i A vector of the (expected) numbers of incorrect responses at each theta value
# @param theta A vector of theta for an item.
# @param D A scaling constant. A value of 1.702 makes the logistic function approximate the normal ogive
# function. Default is 1.
# @param use.aprior,use.bprior,use.gprior Logical values indicating whether the priors in aprior, bprior, and
# gprior are added to the log-likelihood of the slope, difficulty, and guessing parameters. The guessing prior is
# used only when the guessing parameter is estimated (fix.g = FALSE).
#
# @return The negative log-likelihood summed over theta values, plus the negative log prior densities of the
# parameters when priors are used
loglike_drm <- function(item_par, r_i, s_i, theta, mod = c("1PLM", "2PLM", "3PLM"), D = 1,
                        fix.a = FALSE, fix.g = FALSE, a.val = 1, g.val = .2, n.1PLM = NULL,
                        aprior = list(dist = "lnorm", params = c(0, 0.5)),
                        bprior = list(dist = "norm", params = c(0.0, 1.0)),
                        gprior = list(dist = "beta", params = c(5, 16)),
                        use.aprior = FALSE, use.bprior = FALSE, use.gprior = TRUE,
                        p_cache = NULL) {
  # `p_cache` (when non-NULL) is the drm() probability matrix for the
  # branch picked below, supplied by make_drm_optim_fns() so that the
  # objective / gradient / hessian share one P(theta) per nlminb point.
  # If NULL, each branch calls drm() directly.
  # compute log-likelihood
  # (1) 1PLM: the slope parameters are constrained to be equal across the 1PLM items
  if (!fix.a & mod == "1PLM") {
    # make vectors of a and b parameters for all 1PLM items
    a <- rep(item_par[1], n.1PLM)
    b <- item_par[-1]

    # compute the negative log-likelihood values for all 1PLM items
    llike <- llike_drm(a = a, b = b, g = 0, r_i = r_i, s_i = s_i, theta = theta, D = D, p_cache = p_cache)

    # when the slope parameter prior is used
    if (use.aprior) {
      ln.aprior <- logprior(
        val = item_par[1], is.aprior = TRUE, D = D, dist = aprior$dist,
        par.1 = aprior$params[1], par.2 = aprior$params[2]
      )
      llike <- llike + ln.aprior
    }

    # when the difficulty parameter prior is used
    if (use.bprior) {
      ln.bprior <- logprior(
        val = item_par[-1], is.aprior = FALSE, D = NULL, dist = bprior$dist,
        par.1 = bprior$params[1], par.2 = bprior$params[2]
      )
      llike <- llike + sum(ln.bprior)
    }
  }

  # (2) 1PLM: the slope parameters are fixed to be a specified value
  if (fix.a & mod == "1PLM") {
    # sum of log-likelihood
    llike <- llike_drm(
      a = a.val, b = item_par, g = 0,
      r_i = r_i, s_i = s_i, theta = theta, D = D, p_cache = p_cache
    )

    # when the difficulty parameter prior is used
    if (use.bprior) {
      ln.bprior <- logprior(
        val = item_par, is.aprior = FALSE, D = NULL, dist = bprior$dist,
        par.1 = bprior$params[1], par.2 = bprior$params[2]
      )
      llike <- llike + ln.bprior
    }
  }

  # (3) 2PLM
  if (mod == "2PLM") {
    # sum of log-likelihood
    llike <- llike_drm(
      a = item_par[1], b = item_par[2], g = 0,
      r_i = r_i, s_i = s_i, theta = theta, D = D, p_cache = p_cache
    )

    # when the slope parameter prior is used
    if (use.aprior) {
      ln.aprior <- logprior(
        val = item_par[1], is.aprior = TRUE, D = D, dist = aprior$dist,
        par.1 = aprior$params[1], par.2 = aprior$params[2]
      )
      llike <- llike + ln.aprior
    }

    # when the difficulty parameter prior is used
    if (use.bprior) {
      ln.bprior <- logprior(
        val = item_par[2], is.aprior = FALSE, D = NULL, dist = bprior$dist,
        par.1 = bprior$params[1], par.2 = bprior$params[2]
      )
      llike <- llike + ln.bprior
    }
  }

  # (4) 3PLM
  if (!fix.g & mod == "3PLM") {
    # sum of log-likelihood
    llike <- llike_drm(
      a = item_par[1], b = item_par[2], g = item_par[3],
      r_i = r_i, s_i = s_i, theta = theta, D = D, p_cache = p_cache
    )

    # when the slope parameter prior is used
    if (use.aprior) {
      ln.aprior <- logprior(
        val = item_par[1], is.aprior = TRUE, D = D, dist = aprior$dist,
        par.1 = aprior$params[1], par.2 = aprior$params[2]
      )
      llike <- llike + ln.aprior
    }

    # when the difficulty parameter prior is used
    if (use.bprior) {
      ln.bprior <- logprior(
        val = item_par[2], is.aprior = FALSE, D = NULL, dist = bprior$dist,
        par.1 = bprior$params[1], par.2 = bprior$params[2]
      )
      llike <- llike + ln.bprior
    }

    # when the guessing parameter prior is used
    if (use.gprior) {
      ln.gprior <- logprior(
        val = item_par[3], is.aprior = FALSE, D = NULL, dist = gprior$dist,
        par.1 = gprior$params[1], par.2 = gprior$params[2]
      )
      llike <- llike + ln.gprior
    }
  }

  # (5) 3PLM: the guessing parameters are fixed to be specified value
  if (fix.g & mod == "3PLM") {
    # sum of log-likelihood
    llike <- llike_drm(
      a = item_par[1], b = item_par[2], g = g.val,
      r_i = r_i, s_i = s_i, theta = theta, D = D, p_cache = p_cache
    )

    # when the slope parameter prior is used
    if (use.aprior) {
      ln.aprior <- logprior(
        val = item_par[1], is.aprior = TRUE, D = D, dist = aprior$dist,
        par.1 = aprior$params[1], par.2 = aprior$params[2]
      )
      llike <- llike + ln.aprior
    }

    # when the difficulty parameter prior is used
    if (use.bprior) {
      ln.bprior <- logprior(
        val = item_par[2], is.aprior = FALSE, D = NULL, dist = bprior$dist,
        par.1 = bprior$params[1], par.2 = bprior$params[2]
      )
      llike <- llike + ln.bprior
    }
  }

  # return a negative loglikelihood value
  -llike
}

# compute a sum of the log-likelihood value for each dichotomous item
llike_drm <- function(a, b, g, r_i, s_i, theta, D = 1, p_cache = NULL) {
  # use the cached probability matrix when supplied; otherwise compute
  # drm() directly (the cache holds the drm() result for the same (a, b, g))
  p <- if (is.null(p_cache)) drm(theta, a = a, b = b, g = g, D = D) else p_cache

  # compute 1 - p
  q <- 1 - p

  # compute the log-likelihood
  log_p <- log(p)
  log_q <- log(q)

  # log-likelihood
  L <- r_i * log_p + s_i * log_q

  # sum of log-likelihood
  llike <- sum(L)

  # return
  llike
}
