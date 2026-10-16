# This function analytically computes the gradient vector of the negative log-likelihood (plus the negative log
# priors) of dichotomous item parameters
grad_item_drm <- function(item_par, f_i, r_i, theta, mod = c("1PLM", "2PLM", "3PLM"), D = 1,
                          fix.a = FALSE, fix.g = TRUE, a.val = 1, g.val = .2, n.1PLM = NULL,
                          aprior = list(dist = "lnorm", params = c(0, 0.5)),
                          bprior = list(dist = "norm", params = c(0.0, 1.0)),
                          gprior = list(dist = "beta", params = c(5, 16)),
                          use.aprior = FALSE,
                          use.bprior = FALSE,
                          use.gprior = TRUE,
                          p_cache = NULL) {
  # `p_cache`: optional drm() probability matrix for the active branch.
  # When supplied by make_drm_optim_fns(), the same P(theta) is reused
  # across loglike_drm / grad_item_drm / hess_item_drm at the same
  # (item_par, theta), which saves two of the three drm() calls per nlminb point.
  # When NULL (e.g. callers outside the optim factory), each branch
  # calls drm() directly.
  # count the number of item parameters to be estimated
  n.par <- length(item_par)

  # (1) 1PLM: the slope parameters are constrained to be equal across the 1PLM items
  if (!fix.a & mod == "1PLM") {
    # make vectors of a and b parameters for all 1PLM items
    a <- rep(item_par[1], n.1PLM)
    b <- item_par[-1]

    # compute the probabilities of correct and incorrect
    # use the cached P matrix when supplied
    p <- if (is.null(p_cache)) drm(theta = theta, a = a, b = b, g = 0, D = D) else p_cache

    # compute the component values
    r_p <- r_i - (f_i * p)
    theta_b <- -Rfast::Outer(x = b, y = theta, oper = "-")

    # compute the gradients of a and bs parameters
    gr_a <- -D * sum(theta_b * r_p)
    gr_b <- (D * a[1]) * Rfast::colsums(r_p)

    # combine all gradients into a vector
    grad <- c(gr_a, gr_b)

    # create a prior gradient vector
    grad.prior <- rep(0, n.par)

    # extract the gradient vector when the slope parameter prior is used
    if (use.aprior) {
      rst.aprior <-
        logprior_deriv(
          val = a[1], is.aprior = TRUE, D = D, dist = aprior$dist,
          par.1 = aprior$params[1], par.2 = aprior$params[2]
        )

      # extract the gradient vector
      grad.prior[1] <- attributes(rst.aprior)$gradient
    }

    # extract the gradient vector when the difficulty parameter prior is used
    if (use.bprior) {
      rst.bprior <-
        logprior_deriv(
          val = b, is.aprior = FALSE, D = NULL, dist = bprior$dist,
          par.1 = bprior$params[1], par.2 = bprior$params[2]
        )

      # extract the gradient vector
      grad.prior[-1] <- attributes(rst.bprior)$gradient
    }
  }

  # (2) 1PLM: the slope parameters are fixed to be a specified value
  if (fix.a & mod == "1PLM") {
    # assign a and b parameters
    a <- a.val
    b <- item_par

    # compute the probabilities of correct and incorrect
    # (reuse the cached P matrix when supplied by the optim factory)
    p <- if (is.null(p_cache)) drm(theta = theta, a = a, b = b, g = 0, D = D) else p_cache

    # compute the component values
    r_p <- r_i - (f_i * p)

    # compute the gradients of a and bs parameters
    gr_b <- -as.numeric(-D * a * sum(r_p))

    # combine all gradients into a vector
    grad <- gr_b

    # create a prior gradient vector
    grad.prior <- rep(0, n.par)

    # extract the gradient vector when the difficulty parameter prior is used
    if (use.bprior) {
      rst.bprior <-
        logprior_deriv(
          val = b, is.aprior = FALSE, D = NULL, dist = bprior$dist,
          par.1 = bprior$params[1], par.2 = bprior$params[2]
        )

      # extract the gradient vector
      grad.prior <- attributes(rst.bprior)$gradient
    }
  }

  # (3) 2PLM
  if (mod == "2PLM") {
    # assign a and b parameters
    a <- item_par[1]
    b <- item_par[2]

    # compute the probabilities of correct and incorrect
    # (reuse the cached P matrix when supplied by the optim factory)
    p <- if (is.null(p_cache)) drm(theta = theta, a = a, b = b, g = 0, D = D) else p_cache

    # compute the component values
    Da <- D * a
    r_p <- r_i - (f_i * p)
    theta_b <- theta - b

    # compute the gradients of a and bs parameters
    gr_a <- -D * sum(theta_b * r_p)
    gr_b <- Da * sum(r_p)

    # combine all gradients into a vector
    grad <- c(gr_a, gr_b)

    # create a prior gradient vector
    grad.prior <- rep(0, n.par)

    # extract the gradient vector when the slope parameter prior is used
    if (use.aprior) {
      rst.aprior <-
        logprior_deriv(
          val = a, is.aprior = TRUE, D = D, dist = aprior$dist,
          par.1 = aprior$params[1], par.2 = aprior$params[2]
        )

      # extract the gradient vector
      grad.prior[1] <- attributes(rst.aprior)$gradient
    }

    # extract the gradient vector when the difficulty parameter prior is used
    if (use.bprior) {
      rst.bprior <-
        logprior_deriv(
          val = b, is.aprior = FALSE, D = NULL, dist = bprior$dist,
          par.1 = bprior$params[1], par.2 = bprior$params[2]
        )

      # extract the gradient vector
      grad.prior[2] <- attributes(rst.bprior)$gradient
    }
  }

  # (4) 3PLM
  if (!fix.g & mod == "3PLM") {
    # assign a, b, g parameters
    a <- item_par[1]
    b <- item_par[2]
    g <- item_par[3]

    # compute the probabilities of correct and incorrect
    # (reuse the cached P matrix when supplied by the optim factory)
    p <- if (is.null(p_cache)) drm(theta = theta, a = a, b = b, g = g, D = D) else p_cache

    # compute the component values
    g_1 <- 1 - g
    p_g <- p - g
    r_p <- r_i - (f_i * p)
    theta_b <- theta - b
    r_p_p <- r_p / p
    u1 <- p_g * r_p_p
    u2 <- D / g_1
    u3 <- a * u2

    # compute the gradients of a and bs parameters
    gr_a <- -u2 * sum(theta_b * u1)
    gr_b <- u3 * sum(u1)
    gr_g <- -(1 / g_1) * sum(r_p_p)

    # combine all gradients into a vector
    grad <- c(gr_a, gr_b, gr_g)

    # create a prior gradient vector
    grad.prior <- rep(0, n.par)

    # extract the gradient vector when the slope parameter prior is used
    if (use.aprior) {
      rst.aprior <-
        logprior_deriv(
          val = a, is.aprior = TRUE, D = D, dist = aprior$dist,
          par.1 = aprior$params[1], par.2 = aprior$params[2]
        )

      # extract the gradient vector
      grad.prior[1] <- attributes(rst.aprior)$gradient
    }

    # extract the gradient vector when the difficulty parameter prior is used
    if (use.bprior) {
      rst.bprior <-
        logprior_deriv(
          val = b, is.aprior = FALSE, D = NULL, dist = bprior$dist,
          par.1 = bprior$params[1], par.2 = bprior$params[2]
        )

      # extract the gradient vector
      grad.prior[2] <- attributes(rst.bprior)$gradient
    }

    # extract the gradient vector when the guessing parameter prior is used
    if (use.gprior) {
      rst.gprior <-
        logprior_deriv(
          val = g, is.aprior = FALSE, D = NULL, dist = gprior$dist,
          par.1 = gprior$params[1], par.2 = gprior$params[2]
        )

      # extract the gradient vector
      grad.prior[3] <- attributes(rst.gprior)$gradient
    }
  }

  # (5) 3PLM: the guessing parameters are fixed to be specified value
  if (fix.g & mod == "3PLM") {
    # assign a, b, g parameters
    a <- item_par[1]
    b <- item_par[2]
    g <- g.val

    # compute the probabilities of correct and incorrect
    # (reuse the cached P matrix when supplied by the optim factory)
    p <- if (is.null(p_cache)) drm(theta = theta, a = a, b = b, g = g, D = D) else p_cache

    # compute the component values
    g_1 <- 1 - g
    p_g <- p - g
    r_p <- r_i - (f_i * p)
    theta_b <- theta - b
    u1 <- p_g * r_p / p
    u2 <- D / g_1
    u3 <- a * u2

    # compute the gradients of a and bs parameters
    gr_a <- -u2 * sum(theta_b * u1)
    gr_b <- u3 * sum(u1)

    # combine all gradients into a vector
    grad <- c(gr_a, gr_b)

    # create a prior gradient vector
    grad.prior <- rep(0, n.par)

    # extract the gradient vector when the slope parameter prior is used
    if (use.aprior) {
      rst.aprior <-
        logprior_deriv(
          val = a, is.aprior = TRUE, D = D, dist = aprior$dist,
          par.1 = aprior$params[1], par.2 = aprior$params[2]
        )

      # extract the gradient vector
      grad.prior[1] <- attributes(rst.aprior)$gradient
    }

    # extract the gradient vector when the difficulty parameter prior is used
    if (use.bprior) {
      rst.bprior <-
        logprior_deriv(
          val = b, is.aprior = FALSE, D = NULL, dist = bprior$dist,
          par.1 = bprior$params[1], par.2 = bprior$params[2]
        )

      # extract the gradient vector
      grad.prior[2] <- attributes(rst.bprior)$gradient
    }
  }

  # add the prior gradient vector
  grad <- grad + grad.prior

  # return results
  grad
}


# This function analytically computes the gradient vector of the negative log-likelihood (plus the negative log
# priors) of polytomous item parameters
grad_item_prm <- function(item_par, r_i, theta, pr.mod, D = 1, fix.a = FALSE, a.val = 1,
                          aprior = list(dist = "lnorm", params = c(0, 0.5)),
                          bprior = list(dist = "norm", params = c(0.0, 1.0)),
                          use.aprior = FALSE,
                          use.bprior = FALSE,
                          prob_cache = NULL) {
  # `prob_cache`, when non-NULL, is produced by make_prm_optim_fns() and
  # carries the pre-computed probability objects for the current item_par:
  #   GRM  -> $allPst (boundary probs from drm()), $P (clipped category probs)
  #   GPCM -> $theta_d, $numer, $denom, $P
  # Using the cache skips the most expensive matrix operations (drm() for GRM;
  # exp(cumsum) + rowsums for GPCM) that loglike_prm already paid for.
  # count the number of item parameters to be estimated
  n.par <- length(item_par)

  ## -------------------------------------------------------------------------
  # compute the gradients
  # (1) GRM
  if (pr.mod == "GRM") {
    # assign a, b parameters
    a <- item_par[1]
    d <- item_par[-1]

    # check the number of step parameters
    m <- length(d)

    # calculate all the probabilities greater than equal to each threshold;
    # reuse cached allPst when supplied
    if (is.null(prob_cache)) {
      allPst <- drm(theta = theta, a = rep(a, m), b = d, g = 0, D = D)
      allQst <- 1 - allPst[, , drop = FALSE]
      P <- (cbind(1, allPst) - cbind(allPst, 0))
      P[P > 9999999999e-10] <- 9999999999e-10
      P[P < 1e-10] <- 1e-10
    } else {
      allPst <- prob_cache$allPst
      allQst <- 1 - allPst[, , drop = FALSE]
      P <- prob_cache$P  # already clipped by the factory
    }

    # compute the component values to get gradients
    Da <- D * a
    pq_st <- allPst * allQst
    q_p_st <- allQst - allPst
    w1 <- D * (-Rfast::Outer(x = d, y = theta, oper = "-"))
    w2 <- w1 * pq_st
    w3 <- cbind(0, w2) - cbind(w2, 0)
    frac_rp <- r_i / P
    w4 <- -Rfast::coldiffs(frac_rp)

    # compute the gradients of a and bs parameters
    gr_a <- -sum(frac_rp * w3)
    gr_b <- -Da * Rfast::colsums(pq_st * w4)

    # combine all gradients into a vector
    grad <- c(gr_a, gr_b)

    # create a prior gradient vector
    grad.prior <- rep(0, n.par)

    # extract the gradient vector when the slope parameter prior is used
    if (use.aprior) {
      rst.aprior <-
        logprior_deriv(
          val = a, is.aprior = TRUE, D = D, dist = aprior$dist,
          par.1 = aprior$params[1], par.2 = aprior$params[2]
        )

      # extract the gradient vector
      grad.prior[1] <- attributes(rst.aprior)$gradient
    }

    # extract the gradient vector when the difficulty parameter prior is used
    if (use.bprior) {
      rst.bprior <-
        logprior_deriv(
          val = d, is.aprior = FALSE, D = NULL, dist = bprior$dist,
          par.1 = bprior$params[1], par.2 = bprior$params[2]
        )

      # extract the gradient vector
      grad.prior[-1] <- attributes(rst.bprior)$gradient
    }
  }

  # (2) GPCM and PCM
  if (pr.mod == "GPCM") {
    if (!fix.a) {
      # For GPCM
      # assign a parameter
      a <- item_par[1]

      # include zero for the step parameter of the first category
      d <- c(0, item_par[-1])

      # check the number of step parameters
      m <- length(d) - 1

      # calculate category probabilities (theta_d, numer, denom, P);
      # reuse cached values when supplied; the entire exp(cumsum) block
      # that loglike_prm already computed is skipped on cache hit.
      # Da is a scalar (D * a) and is not cached: negligible to recompute.
      Da <- D * a
      if (is.null(prob_cache)) {
        theta_d <- (Rfast::Outer(x = theta, y = d, oper = "-"))
        z <- Da * theta_d
        cumsum_z <- t(Rfast::colCumSums(z))
        # subtract the row maximum so that exp() and the powers of the denominator stay finite
        if (any(abs(cumsum_z) > 100)) {
          cumsum_z <- cumsum_z - Rfast::rowMaxs(cumsum_z, value = TRUE)
        }
        numer <- exp(cumsum_z) # numerator
        denom <- Rfast::rowsums(numer)
        P <- (numer / denom)
        # keep the probabilities above zero so that their logarithms are finite
        P[P < .Machine$double.xmin] <- .Machine$double.xmin
      } else {
        theta_d <- prob_cache$theta_d
        numer <- prob_cache$numer
        denom <- prob_cache$denom
        P <- prob_cache$P
      }

      # compute the component values to get a gradient vector
      dsmat <- array(0, c((m + 1), (m + 1)))
      dsmat[-1, -1][upper.tri(x = dsmat[-1, -1], diag = TRUE)] <- 1
      frac_rp <- r_i / P
      w1 <- D * theta_d
      w1cum <- t(Rfast::colCumSums(w1))
      denom2 <- denom^2
      d1a_denom <- Rfast::rowsums(numer * w1cum)
      deriv_Pa <- (numer * (w1cum * denom - d1a_denom)) / denom2
      d1b_denom <- tcrossprod(x = (-Da * numer), y = dsmat)
      denom_vec <- cbind(denom)
      DaDenom_vec <- Da * denom_vec

      # compute the gradients of a and bs parameters
      gr_a <- -sum(frac_rp * deriv_Pa)
      gr_b <- c()
      for (k in 1:m) {
        gr_b[k] <- -sum(frac_rp * (-numer * (DaDenom_vec %*% dsmat[k + 1, , drop = FALSE] + d1b_denom[, k + 1]) / denom2))
      }

      # combine all gradients into a vector
      grad <- c(gr_a, gr_b)

      # create a prior gradient vector
      grad.prior <- rep(0, n.par)

      # extract the gradient vector when the slope parameter prior is used
      if (use.aprior) {
        rst.aprior <-
          logprior_deriv(
            val = a, is.aprior = TRUE, D = D, dist = aprior$dist,
            par.1 = aprior$params[1], par.2 = aprior$params[2]
          )

        # extract the gradient vector
        grad.prior[1] <- attributes(rst.aprior)$gradient
      }

      # extract the gradient vector when the difficulty parameter prior is used
      if (use.bprior) {
        rst.bprior <-
          logprior_deriv(
            val = d[-1], is.aprior = FALSE, D = NULL, dist = bprior$dist,
            par.1 = bprior$params[1], par.2 = bprior$params[2]
          )

        # extract the gradient vector
        grad.prior[-1] <- attributes(rst.bprior)$gradient
      }
    } else {
      # for PCM
      # assign a parameter
      a <- a.val

      # include zero for the step parameter of the first category
      d <- c(0, item_par)

      # check the number of step parameters
      m <- length(d) - 1

      # calculate category probabilities (same cache strategy as GPCM)
      Da <- D * a
      if (is.null(prob_cache)) {
        theta_d <- (Rfast::Outer(x = theta, y = d, oper = "-"))
        z <- Da * theta_d
        cumsum_z <- t(Rfast::colCumSums(z))
        # subtract the row maximum so that exp() and the powers of the denominator stay finite
        if (any(abs(cumsum_z) > 100)) {
          cumsum_z <- cumsum_z - Rfast::rowMaxs(cumsum_z, value = TRUE)
        }
        numer <- exp(cumsum_z) # numerator
        denom <- Rfast::rowsums(numer)
        P <- (numer / denom)
        # keep the probabilities above zero so that their logarithms are finite
        P[P < .Machine$double.xmin] <- .Machine$double.xmin
      } else {
        theta_d <- prob_cache$theta_d
        numer <- prob_cache$numer
        denom <- prob_cache$denom
        P <- prob_cache$P
      }

      # compute the component values to get a gradient vector
      dsmat <- array(0, c((m + 1), (m + 1)))
      dsmat[-1, -1][upper.tri(x = dsmat[-1, -1], diag = TRUE)] <- 1
      frac_rp <- r_i / P
      denom2 <- denom^2
      d1b_denom <- tcrossprod(x = (-Da * numer), y = dsmat)
      denom_vec <- cbind(denom)
      DaDenom_vec <- Da * denom_vec

      # compute the gradients of a and bs parameters
      gr_b <- c()
      for (k in 1:m) {
        gr_b[k] <- -sum(frac_rp * (-numer * (DaDenom_vec %*% dsmat[k + 1, , drop = FALSE] + d1b_denom[, k + 1]) / denom2))
      }

      # combine all gradients into a vector
      grad <- gr_b

      # create a prior gradient vector
      grad.prior <- rep(0, n.par)

      # extract the gradient vector when the difficulty parameter prior is used
      if (use.bprior) {
        rst.bprior <-
          logprior_deriv(
            val = d[-1], is.aprior = FALSE, D = NULL, dist = bprior$dist,
            par.1 = bprior$params[1], par.2 = bprior$params[2]
          )

        # extract the gradient vector
        grad.prior <- attributes(rst.bprior)$gradient
      }
    }
  }

  # add the prior gradient vector
  grad <- grad + grad.prior

  # return results
  grad
}



# This function returns, for each row of theta, the gradient of the negative log-likelihood of an item with respect
# to its parameters; info_xpd() uses it to build the cross-product information matrix
grad_llike <- function(item_par, f_i, r_i, theta, n.theta, mod = c("1PLM", "2PLM", "3PLM", "GRM", "GPCM"),
                       D = 1, fix.a.gpcm = FALSE, fix.g = FALSE, a.val.1pl = 1, a.val.gpcm = 1,
                       g.val = .2, n.1PLM = NULL) {
  # for dichotomous models
  if (mod %in% c("1PLM", "2PLM", "3PLM")) {
    if (!is.null(n.1PLM)) {
      # 1PLM: when the item slope parameters are not constrained to be equal across all items
      # compute the gradient vectors
      grad <- grad_item_drm_se(
        item_par = item_par, f_i = f_i, r_i = r_i, theta = theta,
        mod = mod, D = D, n.1PLM = n.1PLM
      )
    } else {
      # for all other dichotomous models
      # compute the gradient vectors
      grad <- grad_item_drm_se(
        item_par = item_par, f_i = f_i, r_i = r_i, theta = theta,
        mod = mod, D = D, fix.g = fix.g,
        a.val = a.val.1pl, g.val = g.val, n.1PLM = NULL
      )
    }
  } else {
    # for polytomous models
    # compute the gradient vectors
    grad <- grad_item_prm_se(
      item_par = item_par, r_i = r_i, theta = theta, n.theta = n.theta,
      pr.mod = mod, D = D, fix.a = fix.a.gpcm, a.val = a.val.gpcm
    )
  }

  # return results
  grad
}


# This function analytically computes a matrix of negative log-likelihood gradients of dichotomous item parameters,
# with one row per theta value; info_xpd() uses it for the cross-product information matrix
grad_item_drm_se <- function(item_par, f_i, r_i, theta, mod = c("1PLM", "2PLM", "3PLM", "DRM"), D = 1,
                             fix.g = TRUE, a.val = 1, g.val = .2, n.1PLM = NULL) {
  # count the number of item parameters to be estimated
  n.par <- length(item_par)

  # (1) 1PLM: the slope parameters are constrained to be equal across the 1PLM items
  if (!is.null(n.1PLM)) {
    # make vectors of a and b parameters for all 1PLM items
    a <- rep(item_par[1], n.1PLM)
    b <- item_par[-1]

    # compute the probabilities of correct and incorrect
    p <- drm(theta = theta, a = a, b = b, g = 0, D = D)

    # compute the component values
    r_p <- r_i - f_i * p
    theta_b <- -Rfast::Outer(x = b, y = theta, oper = "-")

    # compute the gradients of a and bs parameters
    gr_a <- -D * Rfast::rowsums(theta_b * r_p)
    gr_b <- (D * a[1]) * r_p

    # create a matrix of gradients across all examinees
    grad <- cbind(gr_a, gr_b)
  }

  # (2) 1PLM with a fixed slope
  if (is.null(n.1PLM) && mod == "1PLM") {
    # assign a and b parameters
    a <- a.val
    b <- item_par

    # compute the probabilities of correct and incorrect
    p <- drm(theta = theta, a = a, b = b, g = 0, D = D)

    # compute the component values
    r_p <- r_i - f_i * p

    # compute the gradients of b parameter
    gr_b <- (D * a) * r_p

    # create a matrix of gradients across all examinees
    grad <- cbind(gr_b)
  }

  # (3) 2PLM
  if (mod == "2PLM") {
    # assign the item parameters
    a <- item_par[1]
    b <- item_par[2]

    # compute the probabilities of correct and incorrect
    p <- drm(theta = theta, a = a, b = b, g = 0, D = D)

    # compute the component values
    Da <- D * a
    r_p <- r_i - f_i * p
    theta_b <- theta - b

    # compute the gradients of a and b parameters
    gr_a <- (-D * theta_b) * r_p
    gr_b <- Da * r_p

    # create a matrix of gradients across all examinees
    grad <- cbind(gr_a, gr_b)
  }

  # (4) 3PLM
  if (!fix.g & mod == "3PLM") {
    # assign the item parameters
    a <- item_par[1]
    b <- item_par[2]
    g <- item_par[3]

    # compute the probabilities of correct and incorrect
    p <- drm(theta = theta, a = a, b = b, g = g, D = D)

    # compute the component values
    g_1 <- 1 - g
    p_g <- p - g
    r_p <- r_i - f_i * p
    theta_b <- theta - b
    r_p_p <- r_p / p
    u1 <- p_g * r_p_p
    u2 <- D / g_1
    u3 <- a * u2

    # compute the gradients of a, b, and g parameters
    gr_a <- (-u2 * theta_b) * u1
    gr_b <- u3 * u1
    gr_g <- -(1 / g_1) * r_p_p

    # create a matrix of gradients across all examinees
    grad <- cbind(gr_a, gr_b, gr_g)
  }

  # (5) 3PLM: the guessing parameters are fixed to be specified value
  if (fix.g & mod == "3PLM") {
    # assign the item parameters
    a <- item_par[1]
    b <- item_par[2]
    g <- g.val

    # compute the probabilities of correct and incorrect
    p <- drm(theta = theta, a = a, b = b, g = g, D = D)

    # compute the component values
    g_1 <- 1 - g
    p_g <- p - g
    r_p <- r_i - f_i * p
    theta_b <- theta - b
    u1 <- p_g * r_p / p
    u2 <- D / g_1
    u3 <- a * u2

    # compute the gradients of a and b parameters
    gr_a <- (-u2 * theta_b) * u1
    gr_b <- u3 * u1

    # create a matrix of gradients across all examinees
    grad <- cbind(gr_a, gr_b)
  }

  # return results
  grad
}

# This function analytically computes a matrix of negative log-likelihood gradients of polytomous item parameters,
# with one row per theta value; info_xpd() uses it for the cross-product information matrix
grad_item_prm_se <- function(item_par, r_i, theta, n.theta, pr.mod, D = 1, fix.a = FALSE, a.val = 1) {
  # count the number of item parameters to be estimated
  n.par <- length(item_par)

  ## -------------------------------------------------------------------------
  # compute the gradients
  # (1) GRM
  if (pr.mod == "GRM") {
    # assign the slope and threshold parameters
    a <- item_par[1]
    d <- item_par[-1]

    # check the number of step parameters
    m <- length(d)

    # calculate all the probabilities greater than equal to each threshold
    allPst <- drm(theta = theta, a = rep(a, m), b = d, g = 0, D = D)
    allQst <- 1 - allPst[, , drop = FALSE]

    # calculate category probabilities
    P <- (cbind(1, allPst) - cbind(allPst, 0))
    P[P > 9999999999e-10] <- 9999999999e-10
    P[P < 1e-10] <- 1e-10

    # compute the component values to get gradients
    Da <- D * a
    pq_st <- allPst * allQst
    q_p_st <- allQst - allPst
    w1 <- D * (-Rfast::Outer(x = d, y = theta, oper = "-"))
    w2 <- w1 * pq_st
    w3 <- cbind(0, w2) - cbind(w2, 0)
    frac_rp <- r_i / P
    w4 <- -Rfast::coldiffs(frac_rp)

    # compute the gradients of a and bs parameters
    gr_a <- -Rfast::rowsums(frac_rp * w3)
    gr_b <- (-Da * pq_st) * w4

    # create a matrix of gradients across all examinees
    grad <- cbind(gr_a, gr_b)
  }

  # (2) GPCM and PCM
  if (pr.mod == "GPCM") {
    if (!fix.a) {
      # For GPCM
      # assign the slope parameter
      a <- item_par[1]

      # include zero for the step parameter of the first category
      d <- c(0, item_par[-1])

      # check the number of step parameters
      m <- length(d) - 1

      # calculate category probabilities
      Da <- D * a
      theta_d <- (Rfast::Outer(x = theta, y = d, oper = "-"))
      z <- Da * theta_d
      cumsum_z <- t(Rfast::colCumSums(z))
      # subtract the row maximum so that exp() and the powers of the denominator stay finite
      if (any(abs(cumsum_z) > 100)) {
        cumsum_z <- cumsum_z - Rfast::rowMaxs(cumsum_z, value = TRUE)
      }
      numer <- exp(cumsum_z) # numerator
      denom <- Rfast::rowsums(numer)
      P <- (numer / denom)
      # keep the probabilities above zero so that their logarithms are finite
      P[P < .Machine$double.xmin] <- .Machine$double.xmin

      # compute the component values to get a gradient vector
      dsmat <- array(0, c((m + 1), (m + 1)))
      dsmat[-1, -1][upper.tri(x = dsmat[-1, -1], diag = TRUE)] <- 1
      frac_rp <- r_i / P
      w1 <- D * theta_d
      w1cum <- t(Rfast::colCumSums(w1))
      denom2 <- denom^2
      d1a_denom <- Rfast::rowsums(numer * w1cum)
      deriv_Pa <- (numer * (w1cum * denom - d1a_denom)) / denom2
      d1b_denom <- tcrossprod(x = (-Da * numer), y = dsmat)
      denom_vec <- cbind(denom)
      DaDenom_vec <- Da * denom_vec

      # compute the gradients of a and bs parameters
      gr_a <- -Rfast::rowsums(frac_rp * deriv_Pa)
      gr_b <- array(0, c(n.theta, m))
      for (k in 1:m) {
        gr_b[, k] <-
          -Rfast::rowsums(frac_rp *
            (-numer * (DaDenom_vec %*%
              dsmat[k + 1, , drop = FALSE] +
              d1b_denom[, k + 1]) / denom2))
      }

      # create a matrix of gradients across all examinees
      grad <- cbind(gr_a, gr_b)
    } else {
      # for PCM
      # assign the fixed slope parameter
      a <- a.val

      # include zero for the step parameter of the first category
      d <- c(0, item_par)

      # check the number of step parameters
      m <- length(d) - 1

      # calculate category probabilities
      Da <- D * a
      theta_d <- (Rfast::Outer(x = theta, y = d, oper = "-"))
      z <- Da * theta_d
      cumsum_z <- t(Rfast::colCumSums(z))
      # subtract the row maximum so that exp() and the powers of the denominator stay finite
      if (any(abs(cumsum_z) > 100)) {
        cumsum_z <- cumsum_z - Rfast::rowMaxs(cumsum_z, value = TRUE)
      }
      numer <- exp(cumsum_z) # numerator
      denom <- Rfast::rowsums(numer)
      P <- (numer / denom)
      # keep the probabilities above zero so that their logarithms are finite
      P[P < .Machine$double.xmin] <- .Machine$double.xmin

      # compute the component values to get a gradient vector
      dsmat <- array(0, c((m + 1), (m + 1)))
      dsmat[-1, -1][upper.tri(x = dsmat[-1, -1], diag = TRUE)] <- 1
      frac_rp <- r_i / P
      denom2 <- denom^2
      d1b_denom <- tcrossprod(x = (-Da * numer), y = dsmat)
      denom_vec <- cbind(denom)
      DaDenom_vec <- Da * denom_vec

      # compute the gradients of a and bs parameters
      gr_b <- array(0, c(n.theta, m))
      for (k in 1:m) {
        gr_b[, k] <-
          -Rfast::rowsums(frac_rp *
            (-numer * (DaDenom_vec %*%
              dsmat[k + 1, , drop = FALSE] +
              d1b_denom[, k + 1]) / denom2))
      }

      # create a matrix of gradients across all examinees
      grad <- gr_b
    }
  }

  # return results
  grad
}
