# This function is an internal function used in the 'est_item' function.
estimation1 <- function(f_i, r_i, s_i, theta, mod = c("1PLM", "2PLM", "3PLM", "GRM", "GPCM"), score.cat, D = 1, nstd,
                        fix.a.1pl = TRUE, fix.a.gpcm = FALSE, fix.g = FALSE, a.val.1pl = 1, a.val.gpcm = 1, g.val = .2, n.1PLM = NULL,
                        aprior = list(dist = "lnorm", params = c(0, 0.5)),
                        bprior = list(dist = "norm", params = c(0.0, 1.0)),
                        gprior = list(dist = "beta", params = c(5, 16)),
                        use.aprior = FALSE, use.bprior = FALSE, use.gprior = TRUE,
                        control, startval = NULL, lower, upper) {
  # build the cached (objective, gradient, hessian) trio once per call
  # for DRM models; see make_drm_optim_fns() for the cache mechanics.
  # n.1PLM only matters for the (!fix.a & mod=="1PLM") branch, so it is NULL otherwise
  if (mod %in% c("1PLM", "2PLM", "3PLM")) {
    drm_fns <- make_drm_optim_fns(
      f_i = f_i, r_i = r_i, s_i = s_i, theta = theta, mod = mod, D = D,
      fix.a = fix.a.1pl, fix.g = fix.g, a.val = a.val.1pl, g.val = g.val,
      n.1PLM = if (!fix.a.1pl & mod == "1PLM") n.1PLM else NULL,
      aprior = aprior, bprior = bprior, gprior = gprior,
      use.aprior = use.aprior, use.bprior = use.bprior, use.gprior = use.gprior
    )
  }

  # estimation of item parameter & standard error
  if (!fix.a.1pl & mod == "1PLM") {
    # estimate the item parameters
    est <-
      suppressWarnings(
        tryCatch(
          {
            stats::nlminb(startval,
              objective = drm_fns$objective,
              gradient = drm_fns$gradient,
              hessian = drm_fns$hessian,
              control = control, lower = lower, upper = upper
            )
          },
          error = function(e) {
            NULL
          }
        )
      )
    if (is.null(est) || est$convergence > 0L) {
      # if error or non-convergence, only use the gradient
      est <- stats::nlminb(startval,
        objective = drm_fns$objective,
        gradient = drm_fns$gradient,
        control = control, lower = lower, upper = upper
      )
    }

    # estimate the standard error of estimates
    # (single call, so no p_cache is passed)
    hess <- hess_item_drm(est$par,
      f_i = f_i, r_i = r_i, theta = theta, mod = mod, D = D,
      fix.a = fix.a.1pl, fix.g = fix.g, a.val = a.val.1pl, g.val = g.val, n.1PLM = n.1PLM,
      aprior = aprior, bprior = bprior, gprior = gprior,
      use.aprior = use.aprior, use.bprior = use.bprior, use.gprior = use.gprior,
      adjust = FALSE
    )
  } else {
    # when the slope parameters are not constrained to be equal across items
    if (mod %in% c("1PLM", "2PLM", "3PLM")) {
      # initial estimation to find better starting values
      # (gradient-only nlminb sweeps with three step.max settings; the
      # cache still helps because nlminb evaluates objective then
      # gradient at each accepted iterate)
      tmp_est <- vector("list", 3)
      tmp_est[[1]] <- stats::nlminb(startval,
        objective = drm_fns$objective,
        gradient = drm_fns$gradient,
        control = list(eval.max = 100, iter.max = 50, trace = 0, step.min = 0.01, step.max = 1), lower = lower, upper = upper
      )
      tmp_est[[2]] <- stats::nlminb(startval,
        objective = drm_fns$objective,
        gradient = drm_fns$gradient,
        control = list(eval.max = 100, iter.max = 50, trace = 0, step.min = 0.01, step.max = 2), lower = lower, upper = upper
      )
      tmp_est[[3]] <- stats::nlminb(startval,
        objective = drm_fns$objective,
        gradient = drm_fns$gradient,
        control = list(eval.max = 100, iter.max = 50, trace = 0, step.min = 0.01, step.max = 3), lower = lower, upper = upper
      )
      tmp_num <- which.min(c(tmp_est[[1]]$objective, tmp_est[[2]]$objective, tmp_est[[3]]$objective))
      startval <- tmp_est[[tmp_num]]$par

      # estimate the item parameters
      est <-
        suppressWarnings(
          tryCatch(
            {
              stats::nlminb(startval,
                objective = drm_fns$objective,
                gradient = drm_fns$gradient,
                hessian = drm_fns$hessian,
                control = control, lower = lower, upper = upper
              )
            },
            error = function(e) {
              NULL
            }
          )
        )

      # second estimation with second alternative values when the first calibration fails
      if (is.null(est) || est$convergence > 0L) {
        tmp_num <- order(c(tmp_est[[1]]$objective, tmp_est[[2]]$objective, tmp_est[[3]]$objective))[2]
        startval <- tmp_est[[tmp_num]]$par
        est <-
          suppressWarnings(
            tryCatch(
              {
                stats::nlminb(startval,
                  objective = drm_fns$objective,
                  gradient = drm_fns$gradient,
                  hessian = drm_fns$hessian,
                  control = control, lower = lower, upper = upper
                )
              },
              error = function(e) {
                NULL
              }
            )
          )
      }

      # third estimation by only using the gradient when the second calibration fails
      if (is.null(est) || est$convergence > 0L) {
        tmp_num <- which.min(c(tmp_est[[1]]$objective, tmp_est[[2]]$objective, tmp_est[[3]]$objective))
        startval <- tmp_est[[tmp_num]]$par

        # if error or non-convergence, only use the gradient
        est <- stats::nlminb(startval,
          objective = drm_fns$objective,
          gradient = drm_fns$gradient,
          control = control, lower = lower, upper = upper
        )
      }

      # estimate the standard error of estimates
      # (single call after optimization; cache not needed)
      hess <- hess_item_drm(est$par,
        f_i = f_i, r_i = r_i, theta = theta, mod = mod, D = D,
        fix.a = fix.a.1pl, fix.g = fix.g, a.val = a.val.1pl, g.val = g.val, n.1PLM = NULL,
        aprior = aprior, bprior = bprior, gprior = gprior,
        use.aprior = use.aprior, use.bprior = use.bprior, use.gprior = use.gprior,
        adjust = FALSE
      )
    } else {
      # keep the starting values for a refit with ordered GRM thresholds
      startval.in <- startval

      # build the cached (objective, gradient, hessian) trio for PRM
      prm_fns <- make_prm_optim_fns(
        r_i = r_i, theta = theta, pr.mod = mod, D = D, nstd = nstd,
        fix.a = fix.a.gpcm, a.val = a.val.gpcm,
        aprior = aprior, bprior = bprior,
        use.aprior = use.aprior, use.bprior = use.bprior
      )

      # initial estimation to find better starting values
      tmp_est <- stats::nlminb(startval,
        objective = prm_fns$objective,
        gradient = prm_fns$gradient,
        control = list(eval.max = 50, iter.max = 30, step.min = 0.01, step.max = 1, trace = 0), lower = lower, upper = upper
      )
      startval <- tmp_est$par

      # estimate the item parameters
      est <-
        suppressWarnings(
          tryCatch(
            {
              stats::nlminb(startval,
                objective = prm_fns$objective,
                gradient = prm_fns$gradient,
                hessian = prm_fns$hessian,
                control = control, lower = lower, upper = upper
              )
            },
            error = function(e) {
              NULL
            }
          )
        )

      # second estimation by only using the gradient when the first calibration fails
      if (is.null(est) || est$convergence > 0L) {
        # if error or non-convergence, only use the gradient
        est <- stats::nlminb(startval,
          objective = prm_fns$objective,
          gradient = prm_fns$gradient,
          control = control, lower = lower, upper = upper
        )
      }

      # refit with ordered thresholds when the GRM thresholds are not in increasing order
      if (mod == "GRM" && isTRUE(is.unsorted(est$par[-1], strictly = TRUE))) {
        est <- est_grm_ordered(
          startval = startval.in, prm_fns = prm_fns, control = control,
          lower = lower, upper = upper
        )
      }

      # estimate the standard error of estimates
      # (single SE call; cache not needed)
      hess <- hess_item_prm(est$par,
        r_i = r_i, theta = theta, pr.mod = mod, D = D,
        nstd = nstd, fix.a = fix.a.gpcm, a.val = a.val.gpcm,
        aprior = aprior, bprior = bprior, use.aprior = use.aprior, use.bprior = use.bprior,
        adjust = FALSE
      )
    }
  }

  # compute the variance-covariance matrix of the item parameter estimates, and
  # check if the hessian matrix can be inversed
  cov_mat <- suppressWarnings(tryCatch(
    {
      solve(hess, tol = 1e-200)
    },
    error = function(e) {
      NULL
    }
  ))

  # compute the standard errors of item parameter estimates
  # also, if the hessian matrix does not have an inverse matrix, then create a matrix with NaNs
  if (is.null(cov_mat)) {
    se <- rep(99999, length(diag(hess)))
    cov_mat <- array(NaN, c(length(startval), length(startval)))
  } else {
    se <- suppressWarnings(sqrt(diag(cov_mat)))
  }

  # prevent showing NaN values of standard errors
  if (any(is.nan(se))) {
    se[is.nan(se)] <- 99999
  }

  # set an upper bound of standard error
  se <- ifelse(se > 99999, 99999, se)

  # negative log-likelihood at the estimates without the log prior densities
  if (mod %in% c("1PLM", "2PLM", "3PLM")) {
    objective.data <- loglike_drm(est$par,
      r_i = r_i, s_i = s_i, theta = theta, mod = mod, D = D,
      fix.a = fix.a.1pl, fix.g = fix.g, a.val = a.val.1pl, g.val = g.val, n.1PLM = n.1PLM,
      use.aprior = FALSE, use.bprior = FALSE, use.gprior = FALSE
    )
  } else {
    objective.data <- loglike_prm(est$par,
      r_i = r_i, theta = theta, pr.mod = mod, D = D, fix.a = fix.a.gpcm, a.val = a.val.gpcm,
      use.aprior = FALSE, use.bprior = FALSE
    )
  }

  # return results
  rst <- list(
    pars = est$par, se = se, covariance = cov_mat, convergence = est$convergence, objective = est$objective,
    objective.data = objective.data
  )
  rst
}


# This function is an internal function used in the 'est_irt' and 'est_mg' functions.
estimation2 <- function(f_i, r_i, s_i, quadpt, mod = c("1PLM", "2PLM", "3PLM", "GRM", "GPCM"), score.cat, D = 1, n.quad,
                        fix.a.1pl = TRUE, fix.a.gpcm = FALSE, fix.g = FALSE, a.val.1pl = 1, a.val.gpcm = 1, g.val = .2, n.1PLM = NULL,
                        aprior = list(dist = "lnorm", params = c(0, 0.5)),
                        bprior = list(dist = "norm", params = c(0.0, 1.0)),
                        gprior = list(dist = "beta", params = c(5, 16)),
                        use.aprior = FALSE, use.bprior = FALSE, use.gprior = TRUE,
                        control, startval = NULL, lower, upper, iter = NULL) {
  # keep the values at the start of the cycle for the case in which the last search fails
  startval.prev <- startval

  # build the cached (objective, gradient, hessian) trio once per call
  # for DRM models; see make_drm_optim_fns() for the cache mechanics.
  if (mod %in% c("1PLM", "2PLM", "3PLM")) {
    drm_fns <- make_drm_optim_fns(
      f_i = f_i, r_i = r_i, s_i = s_i, theta = quadpt,
      mod = mod, D = D,
      fix.a = fix.a.1pl, fix.g = fix.g, a.val = a.val.1pl, g.val = g.val,
      n.1PLM = if (!fix.a.1pl & mod == "1PLM") n.1PLM else NULL,
      aprior = aprior, bprior = bprior, gprior = gprior,
      use.aprior = use.aprior, use.bprior = use.bprior, use.gprior = use.gprior
    )
  }

  # estimation of item parameters
  if (!fix.a.1pl & mod == "1PLM") {
    # when the item slope parameters are constrained to be equal across all 1PLM items
    # estimate the item parameters
    est <-
      suppressWarnings(
        tryCatch(
          {
            stats::nlminb(startval,
              objective = drm_fns$objective,
              gradient = drm_fns$gradient,
              hessian = drm_fns$hessian,
              control = control, lower = lower, upper = upper
            )
          },
          error = function(e) {
            NULL
          }
        )
      )
    if (is.null(est) || est$convergence > 0L) {
      # if error or non-convergence, only use the gradient
      est <- tryCatch(
        stats::nlminb(startval,
          objective = drm_fns$objective,
          gradient = drm_fns$gradient,
          control = control, lower = lower, upper = upper
        ),
        error = function(e) {
          NULL
        }
      )
      # keep the previous values and mark the item as not converged when the search fails
      if (is.null(est)) est <- list(par = startval.prev, convergence = 1L, objective = NA_real_)
    }
  } else {
    # when the slope parameters are not constrained to be equal across items
    if (mod %in% c("1PLM", "2PLM", "3PLM")) {
      # in the first two EM cycles, run short gradient-only searches to find better starting values
      if (iter < 3) {
        # a failed search is kept as NULL
        tmp_est <- lapply(c(1, 2, 3), function(step.max) {
          tryCatch(
            stats::nlminb(startval,
              objective = drm_fns$objective,
              gradient = drm_fns$gradient,
              control = list(
                eval.max = 40 - (iter * 10), iter.max = 25 - (iter * 5),
                trace = 0, step.min = 0.01, step.max = step.max
              ),
              lower = lower, upper = upper
            ),
            error = function(e) NULL
          )
        })
        # use the best successful short search, or keep the starting values when all of them fail
        tmp_ok <- !vapply(tmp_est, is.null, logical(1))
        tmp_num <- which.min(vapply(tmp_est[tmp_ok], function(e) e$objective, numeric(1)))
        if (length(tmp_num) == 1L) startval <- tmp_est[tmp_ok][[tmp_num]]$par
      }

      # estimate the item parameters
      est <- tryCatch(
        {
          stats::nlminb(startval,
            objective = drm_fns$objective,
            gradient = drm_fns$gradient,
            hessian = drm_fns$hessian,
            control = control, lower = lower, upper = upper
          )
        },
        error = function(e) {
          NULL
        }
      )
      if (is.null(est) || est$convergence > 0L) {
        # if error or non-convergence, only use the gradient
        est <- tryCatch(
          stats::nlminb(startval,
            objective = drm_fns$objective,
            gradient = drm_fns$gradient,
            control = control, lower = lower, upper = upper
          ),
          error = function(e) {
            NULL
          }
        )
        # keep the previous values and mark the item as not converged when the search fails
        if (is.null(est)) est <- list(par = startval.prev, convergence = 1L, objective = NA_real_)
      }
    } else {
      # keep the starting values for a refit with ordered GRM thresholds
      startval.in <- startval

      # build the cached (objective, gradient, hessian) trio for PRM
      prm_fns <- make_prm_optim_fns(
        r_i = r_i, theta = quadpt, pr.mod = mod, D = D, nstd = n.quad,
        fix.a = fix.a.gpcm, a.val = a.val.gpcm,
        aprior = aprior, bprior = bprior,
        use.aprior = use.aprior, use.bprior = use.bprior
      )

      # in the first two EM cycles, run short gradient-only searches to find better starting values
      if (iter < 3) {
        # a failed search is kept as NULL
        tmp_est <- lapply(c(1, 2, 3), function(step.max) {
          tryCatch(
            stats::nlminb(startval,
              objective = prm_fns$objective,
              gradient = prm_fns$gradient,
              control = list(
                eval.max = 40 - (iter * 10), iter.max = 25 - (iter * 5),
                step.min = 0.01, step.max = step.max, trace = 0
              ),
              lower = lower, upper = upper
            ),
            error = function(e) NULL
          )
        })
        # use the best successful short search, or keep the starting values when all of them fail
        tmp_ok <- !vapply(tmp_est, is.null, logical(1))
        tmp_num <- which.min(vapply(tmp_est[tmp_ok], function(e) e$objective, numeric(1)))
        if (length(tmp_num) == 1L) startval <- tmp_est[tmp_ok][[tmp_num]]$par
      }

      # estimate the item parameters
      est <- tryCatch(
        {
          stats::nlminb(startval,
            objective = prm_fns$objective,
            gradient = prm_fns$gradient,
            hessian = prm_fns$hessian,
            control = control, lower = lower, upper = upper
          )
        },
        error = function(e) {
          NULL
        }
      )
      if (is.null(est) || est$convergence > 0L) {
        # if error or non-convergence, only use the gradient
        est <- tryCatch(
          stats::nlminb(startval,
            objective = prm_fns$objective,
            gradient = prm_fns$gradient,
            control = control, lower = lower, upper = upper
          ),
          error = function(e) {
            NULL
          }
        )
        # keep the previous values and mark the item as not converged when the search fails
        if (is.null(est)) est <- list(par = startval.prev, convergence = 1L, objective = NA_real_)
      }

      # refit with ordered thresholds when the GRM thresholds are not in increasing order
      if (mod == "GRM" && isTRUE(is.unsorted(est$par[-1], strictly = TRUE))) {
        est <- est_grm_ordered(
          startval = startval.in, prm_fns = prm_fns, control = control,
          lower = lower, upper = upper
        )
      }
    }
  }

  # return results
  rst <- list(pars = est$par, convergence = est$convergence, objective = est$objective)

  rst
}


# This function re-estimates the parameters of a GRM item with its thresholds
# kept in increasing order.  The thresholds are written as b_1 and the logs of
# the gaps b_k - b_(k-1), so every point visited by nlminb has ordered
# thresholds.  It is used only when the unconstrained fit returns thresholds
# that are not ordered: grm() truncates the negative category probabilities of
# such a fit at 1e-10, so the category probabilities no longer sum to one and
# the objective is not a valid likelihood there.
est_grm_ordered <- function(startval, prm_fns, control, lower, upper) {
  # number of thresholds
  m <- length(startval) - 1L

  # smallest gap and largest log gap allowed between adjacent thresholds
  min.gap <- 1e-8
  max.lgap <- 20

  # starting values on the transformed scale (sorted thresholds, positive gaps)
  b0 <- sort(startval[-1])
  phi0 <- c(startval[1], b0[1], log(pmax(diff(b0), min.gap)))

  # map the transformed parameters (a, b_1, log gaps) back to (a, b_1, ..., b_m)
  to_par <- function(phi) c(phi[1], cumsum(c(phi[2], exp(phi[-(1:2)]))))

  # Jacobian d(a, b_1, ..., b_m) / d(a, b_1, log gaps)
  jacob <- function(phi) {
    J <- diag(m + 1L)

    # every threshold moves one-to-one with b_1
    J[2:(m + 1L), 2] <- 1

    # threshold k depends on the log gap j through exp(log gap j) when k >= j
    for (j in seq_len(m)[-1]) J[(j + 1L):(m + 1L), j + 1L] <- exp(phi[j + 1L])
    J
  }

  # objective, gradient, and hessian on the transformed scale (chain rule)
  obj <- function(phi) prm_fns$objective(to_par(phi))
  grd <- function(phi) as.vector(crossprod(jacob(phi), prm_fns$gradient(to_par(phi))))
  hss <- function(phi) {
    J <- jacob(phi)
    gb <- prm_fns$gradient(to_par(phi))[-1]

    # second-derivative term of the transformation: exp(log gap j) * sum_{k >= j} g_k
    curv <- c(0, 0, exp(phi[-(1:2)]) * rev(cumsum(rev(gb)))[-1])
    crossprod(J, prm_fns$hessian(to_par(phi)) %*% J) + diag(curv, m + 1L)
  }

  # bounds on the transformed scale (the slope keeps its own bounds)
  lower.phi <- c(lower[1], -Inf, rep(log(min.gap), m - 1L))
  upper.phi <- c(upper[1], Inf, rep(max.lgap, m - 1L))

  # estimate with the hessian first, and only with the gradient when that fails
  est <- tryCatch(
    {
      stats::nlminb(phi0,
        objective = obj, gradient = grd, hessian = hss,
        control = control, lower = lower.phi, upper = upper.phi
      )
    },
    error = function(e) {
      NULL
    }
  )
  if (is.null(est) || est$convergence > 0L) {
    est <- stats::nlminb(phi0,
      objective = obj, gradient = grd,
      control = control, lower = lower.phi, upper = upper.phi
    )
  }

  # return the estimates on the original scale
  list(par = to_par(est$par), convergence = est$convergence, objective = est$objective)
}
