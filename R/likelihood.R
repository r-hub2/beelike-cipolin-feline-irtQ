# This function computes the matrix of likelihoods L(theta_q | response_i)
# for each examinee at every quadrature point.  It is the core ingredient
# of the E-step in the MMLE-EM algorithm: posterior() weights this matrix
# by the prior to obtain the posterior density of ability per examinee.
#
# Two test components contribute additively in log-space:
#   (i)  dichotomous items (DRM): 1PLM / 2PLM / 3PLM via drm()
#   (ii) polytomous items (PRM):  GRM / GPCM via prm()
#
# Args:
#   elm_item : list returned by breakdown(); contains $pars (one row per item),
#              $cats, and $model
#   idx.drm  : integer vector of row indices in elm_item$pars for DRM items;
#              NULL if the test has no DRM items
#   idx.prm  : integer vector of row indices for PRM items; NULL if none
#   data_drm : sparse cbind(data_drm_q, data_drm_p) where each item
#              contributes a 'q' (incorrect = 1) and 'p' (correct = 1)
#              indicator column; rows are examinees (NA responses appear
#              as zero in both columns so they contribute nothing)
#   data_prm : sparse cbind of category-indicator columns across PRM items
#              (one column per (item, score-category) pair)
#   theta    : numeric vector of quadrature points
#   D        : scaling constant (1 for the logistic metric; 1.702 approximates the normal ogive)
#
# Returns:
#   list with two elements:
#     $L      : nstd x ntheta numeric matrix of likelihoods.  For an examinee
#               whose likelihoods sum to less than exp(-700) over the grid (a
#               long test whose likelihood underflows in exp()), the row is
#               divided by exp(row maximum of the log-likelihood) so that its
#               largest element is 1.
#     $offset : numeric vector of length nstd holding the log of that divisor
#               (0 for every row that is not rescaled).  The log-likelihood of
#               examinee i is log(L[i, ]) + offset[i]; posterior() is
#               unaffected because the divisor cancels within a row.
#' @importFrom Matrix tcrossprod
likelihood <- function(elm_item, idx.drm = NULL, idx.prm = NULL,
                       data_drm = NULL, data_prm = NULL, theta, D = 1) {

  ## --------------------------------------------------------------------
  ## (1) Log-likelihood contribution from DRM (dichotomous) items
  ## --------------------------------------------------------------------
  if (!is.null(idx.drm)) {

    # P(correct response) at each (theta, item) pair: ntheta x n.drm matrix
    # produced by drm() vectorized over theta and items simultaneously
    ps <- drm(
      theta = theta, a = elm_item$pars[idx.drm, 1], b = elm_item$pars[idx.drm, 2],
      g = elm_item$pars[idx.drm, 3], D = D
    )

    # P(incorrect response) = 1 - P(correct), same shape
    qs <- 1 - ps

    # element-wise log of [Q | P] stacked column-wise (ntheta x 2*n.drm);
    # the column order (qs first, then ps) must match data_drm, which is
    # cbind(data_drm_q, data_drm_p) in the caller
    log_qps <- log(cbind(qs, ps))

    # log-likelihood for DRM items per (examinee, quadrature point):
    # tcrossprod(data_drm, log_qps) computes data_drm %*% t(log_qps),
    # producing an nstd x ntheta matrix where entry (i, q) is the sum of
    # log p(observed response_i | theta_q) over all DRM items.  data_drm
    # is a sparse Matrix, so this exploits zero entries (incl. missing
    # responses encoded as 0/0)
    llike_drm <- Matrix::tcrossprod(x = data_drm, y = log_qps)

  } else {

    # no DRM items: contribute zero to the additive log-likelihood
    # (integer 0 broadcasts harmlessly with the PRM term below)
    llike_drm <- 0L
  }

  ## --------------------------------------------------------------------
  ## (2) Log-likelihood contribution from PRM (polytomous) items
  ## --------------------------------------------------------------------
  if (!is.null(idx.prm)) {

    # number of polytomous items in the test
    n.prm <- length(idx.prm)

    # accumulator for per-item category-probability matrices; preallocated
    # to avoid quadratic copying when growing inside the loop
    prob.prm <- vector("list", n.prm)

    # category probabilities cannot be computed in a single vectorized
    # call because items can have different numbers of score categories
    # (and different models, GRM vs GPCM), so the loop runs item by item
    for (k in 1:n.prm) {

      # item-k parameter row, dropping NA padding columns (unused threshold
      # slots that exist when this item has fewer categories than the
      # widest item in the test bank)
      par.tmp <- stats::na.exclude(elm_item$pars[idx.prm[k], ])

      # P(category j | theta_q) at each quadrature point: ntheta x cats[k]
      # matrix; prm() dispatches on pr.model to GRM or GPCM internally
      prob.prm[[k]] <- prm(
        theta = theta, a = par.tmp[1], d = par.tmp[-1],
        D = D, pr.model = elm_item$model[idx.prm[k]]
      )
    }

    # concatenate per-item category-probability matrices column-wise into a
    # single ntheta x sum(cats) matrix; column order matches data_prm
    prob.prm <- do.call(what = "cbind", prob.prm)

    # element-wise log; grm() bounds the probabilities to [1e-10, 1 - 1e-10]
    # and gpcm() bounds them below by the smallest positive normalized double,
    # so log() never produces -Inf or NaN here
    log_prob.prm <- log(prob.prm)

    # log-likelihood for PRM items per (examinee, quadrature point):
    # tcrossprod(data_prm, log_prob.prm) sums log P(observed category)
    # across all PRM items.  data_prm is the sparse category-indicator
    # matrix (one column per (item, category)) so only the chosen
    # categories contribute per examinee
    llike_prm <- Matrix::tcrossprod(x = data_prm, y = log_prob.prm)

  } else {

    # no PRM items: contribute zero
    llike_prm <- 0L
  }

  ## --------------------------------------------------------------------
  ## (3) Combine the two components and convert to the likelihood scale
  ## --------------------------------------------------------------------

  # total log-likelihood = sum of DRM and PRM contributions; one of them
  # may be the integer 0 placeholder, which broadcasts elementwise.
  # base::as.matrix() densifies the sparse-Matrix output of tcrossprod()
  # so that exp() returns a plain numeric matrix usable downstream.
  llike <- base::as.matrix(llike_drm + llike_prm)

  # convert to the likelihood scale
  L <- exp(llike)

  # log-scale offset of each examinee (0 unless the row is rescaled below)
  offset <- numeric(nrow(L))

  # find the rows whose likelihoods underflow: their sum over the grid is below
  # exp(-700), so every element is close to the double underflow limit or below it
  idx.uf <- which(rowSums(L) < exp(-700))

  # recompute those rows only after subtracting the row maximum of the
  # log-likelihood, and keep that maximum as the offset; all other rows keep
  # exactly the same likelihoods as exp(llike)
  if (length(idx.uf) > 0L) {
    row.max <- apply(llike[idx.uf, , drop = FALSE], 1, max)

    # a row with no finite log-likelihood cannot be rescaled
    idx.ok <- is.finite(row.max)
    idx.uf <- idx.uf[idx.ok]
    row.max <- row.max[idx.ok]
    L[idx.uf, ] <- exp(llike[idx.uf, , drop = FALSE] - row.max)
    offset[idx.uf] <- row.max
  }

  # return the likelihood matrix and the log-scale offsets; callers that
  # need the marginal log-likelihood add the offsets back (Estep, est_irt,
  # est_mg), while posterior() can use $L directly (Estep, pcd2)
  list(L = L, offset = offset)
}
