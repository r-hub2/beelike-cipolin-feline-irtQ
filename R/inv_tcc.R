# Find the ability estimate for each sum score by inverting the test
# characteristic curve (TCC). This function modifies irt.eq.tse() from the
# SNSequate package (Gonzalez, 2014). The root for each sum score is found with
# bisection(), a modified version of bisection() from the cmna package
# (Howard, 2017).
#
# References
# Gonzalez, J. (2014). SNSequate: Standard and nonstandard statistical models
#   and methods for test equating. Journal of Statistical Software, 59(7), 1-30.
# Howard, J. P. (2017). Computational methods for numerical analysis with R.
#   New York: Chapman and Hall/CRC.
#
#' @importFrom Rfast colsums
inv_tcc <- function(x, data, D = 1, intpol = TRUE, range.tcc = c(-7, 7), tol = 1e-4, max.it = 500) {
  
  # transform a data set to matrix
  data <- data.matrix(data)
  
  # check missing data
  # replace NAs with 0
  na.lg <- is.na(data)
  if (any(na.lg)) {
    data[na.lg] <- 0
    memo <- "Missing responses are replaced with 0."
    warning(memo, call. = FALSE)
  }

  # find the ability estimates for all possible observed sum scores
  tab <- inv_tcc_nr(
    x = x, D = D, intpol = intpol, range.tcc = range.tcc, tol = tol,
    max.it = max.it
  )

  # extract the observed sum scores
  obs.score <- tab$sum.score

  # extract the ability estimates without names
  thetas <- unname(tab$est.theta)

  # flag the observed sum scores that have an ability estimate
  obs2theta.lg <- !is.na(thetas)

  # a vector of SEs corresponding to all thetas
  se.theta <- rep(NA, length(obs.score))

  # a vector of theta excluding NAs
  thetas.nona <- thetas[obs2theta.lg]

  # compute the conditional sum score distribution at each estimated theta
  # with the Lord-Wingersky recursion
  if (length(thetas.nona) > 1) {
    lkhd <- lwrc(x = x, theta = thetas.nona, D = D)

    # keep only the scores that have an estimate and rescale each column to sum to one
    if (!all(obs2theta.lg)) {
      # drop the rows of the scores without an ability estimate
      lkhd <- lkhd[obs2theta.lg, , drop = FALSE]

      # rescale the conditional probabilities over the remaining scores
      lkhd <- t(t(lkhd) / Rfast::colsums(lkhd))
    }

    # calculate the standard error of ability estimates
    mu <- Rfast::colsums(lkhd * thetas.nona)
    se.theta[obs2theta.lg] <- sqrt(Rfast::colsums(lkhd * thetas.nona^2) - mu^2)
  }

  # assign the ability estimates and SEs to examinees
  names(thetas) <- obs.score
  names(se.theta) <- obs.score
  sumScore <- Rfast::rowsums(data)
  est_score <- thetas[as.character(sumScore)]
  est_se <- se.theta[as.character(sumScore)]
  score_table <-
    data.frame(
      sum.score = obs.score,
      est.theta = thetas, se.theta = se.theta,
      stringsAsFactors = FALSE
    )
  rownames(score_table) <- NULL

  # create a data frame for the estimated scores
  est.par <- data.frame(
    sum.score = sumScore,
    est.theta = est_score,
    se.theta = est_se
  )
  rownames(est.par) <- NULL

  # return the results
  rst <- list(
    est.par = est.par,
    score.table = score_table
  )
  rst
}


# This function finds the ability estimates (thetas) corresponding
# to the possible observed scores using the test characteristic curve (TCC).
# No response data set needs to be provided.
# This function is used by inv_tcc() and reval_mst().
inv_tcc_nr <- function(x, D = 1, intpol = TRUE, range.tcc = c(-7, 7),
                       tol = 1e-4, max.it = 500) {

  # break down the item metadata into several elements
  elm_item <- breakdown(x)

  ## extract score category information
  cats <- elm_item$cats

  # classify the items into DRM and PRM item groups
  idx.item <- idxfinder(elm_item)
  idx.drm <- idx.item$idx.drm

  # set all available observed sum scores
  obs.score <- 0:sum(cats - 1)

  # possible maximum total observed sum score
  max.obs <- max(obs.score)

  # a set of indexes of observed sum scores whose corresponding
  # theta scores will be found
  obs2theta.lg <- logical(length(obs.score))

  # admissible sum scores for the root search
  # when DRM items exist
  # \sum c_j < \tau < K, eq (6.19), Kolen & Brennan (2004, p.176)
  if (!is.null(idx.drm)) {
    g_sum <- sum(elm_item$pars[idx.drm, 3])
    if (max.obs == 1) {
      idx.score <- NULL
    } else if (ceiling(g_sum) == g_sum) {
      # when g_sum is an integer, the smallest sum score with a root
      # is g_sum + 1
      idx.score <- ((ceiling(g_sum) + 1):(max.obs - 1) + 1)
    } else {
      idx.score <- (ceiling(g_sum):(max.obs - 1) + 1)
    }
  } else {
    g_sum <- 0
    idx.score <- ((g_sum + 1):(max.obs - 1) + 1)
  }

  # update the score index
  obs2theta.lg[idx.score] <- TRUE

  # a vector of thetas corresponding to the observed sum scores
  thetas <- rep(NA, length(obs.score))

  # a function to estimate a theta corresponding to the observed sum score
  # fun(theta) = tau - \sum p(theta), eq (6.21), Kolen & Brennan (2004, p.177)
  f.o2t <- function(theta, tau) {
    tau - trace(elm_item = elm_item, theta = theta, D = D, tcc = TRUE)$tcc
  }

  # compute tcc values at discrete theta nodes
  theta.nodes <- seq(-20, 20, 0.01)
  tcc.vals <- trace(elm_item = elm_item, theta = theta.nodes, D = D, tcc = TRUE)$tcc

  # find the thetas for all admissible sum scores
  th4obs <- vapply(
    X = obs.score[obs2theta.lg], FUN.VALUE = numeric(1), USE.NAMES = FALSE,
    FUN = function(x) {
      loc.node <- which(diff(sign(x - tcc.vals)) != 0)

      # stop with a clear message when the TCC does not reach the sum score
      if (length(loc.node) == 0L) {
        stop(
          "The TCC does not reach the sum score ", x,
          " within theta in [-20, 20]; INV.TCC scoring is not available.",
          call. = FALSE
        )
      }
      bd <- theta.nodes[c(loc.node, loc.node + 1)]
      bisection(
        .fun = f.o2t, tau = x, lb = bd[1], ub = bd[2], tol = tol,
        max.it = max.it
      )$root
    }
  )
  thetas[obs2theta.lg] <- th4obs

  # if intpol = TRUE, interpolate linearly for sum scores less than or equal
  # to g_sum and assign range.tcc[2] to the maximum sum score
  if (intpol) {
    if (max.obs == 1) {
      # the test has a single dichotomous item
      thetas <- range.tcc
      obs2theta.lg <- c(TRUE, TRUE)
    } else if (range.tcc[1] > th4obs[1]) {
      memo <- paste0(
        "A lower bound of theta must be less than ", round(th4obs[1], 3), "\n",
        "Set a lower bound below this value in the 'range.tcc' argument."
      )
      warning(memo, call. = FALSE)
    } else if (range.tcc[2] < utils::tail(th4obs, 1)) {
      memo <- paste0(
        "An upper bound of theta must be greater than ", round(utils::tail(th4obs, 1), 3), "\n",
        "Set an upper bound above this value in the 'range.tcc' argument."
      )
      warning(memo, call. = FALSE)
    } else {
      if (g_sum > 0) {
        lessg.score <- 0:(obs.score[idx.score[1]] - 1)
        slope <- obs.score[obs2theta.lg][1] / (th4obs[1] - range.tcc[1])
        intercept <- -slope * range.tcc[1]
        lessg.theta <- (lessg.score - intercept) / slope
        thetas[!obs2theta.lg] <- c(lessg.theta, range.tcc[2])
        thetas[is.nan(thetas)] <- range.tcc[1]
      } else {
        thetas[!obs2theta.lg] <- range.tcc
      }
      obs2theta.lg[!obs2theta.lg] <- TRUE
    }
  }

  # Table of the estimated thetas corresponding to each observed score
  names(thetas) <- obs.score
  score_table <-
    data.frame(
      sum.score = obs.score,
      est.theta = thetas,
      stringsAsFactors = FALSE
    )
  rownames(score_table) <- NULL

  # return the results
  score_table
}
