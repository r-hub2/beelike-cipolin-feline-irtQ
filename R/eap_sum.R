# Compute EAP summed scores
#
# This function computes the expected a posteriori (EAP) estimate of ability
# for each sum score (Thissen et al., 1995; Thissen & Orlando, 2001). The
# estimate is the mean of the posterior distribution of theta given the sum
# score, and the standard error is the posterior standard deviation. The
# likelihood of each sum score at each quadrature point is computed with the
# Lord-Wingersky recursion in lwrc().
#
# Arguments
# x: item metadata (see shape_df()).
# data: a matrix of item responses (rows: examinees, columns: items). Missing
#   responses are replaced with 0.
# norm.prior: mean and standard deviation of the normal prior used to generate
#   Gaussian quadrature points and weights when weights = NULL.
# nquad: number of Gaussian quadrature points. Default is 41.
# weights: a two-column matrix or data frame of quadrature points and weights.
# D: scaling constant (1.702 approximates the normal ogive). Default is 1.
#
# Value
# A list with est.par (sum.score, est.theta, and se.theta for each examinee)
# and score.table (the same columns for every possible sum score).
#
#' @importFrom Rfast colsums rowsums
eap_sum <- function(x, data, norm.prior = c(0, 1), nquad = 41, weights = NULL, D = 1) {
  
  ## ------------------------------------------------------------------------------------------------
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

  # generate quad nodes and weights
  if (is.null(weights)) {
    weights <- gen.weight(n = nquad, dist = "norm", mu = norm.prior[1], sigma = norm.prior[2])
  } else {
    weights <- data.frame(weights)
  }

  # compute the sum score likelihoods with the Lord-Wingersky recursion
  lkhd <- lwrc(x = x, theta = weights[, 1], prob = NULL, D = D)

  # compute the posterior mean and standard deviation for each sum score
  ss.prob <- c(lkhd %*% weights[, 2])
  post <- t((t(lkhd) * weights[, 2])) / ss.prob
  tr_post <- t(post)
  eap.est <- Rfast::colsums(tr_post * weights[, 1])
  eap.est2 <- Rfast::colsums(tr_post * weights[, 1]^2)
  se.est <- sqrt(eap.est2 - eap.est^2)

  # assign the EAP summed scores to each examinee
  obs.score <- 0:(length(eap.est) - 1)
  names(eap.est) <- obs.score
  names(se.est) <- obs.score
  sumScore <- Rfast::rowsums(data)
  est_score <- eap.est[as.character(sumScore)]
  est_se <- se.est[as.character(sumScore)]
  score_table <-
    data.frame(
      sum.score = obs.score,
      est.theta = eap.est, se.theta = se.est,
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

  # return results
  rst <- list(
    est.par = est.par,
    score.table = score_table
  )
  rst
}
