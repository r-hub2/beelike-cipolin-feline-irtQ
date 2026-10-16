#' Asymptotic Variance-Covariance Matrices of Item Parameter Estimates
#'
#' This function computes the analytical asymptotic variance-covariance matrices
#' of item parameter estimates for dichotomous and polytomous IRT models without
#' examinee response data. For each item, the matrix is the inverse of the
#' expected (Fisher) information of the item parameters, which is integrated
#' over an assumed ability distribution and multiplied by the calibration sample
#' size (Thissen & Wainer, 1982; Li & Lissitz, 2004). The square roots of the
#' diagonal elements are the asymptotic standard errors of the maximum
#' likelihood estimates.
#'
#' @param x A data frame containing item metadata (item IDs, numbers of score
#'   categories, IRT models, and item parameters). The item parameters are
#'   treated as the true values at which the information is evaluated; for 1PLM
#'   items and for PCM items specified in `pcm.loc`, the slope in `x` is held
#'   fixed. See [irtQ::est_irt()] or [irtQ::simdat()] for the format of the item
#'   metadata, which can be created with [irtQ::shape_df()].
#' @param D A scaling constant used in IRT models to make the logistic function
#'   closely approximate the normal ogive function. A value of 1.702 is commonly
#'   used for this purpose. Default is 1.
#' @param nstd An integer, or a vector of integers with one value per item in
#'   `x`, giving the number of examinees used to calibrate the items. The
#'   information of each item is multiplied by its sample size, so the standard
#'   errors are proportional to 1 / sqrt(`nstd`). Default is 1000. See
#'   **Details**.
#' @param pcm.loc A vector of integers giving the positions (rows of `x`) of
#'   GPCM items that were calibrated with the partial credit model (PCM). For
#'   these items, the slope is held at the value in `x` and the
#'   variance-covariance matrices are computed only for the threshold
#'   parameters. Positions of items that are not GPCM items are ignored. Default
#'   is `NULL`. See **Details**.
#' @param norm.prior A numeric vector of length two specifying the mean and
#'   standard deviation of the normal ability distribution over which the
#'   information is integrated. Used to generate the Gaussian quadrature points
#'   and weights when `weights = NULL`. Default is `c(0, 1)`.
#' @param nquad An integer specifying the number of Gaussian quadrature points
#'   generated from the normal distribution given in `norm.prior`. Used only
#'   when `weights = NULL`. Default is 41.
#' @param weights An optional two-column data frame or matrix of quadrature
#'   points (first column) and weights (second column) that defines the ability
#'   distribution. The weights should sum to 1, because they are not normalized.
#'   If `weights` is supplied, `norm.prior` and `nquad` are ignored. Such a data
#'   frame can be created with [irtQ::gen.weight()]. Default is `NULL`.
#'
#' @details
#' For each item, the asymptotic variance-covariance matrix is the inverse of
#' the expected (Fisher) information matrix of the item parameters, which is
#' integrated over the ability distribution with quadrature and multiplied by
#' the sample size (Thissen & Wainer, 1982; Li & Lissitz, 2004). For an item
#' with \eqn{K} score categories, the element of the information matrix for a
#' pair of item parameters \eqn{\xi_s} and \eqn{\xi_t} is
#'
#' \deqn{I(\xi_s, \xi_t) = N \sum_{q=1}^{Q} A(X_q) \sum_{k=1}^{K}
#'   \frac{1}{P_k(X_q)} \frac{\partial P_k(X_q)}{\partial \xi_s}
#'   \frac{\partial P_k(X_q)}{\partial \xi_t},}{I(xi_s, xi_t) = N sum_q A(X_q)
#'   sum_k [1 / P_k(X_q)] [d P_k(X_q) / d xi_s] [d P_k(X_q) / d xi_t],}
#'
#' where \eqn{N} is the sample size (`nstd`), \eqn{X_q} and \eqn{A(X_q)} are
#' the \eqn{Q} quadrature points and weights that represent the ability
#' distribution (`norm.prior` and `nquad`, or `weights`), and \eqn{P_k(X_q)} is
#' the probability of a response in category \eqn{k} at \eqn{X_q}. This is the
#' information of a categorical (multinomial) response. For a dichotomous item,
#' it reduces to the familiar form in which each term is divided by
#' \eqn{P(X_q) [1 - P(X_q)]} (Li & Lissitz, 2004, Equation A4).
#'
#' For polytomous items, the information above differs from Equation A8 of Li
#' and Lissitz (2004). Their Equation A6 writes the likelihood of a response
#' (and Equation A7 its logarithm) as
#' \eqn{\prod_{k} P_k^{u_k} (1 - P_k)^{1 - u_k}}{prod_k
#' P_k^u_k (1 - P_k)^(1 - u_k)}, where \eqn{u_k}
#' is 1 for the observed category and 0 otherwise. The \eqn{K} category
#' indicators are thus treated as independent Bernoulli variables, and
#' Equation A8 accordingly divides each term by \eqn{P_k (1 - P_k)} instead of
#' \eqn{P_k}. However, the indicators of one response sum to 1, so the response
#' follows a multinomial distribution with likelihood
#' \eqn{\prod_{k} P_k^{u_k}}{prod_k P_k^u_k},
#' whose information divides each term by \eqn{P_k} only. Because
#' \eqn{1 - P_k < 1}, Equation A8 overstates the information and understates the
#' standard errors. For example, for an item with two categories, Equation A8
#' gives exactly twice the information of their Equation A4 for dichotomous
#' items. This function uses the multinomial information.
#'
#' The multinomial information was checked by a simulation that matches the
#' setting of this function: abilities drawn from N(0, 1) and treated as known,
#' maximum likelihood estimation, a sample size of 1,000, and 1,000
#' replications for GRM and GPCM items with three to five categories. The
#' standard errors from the multinomial information agreed with the standard
#' deviations of the estimates across replications (ratios of 0.96 to 1.08),
#' whereas those from Equation A8 were about 10 to 23 percent smaller (ratios
#' of 0.77 to 0.90). For a dichotomous (2PLM) item, the two forms give the
#' same results.
#'
#' The matrices are derived under the assumptions that the IRT model fits the
#' data, that the item parameters are estimated by maximum likelihood, and that
#' the ability distribution is known (Li & Lissitz, 2004). Because the
#' abilities are treated as known, the uncertainty in the abilities is ignored,
#' and prior distributions on the item parameters are not taken into account.
#' The resulting standard errors are therefore generally considered lower bounds
#' of the standard errors of an actual calibration (Thissen & Wainer, 1982).
#' They may serve as approximations for evaluating the precision of item
#' parameter estimates when empirical standard errors are not reported.
#'
#' If the item parameters provided in the `x` argument were calibrated using
#' different sample sizes, a corresponding vector of sample sizes must be
#' specified via the `nstd` argument. For example, suppose you wish to compute
#' the variance-covariance matrices of five 3PLM items that were calibrated
#' using 500, 600, 1,000, 2,000, and 700 examinees, respectively. In this case,
#' set `nstd = c(500, 600, 1000, 2000, 700)`.
#'
#' Since the item metadata allows only `"GPCM"` to denote both the partial
#' credit model (PCM) and the generalized partial credit model (GPCM), PCM items
#' must be explicitly identified using the `pcm.loc` argument. For these items,
#' the slope is not treated as an estimated parameter: it is held at the value
#' in `x`, and only the threshold parameters enter the information matrix. For
#' instance, if the last two of five polytomous items were calibrated using the
#' PCM, specify `pcm.loc = c(4, 5)`. A GRM or GPCM item with two score
#' categories is treated as a 2PLM item, and a PCM item with two score
#' categories as a 1PLM item with a fixed slope.
#'
#' @return A named list with the following two components:
#'
#' - `cov`: A named list of variance-covariance matrices, one per item. Each
#' matrix contains only the estimated parameters of the item, and its row and
#' column names are the parameter columns of `x`: `par.2` (difficulty) for the
#' 1PLM; `par.1` and `par.2` (slope and difficulty) for the 2PLM; `par.1` to
#' `par.3` (slope, difficulty, and guessing) for the 3PLM; `par.1` to `par.K`
#' (slope and K - 1 thresholds) for GRM and GPCM items with K score categories;
#' and `par.2` to `par.K` (thresholds only) for PCM items specified in
#' `pcm.loc`. Covariances between parameters of different items are not
#' computed. If the information matrix of an item is singular, a warning is
#' issued and the matrix is filled with `NA`.
#'
#' - `se`: A named list of vectors containing the asymptotic standard errors
#' (ASEs) of the item parameter estimates, computed as the square roots of the
#' diagonal elements of each corresponding variance-covariance matrix in `cov`.
#'
#' The names of the list elements in both `cov` and `se` correspond to the item
#' identifiers (e.g., item names or labels) as given in the first column of the
#' input `x`.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @references Li, Y. H., & Lissitz, R. W. (2004). Applications of the
#'   analytically derived asymptotic standard errors of item response theory
#'   item parameter estimates. *Journal of Educational Measurement, 41*(2),
#'   85-117. \doi{10.1111/j.1745-3984.2004.tb01109.x}.
#'
#'   Thissen, D., & Wainer, H. (1982). Some standard errors in item response
#'   theory. *Psychometrika, 47*(4), 397-412. \doi{10.1007/BF02293705}.
#'
#' @seealso [irtQ::est_irt()], [irtQ::simdat()], [irtQ::shape_df()],
#'   [irtQ::gen.weight()]
#'
#' @export
#'
#' @examples
#' # Example using a "-prm.txt" file exported from flexMIRT
#' flex_prm <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
#'
#' # Select the first two dichotomous items and the last polytomous item
#' x <- bring.flexmirt(file = flex_prm, "par")$Group1$full_df[c(1:2, 55), ]
#'
#' # Compute the variance-covariance matrices assuming a sample size of 2,000
#' covirt(x, D = 1, nstd = 2000, norm.prior = c(0, 1), nquad = 41)
#'
covirt <- function(x, D = 1, nstd = 1000, pcm.loc = NULL, norm.prior = c(0, 1), nquad = 41, weights = NULL) {
  # confirm and correct all item metadata information
  x <- confirm_df(x)

  # break down the item metadata into several components
  elm_list <-
    purrr::map(.x = 1:nrow(x), ~ {
      breakdown(x[.x, ])
    })

  # specify the locations of PCM items
  pcm.a.logic <- rep(FALSE, nrow(x))
  if (!is.null(pcm.loc)) {
    pcm.a.logic[pcm.loc] <- TRUE
  }

  # generate quadrature points and weights
  if (is.null(weights)) {
    wts <- gen.weight(n = nquad, dist = "norm", mu = norm.prior[1], sigma = norm.prior[2])
  } else {
    wts <- data.frame(weights)
  }

  # sample size
  if (length(nstd) == 1L) nstd <- rep(nstd, nrow(x))
  if (length(nstd) != nrow(x)) {
    stop("Sample size must be an integer value or a vector with a length of total items.", call. = FALSE)
  }

  # compute the variance-covariance matrix of each item
  cov_par <- purrr::pmap(
    .l = list(x = elm_list, y = nstd, z = pcm.a.logic),
    .f = function(x, y, z) {
      cov_mat(elm_item = x, D = D, nstd = y, wts = wts, pcm.a = z)
    }
  )
  se_par <- purrr::map(cov_par, .f = function(x) sqrt(diag(x)))
  names(cov_par) <- x[, 1]
  names(se_par) <- x[, 1]

  # return results
  rst <- list(cov = cov_par, se = se_par)
  rst
}


# This function computes the analytical variance-covariance matrix of the parameter estimates of an item
cov_mat <- function(elm_item, D = 1, nstd = 1000, wts, pcm.a = FALSE) {
  info_mat <- integrand(elm_item, theta = wts[, 1], dens = wts[, 2], D = D, pcm.a = pcm.a)

  # invert the information matrix multiplied by the sample size
  rst <- tryCatch(
    {
      solve(info_mat * nstd)
    },
    error = function(e) {
      NULL
    }
  )

  # warn and return NA with the dimnames of the information matrix when the matrix is singular
  if (is.null(rst)) {
    warning(
      "The information matrix of item ", elm_item$id, " is singular, so NA is returned.",
      call. = FALSE
    )
    rst <- matrix(NA_real_, nrow = nrow(info_mat), ncol = ncol(info_mat), dimnames = dimnames(info_mat))
  }

  rst
}


# This function computes the expected information matrix of the parameters of an item
# by quadrature over the ability distribution
#' @importFrom stats na.exclude
#' @importFrom Rfast rowprods
integrand <- function(elm_item, theta, dens, D = 1, pcm.a = FALSE) {
  # assign items into DRM and PRM groups
  idx.all <- idxfinder(elm_item)
  idx.drm <- idx.all$idx.drm
  idx.prm <- idx.all$idx.prm

  # For a DRM item
  if (!is.null(idx.drm)) {
    # extract required information (par.3 is NA for 1PLM/2PLM; substitute 0)
    item_par <- c(
      elm_item$pars[idx.drm, 1],
      elm_item$pars[idx.drm, 2],
      elm_item$pars[idx.drm, 3]
    )
    item_par[is.na(item_par)] <- 0
    cats <- elm_item$cats[idx.drm]
    model <- elm_item$model[idx.drm]

    # treat a GRM or GPCM item with two categories as a 2PLM item (1PLM item for a PCM item)
    if (model %in% c("GRM", "GPCM")) {
      model <- if (model == "GPCM" && pcm.a) "1PLM" else "2PLM"
    }

    # compute a product of the probabilities p and q
    pq <- Rfast::rowprods(trace(elm_item = elm_item, theta = theta, D = D, tcc = FALSE)$prob.cats[[1]])

    # replace zero products with 1e-20 to avoid division by zero
    pq[pq == 0L] <- 1e-20

    # create the gradient for score category probability
    funList <- equation_scocat(model = model, cats = NULL, hessian = FALSE, type = "item")

    # create a list containing the arguments to be used in the equation function
    args.pars <- list()
    for (i in 1:3) {
      args.pars[[i]] <- item_par[i]
    }
    args.list <- args.pars
    args.list$theta <- theta
    args.list$D <- D

    # compute the score category information
    # select a function for each score category
    fun.tmp <- funList[[1]]

    # implement the function for each score category
    tmp <- do.call("fun.tmp", args.list, envir = environment())

    # compute the kernel of integration
    p.d1 <- attributes(tmp)$gradient
    info.list <- purrr::map(
      .x = 1:length(theta),
      .f = function(i) (1 / pq[i]) * outer(p.d1[i, ], p.d1[i, ]) * dens[i]
    )
    info.mat <- Reduce(f = "+", x = info.list)
  }

  # For a PRM item
  if (!is.null(idx.prm)) {
    # extract required information
    cats <- elm_item$cats[idx.prm]
    model <- elm_item$model[idx.prm]
    item_par <- stats::na.exclude(elm_item$pars[idx.prm, ])

    # compute the score category probabilities; each category term of the information is divided by P_k
    ps <- trace(elm_item = elm_item, theta = theta, D = D, tcc = FALSE)$prob.cats[[1]]

    # replace zero probabilities with 1e-20 to avoid division by zero
    ps[ps == 0L] <- 1e-20

    # create the gradient for score category probability
    funList <- equation_scocat(model = model, cats = cats, fix.a.gpcm = pcm.a, hessian = FALSE, type = "item")

    # create a list containing the arguments to be used in the equation function
    args.pars <- list()
    for (i in 1:cats) {
      args.pars[[i]] <- item_par[i]
    }
    args.list <- args.pars
    args.list$theta <- theta
    args.list$D <- D

    # compute the score category information
    info.mat <- vector("list", cats)
    for (i in 1:cats) {
      # select a function for each score category
      fun.tmp <- funList[[i]]

      # implement the function for each score category
      tmp <- do.call("fun.tmp", args.list, envir = environment())

      # compute the kernel of integration (the multinomial information of each category)
      ps.d1 <- attributes(tmp)$gradient
      info.list <- purrr::map(
        .x = 1:length(theta),
        .f = function(k) (1 / ps[k, i]) * outer(ps.d1[k, ], ps.d1[k, ]) * dens[k]
      )
      info.mat[[i]] <- Reduce(f = "+", x = info.list)
    }

    info.mat <- Reduce(f = "+", x = info.mat)
  }

  info.mat
}
