#' S-X2 Fit Statistic
#'
#' @description Computes the \eqn{S\text{-}X^2} item fit statistic proposed by
#'   Orlando and Thissen (2000, 2003). This statistic evaluates the fit of IRT
#'   models by comparing observed and expected item response frequencies across
#'   summed score groups.
#'
#' @inheritParams est_score
#' @param data A matrix of examinees' item responses, with rows for examinees
#'   and columns for items in the same order as the items in `x`. Responses to
#'   an item with K score categories must be integers from 0 to K - 1. Missing
#'   responses must be coded as `NA`; they are replaced with 0 before the
#'   statistics are computed, with a warning.
#' @param alpha A numeric value specifying the significance level (\eqn{\alpha})
#'   for the hypothesis test associated with the \eqn{S\text{-}X^2} statistic.
#'   Default is 0.05.
#' @param min.collapse A numeric value giving the minimum expected frequency
#'   per cell. Cells with smaller expected frequencies are merged with adjacent
#'   cells. Default is 1. See **Details**.
#' @param norm.prior A numeric vector of length two giving the mean and standard
#'   deviation of the normal latent ability distribution used to generate the
#'   quadrature points and weights. Ignored when `weights` is supplied. For an
#'   `est_irt` object, see `weights`. Default is `c(0, 1)`.
#' @param nquad An integer specifying the number of Gaussian quadrature points
#'   used to approximate the normal latent ability distribution. Ignored when
#'   `weights` is supplied. For an `est_irt` object, see `weights`. Default is
#'   30.
#' @param weights A two-column matrix or data frame containing the quadrature
#'   points (first column) and their corresponding weights (second column) for
#'   the latent ability distribution. If omitted, default values are generated
#'   using [irtQ::gen.weight()] according to the `norm.prior` and `nquad`
#'   arguments.
#'
#'   For an `est_irt` object, when `weights`, `norm.prior`, and `nquad` are all
#'   omitted, the latent distribution stored in the object (`x$weights`), which
#'   was used or estimated in the calibration, is used. Otherwise, `weights` is
#'   used if given, and a normal distribution given by `norm.prior` and `nquad`
#'   is used if not. For an `est_item` object, which stores no latent
#'   distribution, the normal distribution is used when `weights` is omitted.
#' @param pcm.loc An optional integer vector giving the row indices of partial
#'   credit model (PCM) items whose slope parameters are fixed. It is used only
#'   to count the item parameters for the degrees of freedom. Default is `NULL`.
#' @param ... Additional arguments passed to or from other methods.
#'
#' @details
#' For item \eqn{i} and summed score group \eqn{s} with \eqn{N_s} examinees,
#' let \eqn{O_{isk}} and \eqn{E_{isk}} be the observed and expected proportions
#' of score category \eqn{k}. The statistic is
#' \eqn{S\text{-}X^2_i = \sum_s \sum_k N_s (O_{isk} - E_{isk})^2 / E_{isk}},
#' which has the form of Orlando and Thissen (2000) for a dichotomous item.
#' The expected proportions are
#' \eqn{E_{isk} = \sum_q P_{ik}(\theta_q) f^{*i}(s - k \mid \theta_q) A(\theta_q)
#' / \sum_q f(s \mid \theta_q) A(\theta_q)},
#' where \eqn{f} and \eqn{f^{*i}} are the summed score likelihoods with and
#' without item \eqn{i}, computed by the Lord-Wingersky recursion, and
#' \eqn{\theta_q} and \eqn{A(\theta_q)} are the quadrature points and weights.
#'
#' The latent ability distribution used to compute the expected proportions is
#' given by `weights`, or by a normal distribution generated from `norm.prior`
#' and `nquad` when `weights` is omitted. For an `est_irt` object, the latent
#' distribution stored in the object is used when `weights`, `norm.prior`, and
#' `nquad` are all omitted.
#'
#' The lowest and highest possible summed scores are excluded. For an item with
#' \eqn{K} score categories, the \eqn{K - 1} lowest remaining summed scores are
#' pooled into one group, and so are the \eqn{K - 1} highest. When the test is
#' so short that the two pooled ranges overlap, each summed score is counted
#' once, in the lower group.
#'
#' The accuracy of the \eqn{\chi^{2}} approximation can be compromised when
#' expected cell frequencies are too small (Orlando & Thissen, 2000). For
#' dichotomous items, Orlando and Thissen (2000) merged adjacent summed score
#' groups so that every expected frequency is at least 1. For polytomous items,
#' this approach can discard too much information, so Kang and Chen (2008)
#' instead merged adjacent score categories *within* each summed
#' score group. [irtQ::sx2_fit()] follows both strategies, working from the ends
#' of the table toward the middle, and `min.collapse` sets the minimum expected
#' frequency.
#'
#' The degrees of freedom are \eqn{G(K - 1) - m - C}, where \eqn{G} is the
#' number of summed score groups after collapsing, \eqn{m} is the number of
#' item parameters, and \eqn{C} is the number of category cells removed by
#' merging. The number of item parameters \eqn{m} is determined by the IRT model
#' of each item: 1 for the 1PLM, 2 for the 2PLM, 3 for the 3PLM, and the number
#' of score categories \eqn{K} for the GRM and GPCM; an item listed in
#' `pcm.loc` has \eqn{m = K - 1}, and an item labeled `"DRM"` is counted as a
#' 3PLM item. This also applies to `est_irt` and `est_item` objects: parameters
#' fixed during estimation (for example, with `fix.g = TRUE`, or the items fixed
#' in FIPC) are still counted, and GPCM items estimated with
#' `fix.a.gpcm = TRUE` are counted as PCM items only when they are given in
#' `pcm.loc`. When no degrees of freedom remain after collapsing (`df` of 0 or
#' less), `crit.val` and `p` are `NA`, and a warning names the item.
#'
#' The test must contain at least two items. Responses in `data` must be
#' integer scores from 0 to the number of score categories minus 1, with one
#' column per item. Character and factor responses are read as numbers, logical
#' responses as 0 and 1, and other responses stop the function with an error.
#' Missing responses in `data` are replaced with 0 (the lowest score category),
#' and a warning is issued.
#'
#' @return
#' A list with the following components:
#' \item{fit_stat}{A data frame with one row per item and the columns `id`,
#' `chisq` (the \eqn{S\text{-}X^2} statistic), `df`, `crit.val` (the critical
#' value at `alpha`), and `p` (the p-value).}
#' \item{item_df}{The item metadata used in the analysis.}
#' \item{exp_freq}{A list of collapsed expected frequency tables, one per item,
#' with summed score groups in rows and score categories in columns.}
#' \item{obs_freq}{A list of collapsed observed frequency tables in the same
#' layout as `exp_freq`.}
#' \item{exp_prob}{A list of matrices of expected proportions within each
#' summed score group.}
#' \item{obs_prop}{A list of matrices of observed proportions within each
#' summed score group.}
#'
#' For a polytomous item, categories merged within a summed score group are
#' stored in the leftmost columns of that row, and the remaining cells are
#' `NA`.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [irtQ::irtfit()], [irtQ::simdat()], [irtQ::shape_df()],
#' [irtQ::est_irt()], [irtQ::est_item()]
#'
#' @references Kang, T., & Chen, T. T. (2008). Performance of the generalized
#' S-X2 item fit index for polytomous IRT models.
#' *Journal of Educational Measurement, 45*(4), 391-406.
#'
#' Orlando, M., & Thissen, D. (2000). Likelihood-based item-fit indices for
#' dichotomous item response theory models.
#' *Applied Psychological Measurement, 24*(1), 50-64.
#'
#' Orlando, M., & Thissen, D. (2003). Further investigation of the performance
#' of S-X2: An item fit index for use with dichotomous item response theory
#' models. *Applied Psychological Measurement, 27*(4), 289-298.
#'
#' @examples
#' ## Example 1: All five polytomous IRT items follow the GRM
#' ## Import the "-prm.txt" output file from flexMIRT
#' flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
#'
#' # Select the item metadata
#' x <- bring.flexmirt(file = flex_sam, "par")$Group1$full_df
#'
#' # Generate examinees' abilities from N(0, 1)
#' set.seed(23)
#' score <- rnorm(500, mean = 0, sd = 1)
#'
#' # Simulate response data
#' data <- simdat(x = x, theta = score, D = 1)
#'
#' \donttest{
#' # Compute fit statistics
#' fit1 <- sx2_fit(x = x, data = data, nquad = 30)
#'
#' # Display fit statistics
#' fit1$fit_stat
#' }
#'
#' ## Example 2: Items 39 and 40 follow the GRM, and items 53, 54, and 55
#' ##            follow the PCM (with slope parameters fixed to 1)
#' # Replace the model names with "GPCM" and
#' # set the slope parameters of items 53-55 to 1
#' x[53:55, 3] <- "GPCM"
#' x[53:55, 4] <- 1
#'
#' # Generate examinees' abilities from N(0, 1)
#' set.seed(25)
#' score <- rnorm(1000, mean = 0, sd = 1)
#'
#' # Simulate response data
#' data <- simdat(x = x, theta = score, D = 1)
#'
#' \donttest{
#' # Compute fit statistics
#' fit2 <- sx2_fit(x = x, data = data, nquad = 30, pcm.loc = 53:55)
#'
#' # Display fit statistics
#' fit2$fit_stat
#' }
#'
#' @export
sx2_fit <- function(x, ...) UseMethod("sx2_fit")


#' @describeIn sx2_fit Default method for computing \eqn{S\text{-}X^{2}} fit
#' statistics from a data frame `x` containing item metadata.
#' @importFrom Rfast rowsums
#' @import dplyr
#' @export
sx2_fit.default <- function(x,
                            data,
                            D = 1,
                            alpha = 0.05,
                            min.collapse = 1,
                            norm.prior = c(0, 1),
                            nquad = 30,
                            weights,
                            pcm.loc = NULL,
                            ...) {
  
  ## ------------------------------------------------------------------------------------------------
  # confirm and correct all item metadata information
  x <- confirm_df(x)

  # stop when the test has fewer than two items
  if (nrow(x) < 2L) {
    stop("'x' must contain at least two items to compute the S-X2 statistic.", call. = FALSE)
  }

  # transform the response data to a numeric matrix and check the responses
  data <- resp_to_matrix(data, x$cats, x$id)

  # check missing data
  # replace NAs with 0
  na.lg <- is.na(data)
  if (any(na.lg)) {
    # find the items that have no response at all
    allmiss <- which(colSums(!na.lg) == 0L)
    data[na.lg] <- 0

    # create the warning message, naming the items without any response
    memo <- "Missing responses are replaced with 0."
    if (length(allmiss) > 0L) {
      memo <- paste0(memo, " Every response is missing for item(s) ", paste(x$id[allmiss], collapse = ", "), ".")
    }
    warning(memo, call. = FALSE)
  }

  # break down the item metadata into several elements
  elm_item <- breakdown(x)

  # classify the items into DRM and PRM item groups
  idx.item <- idxfinder(elm_item)
  idx.drm <- idx.item$idx.drm
  idx.prm <- idx.item$idx.prm

  ## ------------------------------------------------------------------
  ## 1. data preparation
  ## ------------------------------------------------------------------
  # keep a copy of the item metadata
  full_df <- x
  
  # compute raw sum scores for all examinees
  rawscore <- Rfast::rowsums(data)
  
  # generate weights and thetas
  if (missing(weights)) {
    wts <- gen.weight(n = nquad, dist = "norm", mu = norm.prior[1], sigma = norm.prior[2])
  } else {
    wts <- data.frame(weights)
    nquad <- nrow(wts)
  }
  
  # select category variable
  cats <- elm_item$cats
  
  # total sum score
  t.score <- sum(cats - 1)
  
  # frequencies of raw sum scores
  score.freq <- as.numeric(table(factor(rawscore, levels = c(0:t.score))))
  
  ## ------------------------------------------------------------------
  ## 2. calculate the likelihoods
  ## ------------------------------------------------------------------
  # compute category probabilities for all items
  prob.cats <- trace(elm_item = elm_item, theta = wts[, 1], D = D, tcc = FALSE)$prob.cats
  
  # compute the summed score likelihoods with the Lord-Wingersky recursion
  lkhd <- t(lwRecurive(prob.cats = prob.cats, cats = cats, n.theta = nquad))
  
  # compute lkhd_noitem for all items via forward-backward LW (2 passes + J convolutions)
  lkhd_noitem <- lwrc_noitem(prob.cats = prob.cats, cats = cats, n.theta = nquad)
  
  ## ------------------------------------------------------------------
  ## 3. prepare the contingency tables
  ## ------------------------------------------------------------------
  # contingency tables of the expected frequencies for all items
  exp_freq <-
    purrr::pmap(
      .l = list(x = cats, y = prob.cats, z = lkhd_noitem),
      .f = function(x, y, z) {
        expFreq(t.score,
                cats = x, prob.cats = y,
                lkhd_noitem = z, lkhd, wts, score.freq
        )
      }
    )
  
  # contingency tables of the observed frequencies for all items
  obs_freq <-
    purrr::map2(
      .x = data.frame(data), .y = cats,
      .f = function(x, y) {
        obsFreq(
          rawscore = rawscore, response = x,
          t.score = t.score, cats = y
        )
      }
    )
  
  ## ------------------------------------------------------------------
  ## 4. collapse cells in the contingency tables
  ## ------------------------------------------------------------------
  # cbind two lists of the expected frequency tables and the observed frequency tables
  ftable_info <- cbind(exp_freq, obs_freq)
  
  # collapse the expected and observed frequency tables
  # (1) DRM items
  if (!is.null(idx.drm)) {
    # select frequency tables of dichotomous items
    if (length(idx.drm) > 1) {
      ftable_info_drm <- ftable_info[idx.drm, ]
    } else {
      ftable_info_drm <- rbind(ftable_info[idx.drm, ])
    }
    
    # collapse the frequency tables for all dichotomous items
    for (i in 1:length(idx.drm)) {
      ftab <- data.frame(ftable_info_drm[i, ])

      # collapse based on the column of the incorrect responses
      tmp1 <- collapse_ftable(x = ftab, col = 1, min.collapse = min.collapse)
      
      # collapse based on the column of the correct responses
      tmp2 <- collapse_ftable(x = tmp1, col = 2, min.collapse = min.collapse)
      
      # replace the frequency tables with the collapsed frequency tables
      ftable_info_drm[i, ][[1]] <- tmp2[, 1:2]
      ftable_info_drm[i, ][[2]] <- tmp2[, 3:4]
    }
    
    ftable_info[idx.drm, ] <- ftable_info_drm
  }
  
  # (2) PRM items
  if (!is.null(idx.prm)) {
    # select frequency tables of polytomous items
    if (length(idx.prm) > 1) {
      ftable_info_plm <- ftable_info[idx.prm, ]
    } else {
      ftable_info_plm <- rbind(ftable_info[idx.prm, ])
    }
    
    # collapse the frequency tables for all polytomous items
    for (i in 1:length(idx.prm)) {
      
      # select the expected and observed frequency tables for the corresponding items
      exp_tmp <- data.frame(ftable_info_plm[[i, 1]])
      obs_tmp <- data.frame(ftable_info_plm[[i, 2]])
      
      # drop summed score groups with no examinees
      if (any(rowSums(exp_tmp) == 0L)) {
        exp_tmp <- exp_tmp[rowSums(exp_tmp) != 0L, ]
        obs_tmp <- obs_tmp[rowSums(obs_tmp) != 0L, ]
      }
      
      # keep the column names before collapsing the frequency tables
      col.name <- colnames(exp_tmp)
      
      # collapse the categories of each summed score group with collapse_ftable_prm()
      out_tables <- collapse_ftable_prm(exp_mat = exp_tmp, obs_mat = obs_tmp,
                                        min.collapse = min.collapse)
      exp_table <- out_tables$exp_table
      obs_table <- out_tables$obs_table
      
      ftable_info_plm[i, ][[1]] <-
        data.frame(bind.fill(exp_table, type = "rbind")) %>%
        stats::setNames(nm = col.name[1:ncol(.)])
      ftable_info_plm[i, ][[2]] <- 
        data.frame(bind.fill(obs_table, type = "rbind")) %>%
        stats::setNames(nm = col.name[1:ncol(.)])
    }
    
    ftable_info[idx.prm, ] <- ftable_info_plm
  }
  
  # replace the two frequency tables with the collapsed two frequency tables
  exp_freq2 <-
    ftable_info[, 1] %>%
    purrr::map(.f = function(x) {
      dplyr::rename_all(x, .funs = list("gsub"), pattern = "exp_freq.", replacement = "")
    })
  
  obs_freq2 <-
    ftable_info[, 2] %>%
    purrr::map(.f = function(x) {
      dplyr::rename_all(x, .funs = list("gsub"), pattern = "obs_freq.", replacement = "")
    })
  
  ## ------------------------------------------------------------------
  ## 5. calculate the fit statistics
  ## ------------------------------------------------------------------
  # check the number of parameters for each item
  model <- full_df$model
  model[pcm.loc] <- "PCM"
  count_prm <- rep(NA, length(model))
  count_prm[model %in% "1PLM"] <- 1
  count_prm[model %in% "2PLM"] <- 2
  count_prm[model %in% c("3PLM", "DRM")] <- 3
  count_prm[model %in% "PCM"] <- full_df[model %in% "PCM", 2] - 1
  count_prm[model %in% "GPCM"] <- full_df[model %in% "GPCM", 2]
  count_prm[model %in% "GRM"] <- full_df[model %in% "GRM", 2]
  
  # compute the fit statistics for all items
  infoList <- list(exp_freq2, obs_freq2, as.list(count_prm))
  fitstat_list <- purrr::pmap(.l = infoList, .f = chisq_stat, alpha = alpha)
  
  # make a data.frame for the fit statistics results
  fit_stat <- list(chisq = NULL, df = NULL, crit.val = NULL)
  for (i in 1:3) {
    fit_stat[[i]] <- purrr::map_dbl(fitstat_list, .f = function(x) x[[i]])
  }
  pval <- purrr::pmap_dbl(
    .l = list(x = fit_stat[[1]], y = fit_stat[[2]]),
    .f = function(x, y) if (y > 0) 1 - stats::pchisq(q = x, df = y, lower.tail = TRUE) else NA_real_
  )

  # warn when an item has no degrees of freedom left for the test
  if (any(fit_stat$df <= 0)) {
    warning(
      "No degrees of freedom remain for the S-X2 statistic of item(s) ",
      paste(full_df$id[fit_stat$df <= 0], collapse = ", "),
      ". Their critical values and p-values are set to NA.",
      call. = FALSE
    )
  }

  fit_stat <- data.frame(id = full_df$id, fit_stat, p = round(pval, 3), stringsAsFactors = FALSE)
  fit_stat$chisq <- round(fit_stat$chisq, 3)
  fit_stat$crit.val <- round(fit_stat$crit.val, 3)
  rownames(fit_stat) <- NULL
  
  # extract the expected and observed proportion tables
  exp_prob <- purrr::map(fitstat_list, .f = function(x) x$exp.prop)
  obs_prop <- purrr::map(fitstat_list, .f = function(x) x$obs.prop)
  
  ## ------------------------------------------------------------------
  ## 6. return the results
  ## ------------------------------------------------------------------
  list(
    fit_stat = fit_stat, item_df = full_df, exp_freq = exp_freq2, obs_freq = obs_freq2,
    exp_prob = exp_prob, obs_prop = obs_prop
  )
}


#' @describeIn sx2_fit Method for an object of class `est_item` created by
#' [irtQ::est_item()]. The item parameter estimates, response data, and `D`
#' are taken from the object.
#' @import dplyr
#' @export
sx2_fit.est_item <- function(x, alpha = 0.05, min.collapse = 1, norm.prior = c(0, 1),
                            nquad = 30, weights, pcm.loc = NULL, ...) {
  # compute the fit statistics with the data and scaling constant stored in the object
  sx2_fit.default(
    x = x$par.est, data = x$data, D = x$scale.D, alpha = alpha,
    min.collapse = min.collapse, norm.prior = norm.prior, nquad = nquad,
    weights = weights, pcm.loc = pcm.loc
  )
}

#' @describeIn sx2_fit Method for an object of class `est_irt` created by
#' [irtQ::est_irt()]. The item parameter estimates, response data, and `D`
#' are taken from the object. When `weights`, `norm.prior`, and `nquad` are all
#' omitted, the latent distribution stored in the object (`x$weights`) is used.
#' @importFrom Rfast rowsums
#' @import dplyr
#' @export
sx2_fit.est_irt <- function(x, alpha = 0.05, min.collapse = 1, norm.prior = c(0, 1),
                            nquad = 30, weights, pcm.loc = NULL, ...) {
  # use the latent distribution stored in the object when no distribution is given
  if (missing(weights) && missing(norm.prior) && missing(nquad)) {
    weights <- x$weights
  }

  # compute the fit statistics with the data and scaling constant stored in the object
  sx2_fit.default(
    x = x$par.est, data = x$data, D = x$scale.D, alpha = alpha,
    min.collapse = min.collapse, norm.prior = norm.prior, nquad = nquad,
    weights = weights, pcm.loc = pcm.loc
  )
}
