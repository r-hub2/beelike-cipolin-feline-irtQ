#' Residual-Based DIF Detection Framework Using Categorical Residuals (RDIF-CR)
#'
#' This function computes three statistics of the residual-based DIF detection
#' framework using categorical residuals (RDIF-CR), \eqn{RDIF_{R}-CR},
#' \eqn{RDIF_{S}-CR}, and \eqn{RDIF_{RS}-CR}, for detecting global differential
#' item functioning (DIF), particularly in polytomously scored items. The RDIF-CR
#' framework evaluates DIF by comparing categorical residual vectors, which are
#' calculated as the difference between a one-hot encoded response vector
#' (with 1 for the selected category and 0 for all others) and the vector of
#' IRT model-predicted probabilities across all score categories.
#' This approach enables fine-grained detection of global DIF patterns at the
#' category level.
#'
#' @inheritParams rdif
#' @param score A numeric vector of examinees' ability estimates (theta values),
#'   in the same row order as `data`. The estimates should be based on the item
#'   parameters in `x` (see **Details**). If `NULL`, abilities are estimated
#'   internally with [irtQ::est_score()] using `method`, `range`, `norm.prior`,
#'   `nquad`, `weights`, and `ncore`. When `purify = TRUE`, abilities are
#'   re-estimated internally at every purification iteration, so a supplied
#'   `score` is used only in the initial analysis. A missing value in `score`
#'   excludes the examinee from the analysis, also during purification, with a
#'   warning. The `est_item` method has no `score` argument; it uses the
#'   abilities stored in the object. Default is `NULL`.
#' @param alpha A numeric value specifying the significance level (\eqn{\alpha})
#'   of the tests. An item is flagged by a statistic when its p-value is less
#'   than or equal to `alpha`. The comparison uses unrounded p-values, although
#'   the p-values reported in `dif_stat` are rounded to four decimal places.
#'   Default is `0.05`.
#' @param purify.by A character string specifying the RDIF-CR statistic used for
#'   purification: `"crdifrs"` for \eqn{RDIF_{RS}-CR}, `"crdifr"` for
#'   \eqn{RDIF_{R}-CR}, or `"crdifs"` for \eqn{RDIF_{S}-CR}. Used only when
#'   `purify = TRUE`. Default is `"crdifrs"`.
#' @param min.resp A single nonnegative number, usually a positive integer,
#'   specifying the minimum number of item responses that an examinee must have
#'   to be used in the analysis. All responses of examinees with fewer than
#'   `min.resp` (but at least one) responses are set to `NA` before the ability
#'   estimation, in the initial analysis and at every purification iteration,
#'   also when `score` is supplied, and these examinees are excluded from the
#'   analysis. A warning reports the number of examinees excluded in the initial
#'   analysis, and one warning at the end of the purification reports the
#'   examinees excluded during it. A large `min.resp` can leave a group without
#'   any ability estimate in a later purification iteration, which stops the
#'   purification with a warning. If `NULL`, no minimum is applied. Default is
#'   `NULL`. See **Details** for more information.
#'
#' @details
#' According to Penfield (2010), differential item functioning (DIF) in
#' polytomously scored items can be conceptualized in two forms: global DIF and
#' net DIF. Global DIF refers to differences between groups in the conditional
#' probabilities of responding in specific score categories, thus offering a
#' fine-grained view of DIF at the category level. In contrast, net DIF
#' summarizes these differences into a single value representing the overall
#' impact of DIF on the item's expected score.
#'
#' The RDIF framework using categorical residuals (RDIF-CR), implemented in
#' [irtQ::crdif()], extends the original residual-based DIF framework proposed
#' by Lim et al. (2022) to detect global DIF in polytomous items (Jung & Lim,
#' 2026; Lim, Malatesta, & Lee, 2024). This framework
#' includes three statistics: \eqn{RDIF_{R}-CR}, \eqn{RDIF_{S}-CR}, and
#' \eqn{RDIF_{RS}-CR}, each designed to capture different aspects of group-level
#' differences in categorical response patterns.
#'
#' To illustrate how the RDIF-CR framework operates, consider an item with five
#' ordered score categories (\eqn{k \in \{0,1,2,3,4\}}). Suppose an examinee with
#' latent ability \eqn{\theta} responds with category 2. The one-hot encoded
#' response vector for this response is \eqn{(0,0,1,0,0)^T}. Assume that the IRT
#' model predicts the category probabilities as
#' \eqn{(0.1, 0.2, 0.4, 0.25, 0.05)^T}, so that the expected item score is
#' \eqn{0(0.1) + 1(0.2) + 2(0.4) + 3(0.25) + 4(0.05) = 1.95}.
#' In the RDIF-CR framework, the categorical residual vector is calculated by
#' subtracting the predicted probability vector from the one-hot response vector,
#' resulting in \eqn{(-0.1, -0.2, 0.6, -0.25, -0.05)^T}.
#'
#' In contrast, net DIF is assessed with a single item score residual, the
#' observed item score minus the expected item score. In this example, the
#' residual would be \eqn{2 - 1.95 = 0.05}. To detect net DIF, use
#' [irtQ::rdif()].
#'
#' \eqn{RDIF_{R}-CR} is the quadratic form of the vector of category-level
#' differences in mean raw residuals between the focal and reference groups,
#' centered at its null mean, in a generalized inverse of its null covariance
#' matrix. \eqn{RDIF_{S}-CR} is the corresponding statistic for the mean squared
#' residuals, and \eqn{RDIF_{RS}-CR} uses the two vectors combined. Under the
#' null hypothesis of no DIF, each statistic asymptotically follows a chi-square
#' distribution whose degrees of freedom equal the rank of the covariance matrix
#' (Moore, 1977). For an item with K score categories (K of 3 or more), the K
#' category residuals of each examinee sum to zero, so the elements of the
#' \eqn{RDIF_{R}-CR} vector also sum to zero and its covariance matrix has rank
#' K - 1. \eqn{RDIF_{R}-CR} therefore has K - 1 degrees of freedom. It is
#' computed after dropping one score category, which gives the same value as the
#' Moore-Penrose generalized inverse of the full covariance matrix, and the same
#' category is dropped from the raw residual part of \eqn{RDIF_{RS}-CR}. The
#' squared residuals are not subject to a constraint that holds for all
#' examinees, so \eqn{RDIF_{S}-CR} has K degrees of freedom, and
#' \eqn{RDIF_{RS}-CR} has 2K - 1. For a dichotomous item, the two category
#' residuals are the negatives of each other, so only the residual for score
#' category 1 is used, and the three statistics have 1, 1, and 2 degrees of
#' freedom; then \eqn{RDIF_{R}-CR} and \eqn{RDIF_{S}-CR} equal the squares of the
#' standardized \eqn{RDIF_{R}} and \eqn{RDIF_{S}} statistics of [irtQ::rdif()],
#' \eqn{RDIF_{RS}-CR} equals \eqn{RDIF_{RS}}, and the p-values and flagged
#' items are identical to those of [irtQ::rdif()].
#'
#' In rare cases, for example when all examinees who responded to an item have
#' the same category probabilities, a covariance matrix has an even lower rank.
#' Then the statistic is computed with the Moore-Penrose generalized inverse and
#' compared with a chi-square distribution whose degrees of freedom equal the
#' rank, and a warning names the item. The columns `df.crdifr`, `df.crdifs`, and
#' `df.crdifrs` of `dif_stat` report the degrees of freedom used.
#'
#' The purification procedure, the `min.resp` argument, the handling of
#' examinees with `NA` ability estimates, the coding of `group`, and the `NaN`
#' results of an item to which no examinee of one group responded follow
#' [irtQ::rdif()] (see its Details). The p-value, not the statistic, is used to
#' choose the item to remove because the degrees of freedom of the RDIF-CR
#' statistics differ between items with different numbers of score categories.
#'
#'   The print method of this class has an argument `digits` that is not used;
#'   the statistics are printed rounded to three decimal places.
#'
#' @return This function returns an object of class `"crdif"`, which is a list
#' with the following five components:
#'
#' \item{no_purify}{A list of sub-objects containing the results of DIF analysis
#' without applying a purification procedure. The sub-objects include:
#'   \describe{
#'     \item{dif_stat}{A data frame with one row per item and the columns `id`,
#'     `crdifr`, `df.crdifr`, `crdifs`, `df.crdifs`, `crdifrs`, `df.crdifrs`,
#'     `p.crdifr`, `p.crdifs`, `p.crdifrs`, `n.ref`, `n.foc`, and `n.total`: the
#'     item ID, each RDIF-CR statistic followed by the degrees of freedom used, the
#'     p-values of the three statistics, the numbers of examinees in the
#'     reference and focal groups who responded to the item, and their sum.
#'     Statistics and p-values are rounded to four decimal places.}
#'     \item{moments}{A list with the elements `mu.crdifr`, `mu.crdifs`,
#'     `mu.crdifrs`, `cov.crdifr`, `cov.crdifs`, and `cov.crdifrs`. Each element
#'     is a list named by item ID that holds, for each item, the null mean
#'     vector or the null covariance matrix of the category-level statistics. For
#'     an item with K score categories (K of 3 or more), the vectors have lengths
#'     K, K, and 2K, and the matrices have dimensions K by K, K by K, and 2K by
#'     2K. For a dichotomous item, only score category 1 is used, so the lengths
#'     are 1, 1, and 2.}
#'     \item{dif_item}{A list with the elements `crdifr`, `crdifs`, and
#'     `crdifrs`, each giving the positions (rows of `x`) of the items flagged by
#'     the corresponding statistic, or `NULL` if no item is flagged.}
#'     \item{score}{A numeric vector of the ability estimates used in the initial
#'     analysis: the supplied `score`, or the internal estimates when
#'     `score = NULL`.}
#'   }
#' }
#'
#' \item{purify}{A logical value indicating whether purification was requested
#' (the `purify` argument).}
#'
#' \item{with_purify}{A list with the results of the purification procedure. All
#' elements are `NULL` when `purify = FALSE`. The elements are:
#'   \describe{
#'     \item{purify.by}{A character string indicating the RDIF-CR statistic used
#'     for purification. Possible values are "crdifr", "crdifs", or "crdifrs".}
#'     \item{dif_stat}{A data frame with the same columns as
#'     `no_purify$dif_stat` and an additional column `n.iter`. For an item
#'     removed during purification, the row reports the statistics from the
#'     iteration in which the item was removed, and `n.iter` is that iteration
#'     minus 1 (0 for an item removed on the basis of the initial analysis). For
#'     the other items, the row reports the statistics from the last iteration,
#'     and `n.iter` is the number of that iteration. If no item is flagged in the
#'     initial analysis, this is `no_purify$dif_stat` with `n.iter = 0`.}
#'     \item{moments}{A list with the same structure as `no_purify$moments`. For
#'     an item removed during purification, it holds the moments from the
#'     iteration in which the item was removed; for the other items, it holds the
#'     moments from the last iteration. If no item is flagged in the initial
#'     analysis, this is `no_purify$moments`.}
#'     \item{dif_item}{A numeric vector of the positions (rows of `x`) of the
#'     items flagged by the `purify.by` statistic, sorted in ascending order. It
#'     contains the items removed during purification and, if `max.iter` is
#'     reached or the purification stops because a group has no ability
#'     estimate left, the items flagged in the last iteration (iteration
#'     `n.iter`). It contains all items if every item is flagged. `NULL` if no
#'     item is flagged in the initial analysis.}
#'     \item{n.iter}{The number of purification iterations performed (0 if no
#'     item is flagged in the initial analysis). The iteration that removes the
#'     last item, or after which a group has no ability estimate left, is not
#'     counted, because the statistics are not recomputed.}
#'     \item{score}{A numeric vector of the ability estimates from the last
#'     iteration counted in `n.iter`. `NULL` if no item is flagged in the initial
#'     analysis, because purification is not carried out; the initial estimates
#'     are then in `no_purify$score`.}
#'     \item{complete}{A logical value. `TRUE` if the procedure stopped because
#'     no remaining item was flagged (including the case in which no item is
#'     flagged in the initial analysis), and `FALSE` if it stopped because
#'     `max.iter` was reached, because every item was flagged, or because a
#'     group had no examinee with an ability estimate left.}
#'   }
#' }
#'
#' \item{alpha}{A numeric value indicating the significance level (\eqn{\alpha})
#' used for the tests.}
#'
#' \item{call}{The matched function call.}
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [irtQ::rdif()], [irtQ::est_irt()], [irtQ::est_item()],
#'  [irtQ::est_score()]
#'
#' @references
#'   Bock, R. D., & Mislevy, R. J. (1982). Adaptive EAP estimation of ability in
#'   a microcomputer environment. *Applied Psychological Measurement, 6*(4),
#'   431-444. \doi{10.1177/014662168200600405}.
#'
#'   Hambleton, R. K., Swaminathan, H., & Rogers, H. J. (1991). *Fundamentals of
#'   item response theory*. Newbury Park, CA: Sage.
#'
#'   Jung, H., & Lim, H. (2026, April). Detecting global and net DIF in
#'   polytomous items using RDIF. Paper presented at the annual meeting of the
#'   National Council on Measurement in Education, Los Angeles, CA.
#'
#'   Lim, H., Choe, E. M., & Han, K. T. (2022). A residual-based differential
#'   item functioning detection framework in item response theory. *Journal of
#'   Educational Measurement, 59*(1), 80-104. \doi{10.1111/jedm.12313}.
#'
#'   Lim, H., Malatesta, J., & Lee, Y. (2024, July). Advancing polytomous DIF
#'   detection with the residual DIF framework. Paper presented at the annual
#'   International Meeting of the Psychometric Society, Prague, Czech Republic.
#'
#'   Moore, D. S. (1977). Generalized inverses, Wald's method, and the
#'   construction of chi-squared tests of fit. *Journal of the American
#'   Statistical Association, 72*(357), 131-137.
#'   \doi{10.1080/01621459.1977.10479921}.
#'
#'   Penfield, R. D. (2010). Distinguishing between net and global DIF in
#'   polytomous items. *Journal of Educational Measurement, 47*(2), 129-149.
#'   \doi{10.1111/j.1745-3984.2010.00105.x}.
#'
#'   Warm, T. A. (1989). Weighted likelihood estimation of ability in item
#'   response theory. *Psychometrika, 54*(3), 427-450.
#'   \doi{10.1007/BF02294627}.
#'
#' @examples
#' \donttest{
#'
#' ############################################################################
#' # This example demonstrates how to detect global DIF in polytomous items
#' # using the RDIF-CR framework implemented in `irtQ::crdif()`.
#' # Simulated response data are generated from 5 GRM items with 4 score
#' # categories. DIF is introduced in the 1st and 5th items.
#' ############################################################################
#'
#' ###############################################
#' # (1) Simulate response data with DIF
#' ###############################################
#'
#' set.seed(1)
#'
#' # Generate ability parameters for 1000 examinees in each group
#' # Reference and focal groups follow N(0, 1.5^2)
#' theta_ref <- rnorm(1000, 0, 1.5)
#' theta_foc <- rnorm(1000, 0, 1.5)
#'
#' # Combine abilities from both groups
#' theta_all <- c(theta_ref, theta_foc)
#'
#' # Define item parameters using `irtQ::shape_df()`
#' # Items 1 and 5 are intentionally modified to exhibit DIF
#' par_ref <- irtQ::shape_df(
#'   par.prm = list(
#'     a = c(1, 1, 1, 2, 2),
#'     d = list(c(-2, 0, 1),
#'              c(-2, 0, 2),
#'              c(-2, 0, 1),
#'              c(-1, 0, 2),
#'              c(-2, 0, 0.5))
#'   ),
#'   cats = 4, model = "GRM"
#' )
#'
#' par_foc <- irtQ::shape_df(
#'   par.prm = list(
#'     a = c(2, 1, 1, 2, 0.5),
#'     d = list(c(-0.5, 0, 0.5),
#'              c(-2, 0, 2),
#'              c(-2, 0, 1),
#'              c(-1, 0, 2),
#'              c(-1.5, -1, 0))
#'   ),
#'   cats = 4, model = "GRM"
#' )
#'
#' # Generate response data
#' resp_ref <- irtQ::simdat(x = par_ref, theta = theta_ref, D = 1)
#' resp_foc <- irtQ::simdat(x = par_foc, theta = theta_foc, D = 1)
#'
#' # Combine response data across groups
#' data <- rbind(resp_ref, resp_foc)
#'
#' ###############################################
#' # (2) Estimate item and ability parameters
#' ###############################################
#'
#' # Estimate GRM item parameters using `irtQ::est_irt()`
#' fit_mod <- irtQ::est_irt(data = data, D = 1, model = "GRM", cats = 4)
#'
#' # Extract estimated item parameters
#' x <- fit_mod$par.est
#'
#' # Estimate ability scores using ML method
#' score <- est_score(x = x, data = data, method = "ML")$est.theta
#'
#' ###############################################
#' # (3) Perform RDIF-CR DIF analysis
#' ###############################################
#'
#' # Define group membership: 1 = focal group
#' group <- c(rep(0, 1000), rep(1, 1000))
#'
#' # (a) DIF detection without purification
#' dif_nopuri <- crdif(
#'   x = x, data = data, score = score,
#'   group = group, focal.name = 1, D = 1, alpha = 0.05
#' )
#' print(dif_nopuri)
#'
#' # (b) DIF detection with purification using RDIF_{RS}-CR (not run: it takes
#' #     about 15 seconds because the abilities of the examinees are
#' #     re-estimated at every iteration)
#' #
#' #   dif_puri <- crdif(
#' #     x = x, data = data, score = score,
#' #     group = group, focal.name = 1, D = 1, alpha = 0.05,
#' #     purify = TRUE, purify.by = "crdifrs"
#' #   )
#' #   print(dif_puri)
#' #
#' #   purify.by = "crdifr" or "crdifs" can be used in the same way.
#'
#' }
#' @export
crdif <- function(x, ...) UseMethod("crdif")

#' @describeIn crdif Default method for a data frame `x` containing item metadata.
#'
#' @export
crdif.default <- function(x,
                          data,
                          score = NULL,
                          group,
                          focal.name,
                          item.skip = NULL,
                          D = 1,
                          alpha = 0.05,
                          missing = NA,
                          purify = FALSE,
                          purify.by = c("crdifrs", "crdifr", "crdifs"),
                          max.iter = 10,
                          min.resp = NULL,
                          method = "ML",
                          range = c(-5, 5),
                          norm.prior = c(0, 1),
                          nquad = 41,
                          weights = NULL,
                          ncore = 1,
                          verbose = TRUE,
                          ...) {

  # match.call
  cl <- match.call()

  # conduct the DIF analysis
  rst <- crdif_main(
    x = x, data = data, score = score, group = group, focal.name = focal.name,
    item.skip = item.skip, D = D, alpha = alpha, missing = missing, purify = purify,
    purify.by = purify.by, max.iter = max.iter, min.resp = min.resp, method = method,
    range = range, norm.prior = norm.prior, nquad = nquad, weights = weights,
    ncore = ncore, verbose = verbose, ...
  )

  # return the DIF detection results
  rst$call <- cl
  rst
}

#' @describeIn crdif Method for an object of class `est_irt` created by
#' [irtQ::est_irt()]. The response data and the scaling factor `D` are taken from
#' the object.
#'
#' @export
#'
crdif.est_irt <- function(x,
                          score = NULL,
                          group,
                          focal.name,
                          item.skip = NULL,
                          alpha = 0.05,
                          missing = NA,
                          purify = FALSE,
                          purify.by = c("crdifrs", "crdifr", "crdifs"),
                          max.iter = 10,
                          min.resp = NULL,
                          method = "ML",
                          range = c(-5, 5),
                          norm.prior = c(0, 1),
                          nquad = 41,
                          weights = NULL,
                          ncore = 1,
                          verbose = TRUE,
                          ...) {

  # match.call
  cl <- match.call()

  # stop when the response data or the scaling factor is passed, since they are taken from the object
  check_obj_dots(list(...))

  # conduct the DIF analysis with the data, scaling factor, and items of the object
  rst <- crdif_main(
    x = x$par.est, data = x$data, score = score, group = group, focal.name = focal.name,
    item.skip = item.skip, D = x$scale.D, alpha = alpha, missing = missing,
    purify = purify, purify.by = purify.by, max.iter = max.iter, min.resp = min.resp,
    method = method, range = range, norm.prior = norm.prior, nquad = nquad,
    weights = weights, ncore = ncore, verbose = verbose, ...
  )

  # return the DIF detection results
  rst$call <- cl
  rst
}

#' @describeIn crdif Method for an object of class `est_item` created by
#' [irtQ::est_item()]. The response data, ability estimates, and the scaling
#' factor `D` are taken from the object.
#'
#' @export
#'
crdif.est_item <- function(x,
                           group,
                           focal.name,
                           item.skip = NULL,
                           alpha = 0.05,
                           missing = NA,
                           purify = FALSE,
                           purify.by = c("crdifrs", "crdifr", "crdifs"),
                           max.iter = 10,
                           min.resp = NULL,
                           method = "ML",
                           range = c(-5, 5),
                           norm.prior = c(0, 1),
                           nquad = 41,
                           weights = NULL,
                           ncore = 1,
                           verbose = TRUE,
                           ...) {

  # match.call
  cl <- match.call()

  # stop when the response data, the ability estimates, or the scaling factor is passed, since they
  # are taken from the object
  check_obj_dots(list(...), obj_score = TRUE)

  # conduct the DIF analysis with the data, ability estimates, scaling factor, and items of the object
  rst <- crdif_main(
    x = x$par.est, data = x$data, score = x$score, group = group, focal.name = focal.name,
    item.skip = item.skip, D = x$scale.D, alpha = alpha, missing = missing,
    purify = purify, purify.by = purify.by, max.iter = max.iter, min.resp = min.resp,
    method = method, range = range, norm.prior = norm.prior, nquad = nquad,
    weights = weights, ncore = ncore, verbose = verbose, ...
  )

  # return the DIF detection results
  rst$call <- cl
  rst
}


# This function conducts the DIF analysis shared by the methods of crdif(): it prepares the
# data and ability estimates, runs the analysis once, and applies the purification procedure
crdif_main <- function(x, data, score, group, focal.name, item.skip, D, alpha, missing,
                       purify, purify.by, max.iter, min.resp, method, range, norm.prior,
                       nquad, weights, ncore, verbose, ...) {

  # confirm and correct all item metadata information
  x <- confirm_df(x)

  # record whether the ability estimates are supplied
  score_given <- !is.null(score)

  # check the inputs and prepare the response data and ability estimates
  prep <- dif_prepare(
    x = x, data = data, score = score, group = group, focal.name = focal.name,
    item.skip = item.skip, D = D, alpha = alpha, missing = missing, purify = purify,
    max.iter = max.iter, min.resp = min.resp, method = method, range = range,
    norm.prior = norm.prior, nquad = nquad, weights = weights, ncore = ncore, fn = "crdif()", ...
  )
  data <- prep$data
  score <- prep$score
  item.skip <- prep$item.skip

  # a) when no purification is set
  # do only one iteration of DIF analysis
  dif_rst <-
    crdif_one(x = x, data = data, score = score, group = group,
              focal.name = focal.name, item.skip = item.skip, D = D,
              alpha = alpha)

  # record the items whose covariance matrix is singular
  singular_id <- dif_rst$singular

  # create two empty lists to contain the results
  no_purify <- list(dif_stat = NULL, moments = NULL, dif_item = NULL, score = NULL)
  with_purify <- list(
    purify.by = NULL, dif_stat = NULL, moments = NULL,
    dif_item = NULL, n.iter = NULL, score = NULL, complete = NULL
  )

  # record the first DIF detection results into the no purification list
  no_purify$dif_stat <- dif_rst$dif_stat
  no_purify$dif_item <- dif_rst$dif_item
  no_purify$moments <- dif_rst$moments
  no_purify$score <- score

  # when purification is used
  if (purify) {
    # verify the criterion for purification
    purify.by <- match.arg(purify.by, c("crdifrs", "crdifr", "crdifs"))

    # create an empty vector and empty data frames
    # to contain the detected DIF items, statistics, and moments
    dif_item <- NULL
    dif_stat <-
      data.frame(
        id = rep(NA_character_, nrow(x)),
        crdifr = NA, df.crdifr = NA,
        crdifs = NA, df.crdifs = NA,
        crdifrs = NA, df.crdifrs = NA,
        p.crdifr = NA, p.crdifs = NA, p.crdifrs = NA,
        n.ref = NA, n.foc = NA, n.total = NA, n.iter = NA,
        stringsAsFactors = FALSE)
    mmt_list <-
      purrr::map(
        .x = dif_rst$moments,
        .f = ~ {
          purrr::map(.x = .x, .f = function(k) {
            NA
          })
        }
      )

    # extract the first DIF analysis results
    # and check if at least one DIF item is detected
    dif_item_tmp <- dif_rst$dif_item[[purify.by]]
    dif_stat_tmp <- dif_rst$dif_stat
    mmt_list_tmp <- no_purify$moments
    log_p_tmp <- dif_rst$log_p[[purify.by]]

    # copy the response data, item meta data, and ability estimates
    x_puri <- x
    data_puri <- data
    score_puri <- score

    # start the iteration if any item is detected as a DIF item
    if (!is.null(dif_item_tmp)) {
      # record unique item numbers
      item_num <- 1:nrow(x)

      # indicate whether the purification stops because every item is flagged
      all_flagged <- FALSE

      # indicate whether the purification stops because a group has no ability estimate left
      no_score <- FALSE

      # print a message
      if (verbose) {
        cat("Purification started...", "\n")
      }

      for (i in 1:max.iter) {
        # print a message
        if (verbose) {
          cat("\r", paste0("Iteration: ", i))
        }

        # find the flagged item with the smallest p-value
        flag_del <- pick_flagged_item(flag_loc = dif_item_tmp, log_p = log_p_tmp)

        # check an item that is deleted
        del_item <- item_num[flag_del]

        # add the deleted item as the DIF item
        dif_item <- c(dif_item, del_item)

        # add the DIF statistics and moments for the detected DIF item
        dif_stat[del_item, 1:13] <- dif_stat_tmp[flag_del, ]
        dif_stat[del_item, 14] <- i - 1
        mmt_list$mu.crdifr[[del_item]] <- mmt_list_tmp$mu.crdifr[[flag_del]]
        mmt_list$mu.crdifs[[del_item]] <- mmt_list_tmp$mu.crdifs[[flag_del]]
        mmt_list$mu.crdifrs[[del_item]] <- mmt_list_tmp$mu.crdifrs[[flag_del]]
        mmt_list$cov.crdifr[[del_item]] <- mmt_list_tmp$cov.crdifr[[flag_del]]
        mmt_list$cov.crdifs[[del_item]] <- mmt_list_tmp$cov.crdifs[[flag_del]]
        mmt_list$cov.crdifrs[[del_item]] <- mmt_list_tmp$cov.crdifrs[[flag_del]]

        # find the examinees who responded to the item to be deleted
        loc_resp <- which(!is.na(data_puri[, flag_del]))

        # refine the leftover items
        item_num <- item_num[-flag_del]

        # stop the purification when no item is left to estimate the abilities
        if (length(item_num) == 0L) {
          all_flagged <- TRUE
          break
        }

        # remove the detected DIF item data which has the smallest p-value from the item metadata
        x_puri <- x_puri[-flag_del, , drop = FALSE]

        # remove the detected DIF item data which has the smallest p-value from the response data
        data_puri <- data_puri[, -flag_del, drop = FALSE]

        # update the locations of the items that should be skipped in the purified data
        if (!is.null(item.skip)) {
          item.skip.puri <- c(1:length(item_num))[item_num %in% item.skip]
        } else {
          item.skip.puri <- NULL
        }

        # if min.resp is not NULL, find the examinees who have the number of responses
        # less than specified value (e.g., 5). Then, replace their all responses with NA
        if (!is.null(min.resp)) {
          n_resp <- rowSums(!is.na(data_puri))
          loc_less <- which(n_resp < min.resp & n_resp > 0)
          data_puri[loc_less, ] <- NA
        }

        # keep the ability estimates of the previous iteration in case the purification stops
        score_prev <- score_puri

        # compute the updated ability estimates after deleting the detected DIF item data;
        # only the examinees who responded to the deleted item are rescored after the first iteration
        score_puri <-
          dif_rescore(
            x = x_puri, data = data_puri, score = score_puri, loc_resp = loc_resp,
            first = (i == 1L), D = D, method = method, range = range,
            norm.prior = norm.prior, nquad = nquad, weights = weights, ncore = ncore,
            quiet_na = TRUE, quiet_ncore = (i > 1L || !score_given), ...
          )

        # stop the purification when the examinees of a group have all lost their ability estimates;
        # the remaining items keep the results of the previous iteration
        if (all(is.na(score_puri[which(group == focal.name)])) ||
            all(is.na(score_puri[which(group != focal.name)]))) {
          no_score <- TRUE
          score_puri <- score_prev
          dif_stat[item_num, 1:13] <- dif_stat_tmp[-flag_del, ]
          dif_stat[item_num, 14] <- i - 1
          mmt_list$mu.crdifr[item_num] <- mmt_list_tmp$mu.crdifr[-flag_del]
          mmt_list$mu.crdifs[item_num] <- mmt_list_tmp$mu.crdifs[-flag_del]
          mmt_list$mu.crdifrs[item_num] <- mmt_list_tmp$mu.crdifrs[-flag_del]
          mmt_list$cov.crdifr[item_num] <- mmt_list_tmp$cov.crdifr[-flag_del]
          mmt_list$cov.crdifs[item_num] <- mmt_list_tmp$cov.crdifs[-flag_del]
          mmt_list$cov.crdifrs[item_num] <- mmt_list_tmp$cov.crdifrs[-flag_del]

          # add the items flagged at the previous iteration, except the item deleted at this one
          dif_item <- c(dif_item, item_num[match(setdiff(dif_item_tmp, flag_del),
                                                 seq_along(log_p_tmp)[-flag_del])])
          break
        }

        # do DIF analysis using the updated ability estimates
        dif_rst_tmp <-
          crdif_one(
            x = x_puri, data = data_puri, score = score_puri, group = group,
            focal.name = focal.name, item.skip = item.skip.puri, D = D,
            alpha = alpha)

        # record the items whose covariance matrix is singular
        singular_id <- union(singular_id, dif_rst_tmp$singular)

        # extract the DIF analysis results
        # and check if at least one DIF item is detected
        dif_item_tmp <- dif_rst_tmp$dif_item[[purify.by]]
        dif_stat_tmp <- dif_rst_tmp$dif_stat
        mmt_list_tmp <- dif_rst_tmp$moments
        log_p_tmp <- dif_rst_tmp$log_p[[purify.by]]

        # check if a further DIF item is flagged
        if (is.null(dif_item_tmp)) {
          # add the DIF statistics for rest of items
          dif_stat[item_num, 1:13] <- dif_stat_tmp
          dif_stat[item_num, 14] <- i
          mmt_list$mu.crdifr[item_num] <- mmt_list_tmp$mu.crdifr
          mmt_list$mu.crdifs[item_num] <- mmt_list_tmp$mu.crdifs
          mmt_list$mu.crdifrs[item_num] <- mmt_list_tmp$mu.crdifrs
          mmt_list$cov.crdifr[item_num] <- mmt_list_tmp$cov.crdifr
          mmt_list$cov.crdifs[item_num] <- mmt_list_tmp$cov.crdifs
          mmt_list$cov.crdifrs[item_num] <- mmt_list_tmp$cov.crdifrs

          break
        }
      }

      # print a message
      if (verbose) {
        cat("", "\n")
      }

      # record the actual number of iterations; the statistics are not recomputed in the
      # iteration that removes the last item
      n_iter <- if (all_flagged || no_score) i - 1L else i

      # if every item is flagged, then, return a warning message
      if (all_flagged) {
        warning("All items were flagged during the purification, so the purification stopped ",
                "with no item left to analyze.", call. = FALSE)
        complete <- FALSE
      } else if (no_score) {
        # if a group has no ability estimate left, then, return a warning message
        warning("The purification stopped because no examinee of a group had enough responses ",
                "to the remaining items to estimate the ability.", call. = FALSE)
        complete <- FALSE
      } else if (max.iter == n_iter & !is.null(dif_item_tmp)) {
        # if the maximum number of iterations is reached before purification is complete,
        # then, return a warning message
        warning("The maximum number of iterations was reached before purification was completed.",
                call. = FALSE)
        complete <- FALSE

        # add flagged DIF item at the last iteration
        dif_item <- c(dif_item, item_num[dif_item_tmp])

        # add the DIF statistics for rest of items
        dif_stat[item_num, 1:13] <- dif_stat_tmp
        dif_stat[item_num, 14] <- i
        mmt_list$mu.crdifr[item_num] <- mmt_list_tmp$mu.crdifr
        mmt_list$mu.crdifs[item_num] <- mmt_list_tmp$mu.crdifs
        mmt_list$mu.crdifrs[item_num] <- mmt_list_tmp$mu.crdifrs
        mmt_list$cov.crdifr[item_num] <- mmt_list_tmp$cov.crdifr
        mmt_list$cov.crdifs[item_num] <- mmt_list_tmp$cov.crdifs
        mmt_list$cov.crdifrs[item_num] <- mmt_list_tmp$cov.crdifrs
      } else {
        complete <- TRUE

        # print a message
        if (verbose) {
          cat("Purification is finished.", "\n")
        }
      }

      # warn once about the examinees who lost their ability estimates during the purification
      warn_purify_excluded(sum(is.na(score_puri) & rowSums(!is.na(data)) > 0))

      # record the final DIF detection results with the purification procedure
      with_purify$purify.by <- purify.by
      with_purify$dif_stat <- dif_stat
      with_purify$moments <- mmt_list
      with_purify$dif_item <- sort(dif_item)
      with_purify$n.iter <- n_iter
      with_purify$score <- score_puri
      with_purify$complete <- complete
    } else {
      # in case when no DIF item is detected from the first DIF analysis results
      with_purify$purify.by <- purify.by
      with_purify$dif_stat <- cbind(no_purify$dif_stat, n.iter = 0)
      with_purify$moments <- no_purify$moments
      with_purify$n.iter <- 0L
      with_purify$complete <- TRUE
    }
  }

  # warn that the generalized inverse is used for the items with a singular covariance matrix
  if (length(singular_id) > 0L) {
    warning("A covariance matrix of the RDIF-CR statistics has a lower rank than expected for ",
            "item(s) ", paste(singular_id, collapse = ", "), ". The statistics of these items ",
            "are computed with a generalized inverse, and the rank of each matrix is used as ",
            "the degrees of freedom.", call. = FALSE)
  }

  # summarize the results
  rst <- list(no_purify = no_purify, purify = purify,
              with_purify = with_purify, alpha = alpha)

  # return the DIF detection results
  class(rst) <- "crdif"
  rst
}


# This function conducts one iteration of DIF analysis using the categorical residual-based statistics
crdif_one <- function(x,
                      data,
                      score,
                      group,
                      focal.name,
                      item.skip = NULL,
                      D = 1,
                      alpha = 0.05) {

  # break down the item metadata into several elements
  elm_item <- breakdown(x)

  # check the unique score categories
  cats <- cats.mod <- elm_item$cats

  # modify the score category for the binary items
  cats.mod[cats.mod == 2] <- 1

  # check the number of items
  nitem <- length(cats)

  # treat the responses of examinees without an ability estimate as missing
  na_score <- is.na(score)
  if (any(na_score)) {
    data[na_score, ] <- NA
    score[na_score] <- 0
  }

  ## ---------------------------------
  # compute the CRDIF statistics
  ## ---------------------------------
  # find the location of examinees for the reference and the focal groups
  loc_ref <- which(group != focal.name)
  loc_foc <- which(group == focal.name)

  # divide the response data into the two group data
  resp_ref <- data[loc_ref, , drop = FALSE]
  resp_foc <- data[loc_foc, , drop = FALSE]

  # check sample size
  n_ref <- Rfast::colsums(!is.na(resp_ref))
  n_foc <- Rfast::colsums(!is.na(resp_foc))

  # check if an item has all missing data for either of two groups
  all_miss <- sort(unique(c(which(n_ref == 0), which(n_foc == 0))))

  # divide the thetas into the two group data
  score_ref <- score[loc_ref]
  score_foc <- score[loc_foc]

  # compute the model probability of score categories
  # and keep the probabilities of the correct answers only for the binary items
  prob_ref <-
    trace(elm_item = elm_item, theta = score_ref, D = D, tcc = FALSE)$prob.cats %>%
    purrr::map2(,
                .y = cats.mod,
                .f = ~ {
                  if (.y == 1) {
                    .x[, 2, drop = FALSE]
                  } else {
                    .x
                  }
                }
    )
  prob_foc <-
    trace(elm_item = elm_item, theta = score_foc, D = D, tcc = FALSE)$prob.cats %>%
    purrr::map2(,
                .y = cats.mod,
                .f = ~ {
                  if (.y == 1) {
                    .x[, 2, drop = FALSE]
                  } else {
                    .x
                  }
                }
    )

  # create a list of the factorized response values across all items
  resp_ref <-
    purrr::map2(
      .x = data.frame(resp_ref, stringsAsFactors = FALSE), .y = cats,
      .f = function(k, m) factor(k, levels = (seq_len(m) - 1))
    )
  resp_foc <-
    purrr::map2(
      .x = data.frame(resp_foc, stringsAsFactors = FALSE), .y = cats,
      .f = function(k, m) factor(k, levels = (seq_len(m) - 1))
    )

  # omit the rows including NA values
  prob_ref <-
    purrr::map2(
      .x = prob_ref, .y = resp_ref,
      .f = ~ {
        .x[!is.na(.y), , drop = FALSE]
      }
    )
  prob_foc <-
    purrr::map2(
      .x = prob_foc, .y = resp_foc,
      .f = ~ {
        .x[!is.na(.y), , drop = FALSE]
      }
    )

  # create a contingency table of score categories for each item
  # and then, transform the table to a matrix format
  # also, omit the rows including NAs
  ref_id <- 1:length(score_ref)
  foc_id <- 1:length(score_foc)
  freqtab_ref <-
    purrr::map2(
      .x = resp_ref, .y = cats,
      .f = function(k, y) {
        tmp.tb <-
          stats::xtabs(~ ref_id + k,
                       na.action = stats::na.omit, addNA = FALSE
          ) %>%
          matrix(ncol = y)
        if (y == 2) {
          tmp.tb[, 2, drop = FALSE]
        } else {
          tmp.tb
        }
      }
    )
  freqtab_foc <-
    purrr::map2(
      .x = resp_foc, .y = cats,
      .f = function(k, y) {
        tmp.tb <-
          stats::xtabs(~ foc_id + k,
                       na.action = stats::na.omit, addNA = FALSE
          ) %>%
          matrix(ncol = y)
        if (y == 2) {
          tmp.tb[, 2, drop = FALSE]
        } else {
          tmp.tb
        }
      }
    )

  # compute the category residuals (one-hot responses minus category probabilities)
  resid_ref <- purrr::map2(.x = freqtab_ref, .y = prob_ref, .f = ~ {
    .x - .y
  })
  resid_foc <- purrr::map2(.x = freqtab_foc, .y = prob_foc, .f = ~ {
    .x - .y
  })

  # compute the residual-based DIF statistic
  rdifr <-
    purrr::map2(
      .x = resid_foc, .y = resid_ref,
      .f = ~ {
        colMeans(.x, na.rm = TRUE) - colMeans(.y, na.rm = TRUE)
      }
    )
  rdifs <-
    purrr::map2(
      .x = resid_foc, .y = resid_ref,
      .f = ~ {
        colMeans(.x^2, na.rm = TRUE) - colMeans(.y^2, na.rm = TRUE)
      }
    )

  # compute the means and variances of the two statistics for the hypothesis testing
  moments <- resid_moments_cr(
    p_ref = prob_ref, p_foc = prob_foc, n_ref = n_ref, n_foc = n_foc,
    resp_ref = resp_ref, resp_foc = resp_foc, cats.mod = cats.mod
  )
  mu_rdifr <- moments$mu.rdifr
  mu_rdifs <- moments$mu.rdifs
  cov_rdifr <- moments$cov.rdifr
  cov_rdifs <- moments$cov.rdifs
  cov_rdifrs <- moments$cov.rdifrs

  # degrees of freedom of the chi-square statistics: the ranks of the covariance matrices,
  # which are K - 1, K, and 2K - 1 for a polytomous item with K score categories
  # and 1, 1, and 2 for a dichotomous item
  df_r <- ifelse(cats.mod > 1, cats.mod - 1, 1)
  df_s <- cats.mod
  df_rs <- df_r + df_s

  # compute the chi-square statistics and their degrees of freedom
  chisq_r <- chisq_s <- chisq_rs <- rep(NaN, nitem)
  reduced_loc <- rep(FALSE, nitem)
  for (i in 1:nitem) {
    if (!i %in% all_miss) {
      # the raw residuals of the score categories sum to zero, so the first category is
      # dropped from the raw residual statistic of a polytomous item; this gives the same
      # quadratic form as the generalized inverse of the full covariance matrix
      loc_r <- if (cats.mod[i] > 1) 2:cats.mod[i] else 1
      loc_rs <- c(loc_r, cats.mod[i] + seq_len(cats.mod[i]))

      # create the vectors of the differences between the statistics and their means
      dev_r_all <- as.numeric(rdifr[[i]] - mu_rdifr[[i]])
      dev_s <- as.numeric(rdifs[[i]] - mu_rdifs[[i]])
      dev_r <- dev_r_all[loc_r]
      dev_rs <- c(dev_r_all, dev_s)[loc_rs]

      # compute the chi-square statistics with the covariance matrices
      qf_r <- quad_form(cov_mat = cov_rdifr[[i]][loc_r, loc_r, drop = FALSE], dev_vec = dev_r,
                        df = length(loc_r))
      qf_s <- quad_form(cov_mat = cov_rdifs[[i]], dev_vec = dev_s, df = length(dev_s))
      qf_rs <- quad_form(cov_mat = cov_rdifrs[[i]][loc_rs, loc_rs, drop = FALSE], dev_vec = dev_rs,
                         df = length(loc_rs))
      chisq_r[i] <- qf_r$stat
      chisq_s[i] <- qf_s$stat
      chisq_rs[i] <- qf_rs$stat

      # use the rank of a singular covariance matrix as the degrees of freedom
      df_r[i] <- qf_r$df
      df_s[i] <- qf_s$df
      df_rs[i] <- qf_rs$df
      reduced_loc[i] <- qf_r$reduced | qf_s$reduced | qf_rs$reduced
    }
  }

  # find the items whose covariance matrix is singular
  singular_id <- x$id[which(reduced_loc & !(seq_len(nitem) %in% item.skip))]

  # calculate p-values for all three statistics
  p_crdifr <- stats::pchisq(chisq_r, df = df_r, lower.tail = FALSE)
  p_crdifs <- stats::pchisq(chisq_s, df = df_s, lower.tail = FALSE)
  p_crdifrs <- stats::pchisq(chisq_rs, df = df_rs, lower.tail = FALSE)

  # compute the log p-values to choose the item to be removed in the purification
  log_p <- list(
    crdifr = stats::pchisq(chisq_r, df = df_r, lower.tail = FALSE, log.p = TRUE),
    crdifs = stats::pchisq(chisq_s, df = df_s, lower.tail = FALSE, log.p = TRUE),
    crdifrs = stats::pchisq(chisq_rs, df = df_rs, lower.tail = FALSE, log.p = TRUE)
  )

  # compute total sample size
  n_total <- n_foc + n_ref

  # create a data frame to contain the results
  stat_df <-
    data.frame(
      id = x$id,
      crdifr = round(chisq_r, 4),
      df.crdifr = df_r,
      crdifs = round(chisq_s, 4),
      df.crdifs = df_s,
      crdifrs = round(chisq_rs, 4),
      df.crdifrs = df_rs,
      p.crdifr = round(p_crdifr, 4),
      p.crdifs = round(p_crdifs, 4),
      p.crdifrs = round(p_crdifrs, 4),
      n.ref = n_ref,
      n.foc = n_foc,
      n.total = n_total,
      stringsAsFactors = FALSE
    )
  rownames(stat_df) <- NULL

  # data frame of the moments
  names(mu_rdifr) <-
    names(mu_rdifs) <-
    names(cov_rdifr) <-
    names(cov_rdifs) <-
    names(cov_rdifrs) <- x$id
  mmt_list <- list(
    mu.crdifr = mu_rdifr,
    mu.crdifs = mu_rdifs,
    mu.crdifrs = purrr::map2(
      .x = mu_rdifr, .y = mu_rdifs,
      .f = ~ {
        c(.x, .y)
      }),
    cov.crdifr = cov_rdifr,
    cov.crdifs = cov_rdifs,
    cov.crdifrs = cov_rdifrs
  )

  # when there are items that should be skipped for the DIF analysis
  # insert NAs to the corresponding results of the items
  if (!is.null(item.skip)) {
    stat_df[item.skip, 2:10] <- NA
    p_crdifr[item.skip] <- NA
    p_crdifs[item.skip] <- NA
    p_crdifrs[item.skip] <- NA
    log_p <- log_p %>%
      purrr::map(.f = ~{
        .x[item.skip] <- NA
        .x})
    mmt_list$mu.crdifr[item.skip] <- NA
    mmt_list$mu.crdifs[item.skip] <- NA
    mmt_list$mu.crdifrs[item.skip] <- NA
    mmt_list$cov.crdifr[item.skip] <- NA
    mmt_list$cov.crdifs[item.skip] <- NA
    mmt_list$cov.crdifrs[item.skip] <- NA
  }

  # find the flagged items
  dif_item_crdifr <- which(p_crdifr <= alpha)
  dif_item_crdifs <- which(p_crdifs <= alpha)
  dif_item_crdifrs <- which(p_crdifrs <= alpha)
  if (length(dif_item_crdifr) == 0) dif_item_crdifr <- NULL
  if (length(dif_item_crdifs) == 0) dif_item_crdifs <- NULL
  if (length(dif_item_crdifrs) == 0) dif_item_crdifrs <- NULL

  # summarize the results
  rst <- list(
    dif_stat = stat_df,
    dif_item = list(
      crdifr = dif_item_crdifr,
      crdifs = dif_item_crdifs,
      crdifrs = dif_item_crdifrs
    ),
    moments = mmt_list,
    alpha = alpha,
    log_p = log_p,
    singular = singular_id
  )

  # return the results
  rst

}

# This function computes the mean and covariance of the RDIF statistics across all category response data
resid_moments_cr <- function(p_ref, p_foc, n_ref, n_foc, resp_ref, resp_foc, cats.mod) {
  # check the number of items
  nitem <- length(cats.mod)

  # check the number of rows and columns in the two-group response data
  nrow_ref <- length(resp_ref[[1]])
  nrow_foc <- length(resp_foc[[1]])

  # compute the mean, variance, and covariance of rdif statistics
  mu_rdifr <-
    mu_rdifs <-
    sigma2_rdifr <-
    sigma2_rdifs <-
    cov_rdifr <-
    cov_rdifs <-
    cov_rs_diag <-
    cov_rs <-
    cov_all <- vector("list", nitem)
  for (i in 1:nitem) {
    # extract the requirements to compute the moments
    p.f <- p_foc[[i]]
    p.r <- p_ref[[i]]
    q.f <- 1 - p.f
    q.r <- 1 - p.r
    p.f2 <- p.f^2
    p.r2 <- p.r^2
    n.f <- n_foc[i]
    n.r <- n_ref[i]
    n.f2 <- n.f^2
    n.r2 <- n.r^2

    # compute the moments of statistics
    # (a) means of the rdifr and rdifs
    mu_rdifr[[i]] <- rep(0, cats.mod[i])
    mu_rdifs[[i]] <-
      purrr::map_dbl(
        .x = 1:cats.mod[i],
        .f = function(j) {
          (sum(p.f[, j] * q.f[, j], na.rm = TRUE) / n.f) -
            (sum(p.r[, j] * q.r[, j], na.rm = TRUE) / n.r)
        }
      )

    # (b) variances of the rdifr and rdifs
    sigma2_rdifr[[i]] <-
      purrr::map_dbl(
        .x = 1:cats.mod[i],
        .f = function(j) {
          (sum(p.f[, j] * q.f[, j], na.rm = TRUE) / n.f2) +
            (sum(p.r[, j] * q.r[, j], na.rm = TRUE) / n.r2)
        }
      )
    sigma2_rdifs[[i]] <-
      purrr::map_dbl(
        .x = 1:cats.mod[i],
        .f = function(j) {
          (sum(p.f[, j] * q.f[, j] * (1 - 2 * p.f[, j])^2, na.rm = TRUE) / n.f2) +
            (sum(p.r[, j] * q.r[, j] * (1 - 2 * p.r[, j])^2, na.rm = TRUE) / n.r2)
        }
      )

    # (c) covariances between the rdifr and rdifs for the same score categories
    cov_rs_diag[[i]] <-
      purrr::map_dbl(
        .x = 1:cats.mod[i],
        .f = function(j) {
          (sum(p.f[, j] * q.f[, j] * (1 - 2 * p.f[, j]), na.rm = TRUE) / n.f2) +
            (sum(p.r[, j] * q.r[, j] * (1 - 2 * p.r[, j]), na.rm = TRUE) / n.r2)
        }
      )

    # (d) covariances of rdifr and rdifs for the different score categories
    if (cats.mod[i] == 1) {
      cov_rdifr[[i]] <- matrix(sigma2_rdifr[[i]])
      cov_rdifs[[i]] <- matrix(sigma2_rdifs[[i]])
      cov_rs[[i]] <- matrix(cov_rs_diag[[i]])
    } else {
      cprod.p.f <- crossprod(p.f)
      cprod.p.r <- crossprod(p.r)
      cprod.p.f2 <- crossprod(p.f2)
      cprod.p.r2 <- crossprod(p.r2)
      cprod.pp2.f <- crossprod(x = p.f, y = p.f2)
      cprod.pp2.r <- crossprod(x = p.r, y = p.r2)
      # note that
      # t(p.f2) %*% p.f == crossprod(x = p.f2, y = p.f)
      # t(p.f) %*% p.f2 == crossprod(x = p.f, y = p.f2)
      # t(cprod.pp2.f) == t(p.f2) %*% p.f

      cov.r.tmp <- -{
        cprod.p.f / n.f2 + cprod.p.r / n.r2
      }
      diag(cov.r.tmp) <- sigma2_rdifr[[i]]
      cov_rdifr[[i]] <- cov.r.tmp

      cov.s.tmp <-
        {
          (-cprod.p.f - 4 * cprod.p.f2 + 2 * (t(cprod.pp2.f) + cprod.pp2.f)) / n.f2
        } + {
          (-cprod.p.r - 4 * cprod.p.r2 + 2 * (t(cprod.pp2.r) + cprod.pp2.r)) / n.r2
        }
      diag(cov.s.tmp) <- sigma2_rdifs[[i]]
      cov_rdifs[[i]] <- cov.s.tmp

      cov.rs.tmp <-
        {
          (-cprod.p.f + 2 * cprod.pp2.f) / n.f2
        } + {
          (-cprod.p.r + 2 * cprod.pp2.r) / n.r2
        }
      diag(cov.rs.tmp) <- cov_rs_diag[[i]]
      cov_rs[[i]] <- cov.rs.tmp
    }

    # (e) covariances of both rdifr and rdifs for the different score categories
    cov_all[[i]] <-
      rbind(
        cbind(cov_rdifr[[i]], cov_rs[[i]]),
        cbind(t(cov_rs[[i]]), cov_rdifs[[i]])
      )
  }

  # return the results
  rst <-
    list(
      mu.rdifr = mu_rdifr, mu.rdifs = mu_rdifs,
      cov.rdifr = cov_rdifr, cov.rdifs = cov_rdifs,
      cov.rdifrs = cov_all
    )
  rst
}
