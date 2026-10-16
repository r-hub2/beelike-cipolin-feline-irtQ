#' Pseudo-Count D2 Method for Item Parameter Drift Detection
#'
#' This function calculates the pseudo-count \eqn{D^{2}} statistic
#' to evaluate item parameter drift, as described by Cappaert et al. (2018) and
#' Stone (2000). The pseudo-count \eqn{D^{2}} statistic is designed to detect
#' item parameter drift efficiently without requiring item recalibration, making
#' it especially valuable in computerized adaptive testing (CAT) environments.
#' This method compares observed and expected response frequencies across
#' quadrature points, which represent latent ability levels. The observed
#' pseudo-frequencies are computed using the posterior distribution of each
#' examinee's ability, whereas the expected frequencies are the response
#' probabilities implied by the item parameters in the item bank (Stone, 2000).
#' The current version supports dichotomous items only.
#'
#' @inheritParams rdif
#' @inheritParams est_irt
#' @param x A data frame containing item metadata (e.g., item parameters, number
#'   of categories, IRT model types, etc.). See [irtQ::est_irt()] or
#'   [irtQ::simdat()] for more details about the item metadata. This data frame
#'   can be easily created using the [irtQ::shape_df()] function. Only
#'   dichotomous items are supported; the function stops if `x` contains GRM or
#'   GPCM items.
#' @param item.skip A numeric vector of item positions (row numbers of `x`), or a
#'   logical vector with one value per item, identifying the items to exclude
#'   from the IPD analysis. Skipped items still contribute to the
#'   posterior ability distributions; only their `pcd2` values are set to `NA`,
#'   and they are never flagged. If `NULL`, all items are analyzed. Default is
#'   `NULL`.
#' @param Quadrature A numeric vector of length two giving the number of
#'   quadrature points and the symmetric bound (absolute value) of the points.
#'   For example, `c(49, 6)` specifies 49 evenly spaced points from -6 to 6.
#'   These points are used to compute each examinee's posterior ability
#'   distribution when `weights = NULL`. Default is `c(49, 6)`.
#' @param weights A two-column matrix or data frame of quadrature points (first
#'   column) and weights (second column) defining the prior ability
#'   distribution used to compute the posterior distributions. It can be created
#'   with [irtQ::gen.weight()]. If `NULL`, a normal prior defined by
#'   `Quadrature`, `group.mean`, and `group.var` is used. Default is `NULL`.
#' @param group.mean A numeric value specifying the mean of the normal prior
#'   ability distribution when `weights = NULL`. Default is 0.
#' @param group.var A positive numeric value specifying the variance of the
#'   normal prior ability distribution when `weights = NULL`. Default is 1.
#' @param crit.val A numeric critical value for the pseudo-count \eqn{D^{2}}
#'   statistic. Items whose statistic is greater than or equal to `crit.val` are
#'   flagged. Because the statistic has no known null distribution, the critical
#'   value is usually obtained by bootstrapping (see Examples). If `NULL`
#'   (default), the statistic is reported but no item is flagged and
#'   purification is not performed.
#' @param min.resp A single nonnegative number, usually a positive integer,
#'   specifying the minimum number of responses an item must have to be
#'   analyzed. All responses of items with fewer than
#'   `min.resp` (but at least one) responses are set to `NA`, so their `pcd2`
#'   value is `NaN` and `N` is 0, and a warning lists those items. If `NULL`, no
#'   minimum is applied. Default is `NULL`.
#' @param purify Logical. Whether to apply the iterative purification
#'   procedure. Purification is performed only when `crit.val` is also supplied.
#'   Default is `FALSE`.
#' @param max.iter A positive integer specifying the maximum number of
#'   purification iterations. Default is 10. If the limit is reached while
#'   flagged items remain, a warning is issued, the items flagged at the last
#'   iteration are added to the flagged set, and the `complete` element of
#'   `with_purify` is `FALSE`. Because one item is removed per iteration,
#'   `max.iter` should be at least the number of items expected to drift.
#'
#' @details
#' The pseudo-count \eqn{D^{2}} statistic quantifies item parameter drift (IPD)
#' by computing, across quadrature points, the weighted squared differences
#' between the observed pseudo-proportion of correct responses and the
#' proportion expected under the item bank parameters, weighted by the
#' pseudo-counts. The observed pseudo-counts are obtained from the posterior
#' distribution of each examinee's ability, and the expected proportions are the
#' probabilities of a correct response computed from the item bank parameters
#' (Stone, 2000).
#'
#' The pseudo-count \eqn{D^{2}} statistic is calculated as:
#' \deqn{
#' D^{2}_{PC} = \sum_{k=1}^{Q} \left( \frac{r_{0k} + r_{1k}}{N}\right)
#' \left( \frac{r_{1k}}{r_{0k} + r_{1k}} - E_{1k} \right)^2
#' }
#'
#' where \eqn{Q} is the number of quadrature points, \eqn{r_{0k}} and
#' \eqn{r_{1k}} are the pseudo-counts of incorrect and correct responses at
#' quadrature point \eqn{k} (the sums of the examinees' posterior probabilities
#' at that point over the examinees with each response), \eqn{E_{1k}} is the
#' probability of a correct response at point \eqn{k} under the item bank
#' parameters, and \eqn{N} is the number of examinees who responded to the item.
#'
#' \strong{Critical Value (`crit.val`)}:
#' The `crit.val` argument specifies the threshold used to flag an item as
#' exhibiting potential parameter drift. If an item's pseudo-count \eqn{D^{2}}
#' value is greater than or equal to this threshold, the item is flagged. If
#' `crit.val = NULL`, the function reports the raw statistic without flagging.
#'
#' \strong{Minimum Response Count (`min.resp`)}:
#' The `min.resp` argument sets a lower bound on the number of responses required
#' for an item to be included in the analysis. Items with fewer responses than
#' `min.resp` are automatically excluded by replacing all their responses
#' with `NA`. This avoids unreliable estimates based on small sample sizes.
#'
#' \strong{Purification Procedure}:
#' Although Cappaert et al. (2018) did not incorporate purification into their
#' method, [irtQ::pcd2()] implements an optional iterative purification process
#' similar to Lim et al. (2022). When `purify = TRUE` and a `crit.val` is
#' provided:
#' \itemize{
#'   \item The procedure begins by identifying items flagged for drift using the initial
#'   pseudo-count \eqn{D^{2}} statistics.
#'   \item In each subsequent iteration, the item with the highest flagged pseudo-count
#'   \eqn{D^{2}} value is removed from the item set, and the statistics are recalculated
#'   using only the remaining items.
#'   \item The process continues until no additional items are flagged or the
#'   number of iterations reaches `max.iter`.
#'   \item If `max.iter` is reached before convergence, a warning is issued, the
#'   items flagged at the last iteration are added to the flagged set, and
#'   `complete` is set to `FALSE`.
#'   \item If every item is flagged, the procedure stops with a warning,
#'   because no item is left to analyze, all items are returned as flagged, and
#'   `complete` is `FALSE`.
#'   \item All flagged items and statistics are saved, and convergence status is reported.
#' }
#' This process reduces the influence of already-flagged items on the posterior
#' distributions used for the remaining items.
#'
#' @return This function returns a list with the following elements:
#'
#' \item{no_purify}{A list containing the results of pseudo-count \eqn{D^{2}} analysis
#' without applying the purification procedure. It includes:
#'   \describe{
#'     \item{ipd_stat}{A data frame summarizing the pseudo-count \eqn{D^{2}} statistics
#'     for all items. The columns include:
#'     `id` (item ID),
#'     `pcd2` (the computed \eqn{D^{2}} value), and
#'     `N` (the number of valid examinee responses per item).}
#'     \item{ipd_item}{A numeric vector of item indices that were flagged as exhibiting
#'     item parameter drift (IPD), based on the specified critical value `crit.val`.
#'     If no items are flagged or `crit.val = NULL`, this is `NULL`.}
#'   }
#' }
#'
#' \item{purify}{The value of the `purify` argument. Purification is actually
#' performed only when `crit.val` is also supplied.}
#'
#' \item{with_purify}{A list containing the results of pseudo-count \eqn{D^{2}} analysis
#' after applying the purification procedure. This list is populated only when both
#' `purify = TRUE` and `crit.val` is not `NULL`; otherwise, all of its elements are `NULL`. It
#' includes:
#'   \describe{
#'     \item{ipd_stat}{A data frame of the pseudo-count \eqn{D^{2}} statistics for all
#'     items, taken for each item from the analysis given in `n.iter`. Columns
#'     include:
#'     `id` (item ID),
#'     `pcd2` (the computed \eqn{D^{2}} value),
#'     `N` (the number of valid responses), and
#'     `n.iter` (for a flagged item, the iteration at which it was flagged, where
#'     0 is the initial analysis; for the other items, the final iteration).}
#'     \item{ipd_item}{A numeric vector of item indices (row positions in `x`)
#'     flagged as IPD items during purification, sorted in ascending order. The
#'     iteration at which each item was flagged is given in the `n.iter` column of
#'     `ipd_stat`. It contains all items if every item is flagged. It is `NULL`
#'     when no item is flagged in the initial analysis.}
#'     \item{n.iter}{An integer indicating the number of purification iterations
#'     completed. The iteration that removes the last item is not counted,
#'     because the statistics are not recomputed.}
#'     \item{complete}{A logical value indicating whether the purification procedure
#'     converged before reaching the maximum number of iterations (`max.iter`).
#'     `FALSE` if the iteration limit was reached before convergence or if every
#'     item was flagged.}
#'   }
#' }
#'
#' \item{crit.val}{A numeric value indicating the critical threshold used to flag
#' items for parameter drift. If not specified by the user, this will be `NULL`.}
#'
#' \item{call}{The matched function call.}
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [irtQ::ripd()], [irtQ::simIPD], [irtQ::gen.weight()]
#'
#' @references
#' Cappaert, K. J., Wen, Y., & Chang, Y. F. (2018). Evaluating CAT-adjusted
#' approaches for suspected item parameter drift detection. *Measurement:
#' Interdisciplinary Research and Perspectives, 16*(4), 226-238.
#' \doi{10.1080/15366367.2018.1511199}.
#'
#'   Lim, H., Choe, E. M., & Han, K. T. (2022). A residual-based differential
#'   item functioning detection framework in item response theory. *Journal of
#'   Educational Measurement, 59*(1), 80-104. \doi{10.1111/jedm.12313}.
#'
#'   Lim, H., & Han, K. T. (2026). IRT residual-based approach to detecting item
#'   parameter drift in CAT. *Journal of Educational and Behavioral Statistics*.
#'   \doi{10.3102/10769986261460852}.
#'
#'   Stone, C. A. (2000). Monte Carlo based null distribution for an alternative
#'   goodness-of-fit test statistic in IRT models. *Journal of Educational
#'   Measurement, 37*(1), 58-75. \doi{10.1111/j.1745-3984.2000.tb01076.x}.
#'
#' @examples
#' ## Example 1: No critical value specified
#' ## Compute the pseudo-count D2 statistics for dichotomous items
#' ## Import the "-prm.txt" output file generated by flexMIRT
#' flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
#'
#' # Extract metadata for the first 30 3PLM items
#' x <- bring.flexmirt(file = flex_sam, "par")$Group1$full_df[1:30, 1:6]
#'
#' # Generate abilities for 500 examinees from N(0, 1)
#' set.seed(25)
#' score <- rnorm(500, mean = 0, sd = 1)
#'
#' # Simulate response data using the item metadata and ability values
#' data <- simdat(x = x, theta = score, D = 1)
#'
#' # Compute the pseudo-count D2 statistics (no purification applied)
#' ps_d2 <- pcd2(x = x, data = data)
#' print(ps_d2)
#'
#' ## Example 2: Applying a critical value with purification
#' # crit.val = 0.002 is an arbitrary value for illustration; no drift was
#' # simulated, so any flagged item is a false positive.
#' # Compute the pseudo-count D2 statistics with purification enabled
#' ps_d2_puri <- pcd2(x = x, data = data, crit.val = 0.002, purify = TRUE)
#' print(ps_d2_puri)
#'
#' \donttest{
#' ## -- Example 3: CAT-based IPD detection using simIPD --------------------------
#' ##
#' ## The pseudo-count D2 statistic has no closed-form null distribution.
#' ## Following Cappaert et al. (2018), the critical value is estimated
#' ## empirically via bootstrap:
#' ##   (1) Select items that are not suspected of IPD (anchor items).
#' ##   (2) Compute D2 for the whole item pool from all responses, and take
#' ##       the D2 values of the anchor items as the null distribution.
#' ##   (3) Resample the null distribution with replacement, take the 95th
#' ##       percentile of each bootstrap sample, and average the percentiles
#' ##       to obtain the critical value.
#' ##   (4) Run pcd2() with this critical value.
#' ##
#' ## The simIPD dataset provides:
#' ##   item_par  : original 360-item pool (3PLM, irtQ format)
#' ##   foc_resp  : focal group CAT responses (N = 3,000; TL = 30; ~92% NA)
#' ##   key_item  : indices of 90 highly-exposed key items (IPD analysis targets)
#' ##   item.skip : indices of 270 non-key items (excluded from analysis)
#' ##   ipd_item  : 18 truly drifted items (ground truth; a & b each -0.5)
#' ## -----------------------------------------------------------------------------
#'
#' data(simIPD)
#'
#' ## -- Step 1. Select anchor items for bootstrap -----------------------------
#' ## Select items that (a) have at least min_resp_anchor observed responses
#' ## and (b) are not suspected of IPD. In practice, exclude the items you
#' ## suspect to have drifted; here the ground truth (simIPD$ipd_item) is
#' ## excluded for illustration.
#' ## The minimum response count follows Lim & Han (2026); it is a practical
#' ## addition to the procedure of Cappaert et al. (2018), who did not set such
#' ## a criterion. In a CAT pool with uneven exposure, items with few responses
#' ## have large D2 values that inflate the tail of the null distribution. For
#' ## the simIPD data, the critical value computed below is about 0.034 when
#' ## all items with at least one response are used, about 0.008 with at least
#' ## 50 responses, about 0.006 with at least 100, and about 0.004 with at
#' ## least 300.
#' min_resp_anchor <- 300
#' n_resp          <- colSums(!is.na(simIPD$foc_resp))   # observed responses per item
#' anchor_items    <- setdiff(which(n_resp >= min_resp_anchor), simIPD$ipd_item)
#'
#' ## -- Step 2. Compute D2 for the whole pool and keep the anchor items -------
#' ## D2 is computed once from all items and responses, so that the posterior
#' ## ability distributions are the same as in the final analysis (Step 4).
#' ## A warning lists the pool items that were never administered.
#' pcd2_all  <- pcd2(
#'   x    = simIPD$item_par,
#'   data = simIPD$foc_resp,
#'   D    = 1.7
#' )
#' pcd2_null <- pcd2_all$no_purify$ipd_stat$pcd2[anchor_items]
#'
#' ## -- Step 3. Bootstrap critical value (Cappaert et al., 2018) --------------
#' ## For each bootstrap iteration: resample the null D2 values with replacement
#' ## (the sample size equals the number of anchor items) and take the 95th
#' ## percentile. The critical value is the mean of these percentiles.
#' ## n_boots = 500: reduced for demo speed; use 10,000 in practice for a
#' ##   stable critical value estimate.
#' set.seed(2024)
#' n_boots <- 500
#'
#' crit_val <- mean(vapply(
#'   X   = seq_len(n_boots),
#'   FUN = function(b) {
#'     stats::quantile(sample(x = pcd2_null, size = length(pcd2_null), replace = TRUE),
#'                     probs = 0.95)
#'   },
#'   FUN.VALUE = numeric(1)
#' ))
#'
#' cat("Bootstrap critical value (alpha = 0.05):", round(crit_val, 6), "\n")
#'
#' ## -- Step 4. Run PCD2 with the bootstrap critical value + purification --
#' ## The statistics of the first analysis (pcd2_result$no_purify) are the same
#' ## as those of pcd2_all for the analyzed items, since the same responses are
#' ## used. The same warning about the items that were never administered
#' ## appears.
#' pcd2_result <- pcd2(
#'   x         = simIPD$item_par,
#'   data      = simIPD$foc_resp,
#'   D         = 1.7,
#'   item.skip = simIPD$item.skip,
#'   crit.val  = crit_val,
#'   purify    = TRUE,
#'   max.iter  = 30
#' )
#'
#' ## -- Step 5. Compare detected items to ground truth --------------------
#' detected <- pcd2_result$with_purify$ipd_item
#' cat("Truly drifted items (ground truth):", simIPD$ipd_item, "\n")
#' cat("PCD2-detected items:               ", detected, "\n")
#' cat("True positives:", sum(detected %in% simIPD$ipd_item), "of",
#'     length(simIPD$ipd_item), "\n")
#' }
#'
#'
#' @import dplyr
#' @export
pcd2 <- function(x,
                 data,
                 D = 1,
                 item.skip = NULL,
                 missing = NA,
                 Quadrature = c(49, 6.0),
                 weights = NULL,
                 group.mean = 0.0,
                 group.var = 1.0,
                 crit.val = NULL,
                 min.resp = NULL,
                 purify = FALSE,
                 max.iter = 10,
                 verbose = TRUE
) {

  # match.call
  cl <- match.call()

  # Transform a data set to matrix
  data <- as.matrix(data)

  # Confirm and correct all item metadata information
  x <- confirm_df(x)

  # Stop when the model includes any polytomous model
  if (any(x$model %in% c("GRM", "GPCM")) | any(x$cats > 2)) {
    stop("The current version only supports dichotomous response data.", call. = FALSE)
  }

  # Recode missing values
  if (!is.na(missing)) {
    data[data == missing] <- NA
  }

  # Transform the response data to a numeric matrix and check the responses
  data <- resp_to_matrix(data, x$cats, x$id)

  # Check the positions of the items to be skipped
  item.skip <- check_item_skip(item.skip, nrow(x))

  # Stop when the minimum number of responses is not a single nonnegative number
  check_min_resp(min.resp)

  # Stop when the maximum number of iterations is not a positive whole number
  check_max_iter(max.iter, purify)

  # Stop when the critical value is not a single nonnegative number
  check_crit_val(crit.val)

  # count the number of item responses across all items
  n.resp <- Rfast::colsums(!is.na(data))

  # check the items which have all missing responses
  loc_allmiss <- which(n.resp == 0L)
  if (length(loc_allmiss) > 0L) {
    memo1 <- paste0("Item(s) ", paste0(loc_allmiss, collapse = ", "), " have no response data.")
    warning(memo1, call. = FALSE)
  }

  # If min.resp is not NULL, find the items which have the number of responses
  # less than the specified value (e.g., 5), and then, replace all responses of
  # those items with NA
  if (!is.null(min.resp)) {
    loc_less <- which(n.resp < min.resp & n.resp > 0)

    # Warn only when at least one item is excluded
    if (length(loc_less) > 0L) {
      data[, loc_less] <- NA
      memo2 <- paste0("Item(s) ", paste0(loc_less, collapse = ", "), " are not analyzed ",
                      "because they do not meet the minimum response count criterion.")
      warning(memo2, call. = FALSE)
    }
  }

  # (1) when no purification is implemented
  # Do only one iteration of IPD analysis
  ipd_rst <-
    pcd2_one(x = x, data = data, D = D, item.skip = item.skip,
             Quadrature = Quadrature, weights = weights, group.mean = group.mean,
             group.var = group.var, crit.val = crit.val)

  # Create two empty lists to contain the results
  no_purify <- list(ipd_stat = NULL, ipd_item = NULL)
  with_purify <-
    list(
      ipd_stat = NULL, ipd_item = NULL, n.iter = NULL, complete = NULL
    )

  # Record the first IPD detection results into the no purification list
  no_purify$ipd_stat <- ipd_rst$ipd_stat
  if(!is.null(ipd_rst$ipd_item)) {no_purify$ipd_item <- ipd_rst$ipd_item}

  # (2) When purification is implemented
  if (!is.null(crit.val) && purify) {

    # Create an empty vector and empty data frames
    # to contain the detected IPD items and statistics
    ipd_item <- NULL
    ipd_stat <-
      data.frame(id = rep(NA_character_, nrow(x)), pcd2 = NA, N = NA,
                 n.iter = NA, stringsAsFactors = FALSE)

    # Extract the first IPD analysis results
    # and check if at least one IPD item is detected
    ipd_item_tmp <- ipd_rst$ipd_item
    ipd_stat_tmp <- ipd_rst$ipd_stat

    # Copy the response data and item meta data
    x_puri <- x
    data_puri <- data

    # Start the iteration if any item is detected as an IPD item
    if (!is.null(ipd_item_tmp)) {

      # Record unique item numbers
      item_num <- 1:nrow(x)

      # Indicate whether the purification stops because every item is flagged
      all_flagged <- FALSE

      # Print a message
      if (verbose) {
        cat("Purification started...", "\n")
      }

      for (i in 1:max.iter) {

        # Print a message
        if (verbose) {
          cat("\r", paste0("Iteration: ", i))
        }

        # Find the flagged item with the largest IPD statistic
        flag_max <- which.max(ipd_stat_tmp$pcd2)

        # Check an item that is deleted
        del_item <- item_num[flag_max]

        # Add the deleted item as the IPD item
        ipd_item <- c(ipd_item, del_item)

        # Add the IPD statistics for the detected IPD item
        ipd_stat[del_item, 1:3] <- ipd_stat_tmp[flag_max, ]
        ipd_stat[del_item, 4] <- i - 1

        # Refine the leftover items
        item_num <- item_num[-flag_max]

        # Stop the purification when no item is left to analyze
        if (length(item_num) == 0L) {
          all_flagged <- TRUE
          break
        }

        # Remove the detected IPD item data which has the largest statistic from the item metadata
        x_puri <- x_puri[-flag_max, , drop = FALSE]

        # Remove the detected IPD item data which has the largest statistic from the response data
        data_puri <- data_puri[, -flag_max, drop = FALSE]

        # Update the locations of the items that should be skipped in the purified data
        if (!is.null(item.skip)) {
          item.skip.puri <- c(1:length(item_num))[item_num %in% item.skip]
        } else {
          item.skip.puri <- NULL
        }

        # Conduct IPD analysis using the purified data
        ipd_rst_tmp <-
          pcd2_one(x = x_puri, data = data_puri, D = D, item.skip = item.skip.puri,
                   Quadrature = Quadrature, weights = weights,
                   group.mean = group.mean, group.var = group.var,
                   crit.val = crit.val)

        # Extract the IPD analysis results
        # and check if at least one IPD item is detected
        ipd_item_tmp <- ipd_rst_tmp$ipd_item
        ipd_stat_tmp <- ipd_rst_tmp$ipd_stat

        # Check if a further IPD item is flagged
        if (is.null(ipd_item_tmp)) {

          # Add the IPD statistics for rest of items
          ipd_stat[item_num, 1:3] <- ipd_stat_tmp
          ipd_stat[item_num, 4] <- i

          break
        }
      }

      # Print a message
      if (verbose) {
        cat("", "\n")
      }

      # Record the actual number of iterations; the statistics are not recomputed in the
      # iteration that removes the last item
      n_iter <- if (all_flagged) i - 1L else i

      # If every item is flagged, then, return a warning message
      if (all_flagged) {
        warning("All items were flagged during the purification, so the purification stopped ",
                "with no item left to analyze.", call. = FALSE)
        complete <- FALSE
      } else if (max.iter == n_iter & !is.null(ipd_item_tmp)) {
        # If the maximum number of iterations is reached before purification is complete,
        # then, return a warning message
        warning("The maximum number of iterations was reached before purification was completed.",
                call. = FALSE)
        complete <- FALSE

        # add flagged IPD item at the last iteration
        ipd_item <- c(ipd_item, item_num[ipd_item_tmp])

        # add the IPD statistics for rest of items
        ipd_stat[item_num, 1:3] <- ipd_stat_tmp
        ipd_stat[item_num, 4] <- i
      } else {
        complete <- TRUE

        # print a message
        if (verbose) {
          cat("Purification is finished.", "\n")
        }
      }

      # Record the final IPD detection results with the purification procedure
      with_purify$ipd_stat <- ipd_stat
      with_purify$ipd_item <- sort(ipd_item)
      with_purify$n.iter <- n_iter
      with_purify$complete <- complete

    } else {

      # In case when no IPD item is detected from the first IPD analysis results
      with_purify$ipd_stat <- cbind(no_purify$ipd_stat, n.iter = 0)
      with_purify$n.iter <- 0
      with_purify$complete <- TRUE

    }
  }

  # Summarize the results
  rst <- list(no_purify = no_purify, purify = purify,
              with_purify = with_purify, crit.val = crit.val)

  # return the IPD detection results
  rst$call <- cl
  rst

}


# Pseudo-count D2 method for a single iteration
pcd2_one <- function(x, data, D = 1, item.skip = NULL,
                     Quadrature = c(49, 6.0), weights = NULL,
                     group.mean = 0.0, group.var = 1.0, crit.val = NULL) {

  # confirm and correct all item metadata information
  x <- confirm_df(x)

  # extract item ids for all items
  id <- x$id

  # extract score categories for all items
  cats <- x$cats

  # count the total number of item in the response data set
  nitem <- ncol(data)

  # check whether the data have the same number of items as indicated in the x argument
  if (nrow(x) != nitem) {
    stop("The number of items included in 'x' and 'data' must be the same.", call. = FALSE)
  }

  # count the number of item responses across all items
  n.resp <- Rfast::colsums(!is.na(data))

  # create the weights of the prior ability distribution when it is not specified
  if (is.null(weights)) {
    # create a vector of quad-points
    quadpt <- seq(-Quadrature[2], Quadrature[2], length.out = Quadrature[1])

    # create a two column data frame to contain the quad-points and weights
    weights <- gen.weight(dist = "norm", mu = group.mean, sigma = sqrt(group.var), theta = quadpt)
  }

  # build the per-item one-hot frequency-category list used downstream
  # by divide_data() and the pseudo-count D^2 statistic (see build_freqcat())
  freq.cat <- build_freqcat(data, cats)

  # break down the item metadata into several elements
  elm_item <- breakdown(x)

  # classify the items into DRM and PRM item groups
  idx.item <- idxfinder(elm_item)
  idx.drm <- idx.item$idx.drm
  idx.prm <- idx.item$idx.prm

  # divide the data set for the mixed-item format
  datlist <- divide_data(data = data, idx.item = idx.item, freq.cat = freq.cat)
  data_drm <- cbind(datlist$data_drm_q, datlist$data_drm_p)
  data_prm <- datlist$data_prm
  data_all <- datlist$data_all

  # delete 'datlist' object
  rm(datlist, envir = environment(), inherits = FALSE)

  # find the columns of the frequency matrix corresponding to all items
  cols.item <- cols4item(nitem, cats, NULL)

  # compute the likelihood matrix
  lmat <-
    likelihood(elm_item,
               idx.drm = idx.drm, idx.prm = idx.prm,
               data_drm = data_drm, data_prm = data_prm, theta = weights[, 1], D = D)$L

  # compute posterior distribution for all examinees
  post_dist <- posterior(likehd = lmat, weights = weights, idx.std = NULL)

  # compute the pseudo-counts (posterior-weighted response counts) of score categories
  # across all items
  freq_exp_all <-
    base::as.matrix(Matrix::crossprod(post_dist, data_all))

  # split the pseudo-counts by item
  freq_exp_item <-
    purrr::map(.x = cols.item$cols.all,
               .f = ~ {freq_exp_all[, .x]})

  # total pseudo-count at each quadrature point for each item
  freq_exp_sum <-
    purrr::map(.x = freq_exp_item,
               .f = ~ {Rfast::rowsums(.x)})

  # compute the model-implied category probabilities at each quadrature point
  exp_prob_item <-
    trace(elm_item = elm_item, theta = weights[, 1], D = D,
          tcc = FALSE)$prob.cats

  # compute the pseudo-count D2 values across all items
  pc_d2 <-
    purrr::pmap_dbl(.l = list(x = freq_exp_item, y = freq_exp_sum,
                              z = n.resp, h = exp_prob_item),
                    .f = function(x, y, z, h) {
                      # contribution of each quadrature point
                      term <- (y / z) * ((x[, 2] / y) - h[, 2])^2

                      # a point without posterior mass adds nothing (this avoids 0 / 0),
                      # while an item without responses (z = 0) keeps its NaN value
                      sum(term[y > 0 | z == 0])
                    })

  # create a data frame to contain the computed pcd2 statistics
  stat_df <- data.frame(id = id, pcd2 = pc_d2, N = n.resp)

  # when there are items that should be skipped for the IPD analysis
  # insert NAs to the corresponding results of the items
  if (!is.null(item.skip)) {
    stat_df[item.skip, 2] <- NA
    pc_d2[item.skip] <- NA
  }

  # find the flagged items based on the provided critical value
  if(is.null(crit.val)) {
    ipd_item <- NULL
  } else {
    ipd_item <- as.numeric(which(pc_d2 >= crit.val))
    if(length(ipd_item) == 0) {ipd_item <- NULL}
  }

  # return the results
  rst <- list(ipd_stat = stat_df, ipd_item = ipd_item, crit.val = crit.val)
  rst

}
