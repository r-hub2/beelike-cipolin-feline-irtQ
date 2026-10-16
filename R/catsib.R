#' CATSIB DIF Detection Procedure
#'
#' This function performs DIF analysis on items using the CATSIB procedure
#' (Nandakumar & Roussos, 2004), a modified version of SIBTEST (Shealy & Stout,
#' 1993). The CATSIB procedure is suitable for computerized adaptive testing
#' (CAT) environments. In CATSIB, examinees are matched on IRT-based ability
#' estimates that have been adjusted using a regression correction method
#' (Shealy & Stout, 1993) to reduce statistical bias in the CATSIB statistic
#' caused by impact.
#'
#' @inheritParams rdif
#' @param x A data frame containing item metadata (e.g., item parameters, number
#'   of categories, IRT model types, etc.). See [irtQ::est_irt()] or
#'   [irtQ::simdat()] for more details about the item metadata. This data frame
#'   can be easily created using the [irtQ::shape_df()] function. `x` is needed
#'   to estimate the abilities internally, that is, when `score = NULL` or
#'   `purify = TRUE`. It can be `NULL` when `score` and `se` are supplied and
#'   `purify = FALSE`; the items are then labeled `item.1`, `item.2`, and so on.
#'   Default is `NULL`.
#' @param data A matrix of examinees' scored responses (0 or 1). Rows represent
#'   examinees and columns represent items, in the order of the rows of `x` when
#'   `x` is supplied.
#' @param score A numeric vector of examinees' ability estimates (theta values),
#'   in the same row order as `data`, with their standard errors in `se`. The
#'   estimates should be based on the item parameters in `x`. If `NULL`,
#'   abilities and standard errors are estimated internally with
#'   [irtQ::est_score()] using `method`, `range`, `norm.prior`, `nquad`,
#'   `weights`, and `ncore`, and `x` must be supplied. When `purify = TRUE`,
#'   abilities and standard errors are re-estimated internally at every
#'   purification iteration, so a supplied `score` is used only in the initial
#'   analysis. A missing value in `score` excludes the examinee from the
#'   analysis, also during purification, with a warning. Default is `NULL`.
#' @param se A numeric vector of the standard errors of the ability estimates in
#'   `score`, in the same order. The standard errors are used to estimate the
#'   reliability of the ability estimates in the regression correction (see
#'   **Details**). It must be supplied when `score` is supplied, and it is
#'   ignored when `score = NULL`. Default is `NULL`.
#' @param D A scaling factor in IRT models, used only when abilities are
#'   estimated internally. Default is `1`.
#' @param n.bin A numeric vector of two positive integers giving the maximum and
#'   minimum numbers of equal-width bins (intervals) of the corrected ability
#'   scale, in this order (the maximum must not be smaller than the minimum).
#'   For each item, the number of bins starts at the maximum and
#'   decreases by one until the condition set by `max.del` is met; if it is not
#'   met with the minimum number, the minimum number is used. Default is
#'   `c(80, 10)`. See **Details**.
#' @param min.binsize An integer of 2 or more specifying the minimum number of
#'   examinees of each group in a bin. A bin is used in the computation of
#'   \eqn{\hat{\beta}} only if it contains at least `min.binsize` examinees of
#'   both the reference and focal groups who responded to the item. Default is
#'   `3`, the value used by Nandakumar and Roussos (2004). See **Details**.
#' @param max.del A numeric value specifying the maximum proportion of the
#'   examinees of each group who responded to the item that may be excluded
#'   because they fall in bins that are not used (see `min.binsize`). The number
#'   of bins is reduced until the excluded proportion is at most `max.del` in
#'   both groups. It must be at least 0 and less than 1. Default is `0.075`.
#'   See **Details**.
#' @param weight.group A character string specifying the target ability
#'   distribution used to weight the bins in \eqn{\hat{\beta}} and its standard
#'   error: `"comb"` for the combined distribution of the reference and focal
#'   groups, `"foc"` for the focal group, or `"ref"` for the reference group. The
#'   weight of a bin is the proportion of the examinees of the target
#'   distribution in the bins used who fall in that bin. Default is `"comb"`.
#' @param item.skip A numeric vector of item positions (columns of `data`, or
#'   rows of `x`), or a logical vector with one value per item, identifying the
#'   items to exclude from the analysis. Skipped items still contribute to the
#'   ability estimates, but their results (`beta`, `se`, `z.beta`, `p`, and the
#'   contingency table) are set to `NA`, and they are never flagged or removed
#'   during purification; their sample sizes (`n.ref`, `n.foc`, and `n.total`)
#'   are still reported, as in [rdif()] and [crdif()]. If `NULL`, all items are
#'   analyzed. Default is `NULL`.
#' @param alpha A numeric value specifying the significance level (\eqn{\alpha})
#'   of the two-sided test of the CATSIB statistic. An item is flagged when its
#'   p-value is less than or equal to `alpha`. The comparison uses unrounded
#'   p-values, although the p-values reported in `dif_stat` are rounded to four
#'   decimal places. Default is `0.05`.
#' @param min.resp A single nonnegative number, usually a positive integer,
#'   specifying the minimum number of item responses that an examinee must have
#'   to be used in the analysis. All responses of examinees with fewer than
#'   `min.resp` (but at least one) responses are set to `NA` before the ability
#'   estimation, in the initial analysis and at every purification iteration,
#'   also when `score` is supplied, and these examinees are excluded from the
#'   analysis. A warning reports the number of examinees excluded in the initial
#'   analysis, and one warning at the end of the purification reports the
#'   examinees excluded during it. A large `min.resp` can leave a group with
#'   fewer than two ability estimates within `range` in a later purification
#'   iteration, which stops the purification with a warning. If `NULL`, no
#'   minimum is applied. Default is `NULL`. See **Details** for more
#'   information.
#' @param method A character string indicating the scoring method used to
#'   estimate the abilities when `score = NULL` and in every purification
#'   iteration. Available options are:
#'   - `"ML"`: Maximum likelihood estimation
#'   - `"MLF"`: Maximum likelihood estimation with fences
#'   - `"WL"`: Weighted likelihood estimation (Warm, 1989)
#'   - `"MAP"`: Maximum a posteriori estimation (Hambleton et al., 1991)
#'   - `"EAP"`: Expected a posteriori estimation (Bock & Mislevy, 1982)
#'
#'   The methods `"EAP.SUM"` and `"INV.TCC"` of [irtQ::est_score()] are not
#'   supported, because they do not give an ability estimate for each response
#'   pattern. `"ML"` or `"WL"` is recommended, because the regression
#'   correction assumes unbiased ability estimates (see **Details**). Default is
#'   `"ML"`.
#' @param range A numeric vector of length two specifying the lower and upper
#'   bounds of the ability scale. It is used when abilities are estimated
#'   internally with the scoring methods `"ML"`, `"MLF"`, `"WL"`, and `"MAP"`. In
#'   addition, ability estimates at or beyond these bounds are excluded when the
#'   mean, the variance, and the error variance of the ability estimates are
#'   computed for the regression correction (see **Details**); this also applies
#'   when `score` is supplied, so `range` should match the bounds used to obtain
#'   `score`. Default is `c(-5, 5)`.
#' @param rho.screen Logical. If `TRUE`, the ability estimates whose squared
#'   standard error is equal to or greater than the variance of the ability
#'   estimates of their group are left out of the computation of the mean, the
#'   variance, and the error variance that define \eqn{\hat{\rho}^2_G} in the
#'   regression correction. These examinees are still matched and used to
#'   compute \eqn{\hat{\beta}}. If `FALSE`, \eqn{\hat{\rho}^2_G} is computed
#'   from all ability estimates within `range`, as in Nandakumar and Roussos
#'   (2004, Equation 7). The screening is not part of the original procedure and
#'   is added in irtQ. Default is `TRUE`. See **Details**.
#' @param ... Additional arguments passed to [irtQ::est_score()] when abilities
#'   are estimated internally, for example `tol`. The `est_score()` arguments
#'   `se`, `max.iter`, and `missing` cannot be passed this way, because
#'   [irtQ::catsib()] uses arguments with these names for its own purposes.
#'
#' @details CATSIB (Nandakumar & Roussos, 2004) estimates the amount of DIF in
#'   an item as \eqn{\hat{\beta} = \sum_{k} (\hat{P}_{R,k} - \hat{P}_{F,k})
#'   \hat{p}_k}, where \eqn{\hat{P}_{R,k}} and \eqn{\hat{P}_{F,k}} are the
#'   proportions of correct responses of the reference and focal groups in bin k
#'   of the corrected ability scale, and \eqn{\hat{p}_k} is the weight of bin k
#'   (see `weight.group`). A positive \eqn{\hat{\beta}} indicates DIF in favor of
#'   the reference group. The standard error of \eqn{\hat{\beta}} is estimated
#'   from the sample variances of the item responses in each bin and group
#'   (Nandakumar & Roussos, 2004, Equation 10), and the test statistic
#'   \eqn{\hat{\beta} / \mathrm{SE}(\hat{\beta})} (`z.beta`) is compared with
#'   the standard normal distribution in a two-sided test. Examinees are matched
#'   on \eqn{\hat{\theta}^{\ast}}, an estimate of the expected value of
#'   \eqn{\theta} given \eqn{\hat{\theta}} obtained by the regression correction
#'   described below. The range of \eqn{\hat{\theta}^{\ast}} is divided into K
#'   equal-width bins, and examinees are classified into the bins by their
#'   \eqn{\hat{\theta}^{\ast}}. A bin is used only if it contains at least
#'   `min.binsize` examinees of both groups who responded to the item;
#'   Nandakumar and Roussos (2004) required at least three, which is the default.
#'
#'   To choose the number of bins K, [irtQ::catsib()] follows the rule of
#'   Nandakumar and Roussos (2004). For each item, K starts at a large value (80
#'   by default), and if more than 7.5% of the examinees of either group who
#'   responded to the item would be excluded because they fall in bins that are
#'   not used, K is reduced by one, until at most 7.5% of each group is excluded.
#'   Nandakumar and Roussos (2004) did not allow fewer than 20 bins, but they
#'   reported that a minimum of 10 bins reduced the estimation bias caused by
#'   excluding examinees from sparse bins. [irtQ::catsib()] therefore uses a
#'   minimum of 10 bins by default; if the condition is not met with the minimum
#'   number, the minimum number is used. The maximum and minimum numbers of bins
#'   are set by `n.bin`, and the maximum proportion of excluded examinees by
#'   `max.del`.
#'
#'   Li and Stout (1996) and Nandakumar and Roussos (2004) used the combined
#'   distribution of the two groups as the target ability distribution, which is
#'   the default of `weight.group`. See Nandakumar and Roussos (2004) for further
#'   details about the CATSIB method.
#'
#'   Because the ability distributions of the reference and focal groups can
#'   differ (impact), matching examinees directly on \eqn{\hat{\theta}} would
#'   bias \eqn{\hat{\beta}} toward DIF in favor of the group with the higher mean
#'   ability. CATSIB therefore applies a regression correction separately in each
#'   group G (Nandakumar & Roussos, 2004, Equation 7): \eqn{\hat{\theta}^{\ast} =
#'   \bar{\hat{\theta}}_G + \hat{\rho}^2_G (\hat{\theta} -
#'   \bar{\hat{\theta}}_G)}, where \eqn{\bar{\hat{\theta}}_G} is the mean of the
#'   ability estimates in group G and \eqn{\hat{\rho}^2_G = 1 -
#'   \hat{\sigma}^2_{e,G} / \hat{\sigma}^2_{\hat{\theta},G}} is the estimated
#'   reliability of the ability estimates in group G. Here
#'   \eqn{\hat{\sigma}^2_{e,G}} is the mean of the squared standard errors (`se`)
#'   and \eqn{\hat{\sigma}^2_{\hat{\theta},G}} is the variance of the ability
#'   estimates in group G. Ability estimates at or beyond a bound of `range` are
#'   excluded when these quantities are computed, because their standard errors
#'   are not meaningful; this also applies when `score` is supplied, so `range`
#'   should match the bounds used to obtain `score`.
#'
#'   Because \eqn{\hat{\sigma}^2_{e,G}} is a mean of squared standard errors, a
#'   few examinees with very large standard errors, typically examinees with very
#'   low or very high scores whose abilities are poorly determined, can make
#'   \eqn{\hat{\rho}^2_G} much smaller than the reliability of the other ability
#'   estimates, or even negative. The corrected ability estimates then shrink too
#'   much toward the group mean, which distorts the matching and can inflate the
#'   Type I error rate, particularly during purification. When `rho.screen =
#'   TRUE` (the default), the examinees whose squared standard error is greater
#'   than or equal to \eqn{\hat{\sigma}^2_{\hat{\theta},G}} of their group are
#'   left out of the computation of the mean, the variance, and the error
#'   variance that define \eqn{\hat{\rho}^2_G}, unless fewer than two
#'   estimates would be left. For these examinees the individual reliability
#'   \eqn{1 - SE^2 / \hat{\sigma}^2_{\hat{\theta},G}} is zero or negative, that
#'   is, their ability estimates carry no information about their abilities
#'   beyond the spread of the group. They still receive corrected ability
#'   estimates and are used in the matching and in the computation of
#'   \eqn{\hat{\beta}}. This screening is not part of the original procedure
#'   (Nandakumar & Roussos, 2004) and is added in irtQ; set `rho.screen = FALSE`
#'   to compute \eqn{\hat{\rho}^2_G} from all ability estimates within `range` as
#'   in the original procedure. Even with the screening, the Type I error rate
#'   of CATSIB can be higher than the nominal level, and in data without DIF the
#'   purification can increase the number of false positives. In addition,
#'   \eqn{\hat{\rho}^2_G} is bounded to the interval \eqn{[0.05, 1]}, and a
#'   warning is issued when its unbounded value is below 0.05. A negative value
#'   would reverse the order of the ability estimates within a group, and a value
#'   near zero would place almost all corrected ability estimates of a group in
#'   one bin. This bound is not part of the original procedure and is added in
#'   irtQ. Even with the bound, the corrected ability estimates of the two groups
#'   may not overlap when the groups differ in ability; then no bin contains
#'   enough examinees of both groups for an item, \eqn{\hat{\beta}}, its standard
#'   error, `z.beta`, and the p-value are `NA`, one warning names the items, and
#'   the items are not flagged.
#'
#'   The regression correction assumes that the ability estimates are unbiased,
#'   that is, that \eqn{\hat{\theta}} equals \eqn{\theta} plus an error with
#'   mean zero. The Bayesian estimates (`"MAP"` and `"EAP"`) are already shrunk
#'   toward the prior mean, so the correction shrinks them a second time. Use
#'   `"ML"` or `"WL"`, especially when the abilities are re-estimated during the
#'   purification.
#'
#'   CATSIB, like SIBTEST (Shealy & Stout, 1993), was designed and evaluated for
#'   uniform DIF (Nandakumar & Roussos, 2004). Because \eqn{\hat{\beta}} sums
#'   signed differences across bins, differences of opposite sign can cancel, so
#'   its sensitivity to nonuniform or mixed DIF is limited, and results for items
#'   suspected of nonuniform DIF should be interpreted with caution.
#'
#'   Nandakumar and Roussos (2004) did not propose a purification procedure for
#'   CATSIB, but [irtQ::catsib()] can apply an iterative purification procedure
#'   in the manner of Lim et al. (2022). It requires `x`. When `purify = TRUE`,
#'   at each iteration the flagged item with the smallest p-value is removed (the
#'   p-values are compared on the log scale), the abilities and their standard
#'   errors are re-estimated from the remaining items with the scoring method
#'   specified in `method`, and the CATSIB statistics of the remaining items are
#'   recomputed. Supplied `score` and `se` are therefore used only in the initial
#'   analysis. The procedure stops when no remaining item is flagged or when
#'   `max.iter` iterations have been performed. If no item is flagged in the
#'   initial analysis, no iteration is performed. If every item is flagged, the
#'   procedure stops with a warning, because no item is left to estimate the
#'   abilities, and `complete` is `FALSE`. The procedure also stops with a
#'   warning and `complete = FALSE` when, after an item is removed, a group has
#'   fewer than two ability estimates within `range`, which the regression
#'   correction requires, or no bin contains enough examinees of both groups
#'   for any remaining item; the remaining items then keep the results of the
#'   previous iteration, and the items flagged in that iteration are added to the
#'   flagged items.
#'
#'   Scoring based on a small number of item responses can lead to large
#'   standard errors, potentially reducing the accuracy of DIF detection with
#'   CATSIB. The `min.resp` argument excludes such examinees. For example, if
#'   `min.resp = 5`, examinees who responded to fewer than five items (but at
#'   least one) have all their responses treated as missing (i.e., `NA`) before
#'   the abilities are estimated, in the initial analysis and at every
#'   purification iteration, and a warning reports their number in the initial
#'   analysis. As a result, their ability estimates are missing, and they are not
#'   used in the computation of the CATSIB statistics. When `score` is supplied,
#'   the scores of these examinees are also treated as missing. During
#'   purification, the count is based on the items that remain after the flagged
#'   items are removed, and an examinee excluded once stays excluded. If
#'   `min.resp = NULL`, a score is computed for any examinee with at least one
#'   valid item response.
#'
#'   Examinees whose ability estimate is `NA`, such as examinees without any item
#'   response or examinees with `NA` in a supplied `score`, are excluded from the
#'   analysis, and a warning reports the number of examinees with item responses
#'   who are excluded because of a missing value in `score`; these examinees stay
#'   excluded during purification. During purification, an examinee whose
#'   remaining responses are all removed with the flagged items, or who has fewer
#'   than `min.resp` responses to the remaining items, has no ability estimate
#'   and is excluded, and one warning at the end reports the number of such
#'   examinees. The sample sizes `n.ref` and `n.foc` count only the examinees in
#'   the bins used to compute \eqn{\hat{\beta}}.
#'
#'   The `group` argument should contain exactly two distinct values, one for the
#'   reference group and one for the focal group identified by `focal.name`, and
#'   its length must equal the number of rows in `data`. Examinees with `NA` in
#'   `group` are excluded from the analysis. The function stops when the
#'   reference group or the focal group has fewer than two ability estimates
#'   within `range`, because the regression correction cannot be computed.
#'
#'   The print method of this class has an argument `digits` that is not used;
#'   the statistics are printed rounded to three decimal places.
#'
#' @return This function returns an object of class `"catsib"`, which is a list
#' with the following five components:
#'
#' \item{no_purify}{A list containing the results of the DIF analysis without
#' applying a purification procedure. This list includes:
#'   \describe{
#'     \item{dif_stat}{A data frame with one row per item and the columns `id`,
#'     `beta`, `se`, `z.beta`, `p`, `n.ref`, `n.foc`, and `n.total`: the item ID,
#'     the CATSIB statistic \eqn{\hat{\beta}}, its standard error, the
#'     standardized statistic, the two-sided p-value, and the numbers of
#'     examinees of the reference and focal groups in the bins used to compute
#'     \eqn{\hat{\beta}}, and their sum. A positive \eqn{\hat{\beta}} indicates
#'     DIF in favor of the reference group. All values are rounded to four
#'     decimal places.}
#'     \item{dif_item}{A numeric vector of the positions (columns of `data`) of
#'     the flagged items, or `NULL` if no item is flagged.}
#'     \item{contingency}{A list with one element per item, named by item ID. Each
#'     element is a data frame with one row per bin used to compute
#'     \eqn{\hat{\beta}} and the columns `bin` (the interval of corrected ability
#'     estimates), `n.ref`, `prop.ref`, `var.ref`, `n.foc`, `prop.foc`, and
#'     `var.foc` (the number of examinees, the proportion of correct responses,
#'     and the sample variance of the item responses of each group in the bin),
#'     `n.total`, `weight` (the weight of the bin), `beta` (the weighted
#'     difference in the proportions of correct responses), and `var.beta` (the
#'     contribution of the bin to the squared standard error), followed by a row
#'     whose `bin` is `Total`, with the totals of `n.ref`, `n.foc`, `n.total`,
#'     `beta`, and `var.beta`. The element is `NA` for a skipped item.}
#'   }
#' }
#'
#' \item{purify}{A logical value indicating whether purification was requested
#' (the `purify` argument).}
#'
#' \item{with_purify}{A list with the results of the purification procedure. All
#' elements are `NULL` when `purify = FALSE`. The elements are:
#'   \describe{
#'     \item{dif_stat}{A data frame with the same columns as
#'     `no_purify$dif_stat` and an additional column `n.iter`. For an item
#'     removed during purification, the row reports the statistics from the
#'     iteration in which the item was removed, and `n.iter` is that iteration
#'     minus 1 (0 for an item removed on the basis of the initial analysis). For
#'     the other items, the row reports the statistics from the last iteration,
#'     and `n.iter` is the number of that iteration. If no item is flagged in the
#'     initial analysis, this is `no_purify$dif_stat` with `n.iter = 0`.}
#'     \item{dif_item}{A numeric vector of the positions (rows of `x`) of the
#'     items flagged by the CATSIB statistic, sorted in ascending order. It
#'     contains the items removed during purification and, if the purification
#'     stops before it is complete (see `complete`) with items left, the items
#'     flagged in the last iteration (iteration `n.iter`). It contains all items
#'     if every item is flagged. `NULL` if no item is flagged in the initial
#'     analysis.}
#'     \item{n.iter}{The number of purification iterations performed (0 if no
#'     item is flagged in the initial analysis). An iteration after which the
#'     statistics cannot be recomputed (see **Details**) is not counted.}
#'     \item{complete}{A logical value. `TRUE` if the procedure stopped because
#'     no remaining item was flagged (including the case in which no item is
#'     flagged in the initial analysis), and `FALSE` if it stopped because
#'     `max.iter` was reached, because every item was flagged, because a group
#'     had fewer than two ability estimates within `range`, or because no
#'     remaining item had a usable bin.}
#'     \item{contingency}{A list with the same structure as
#'     `no_purify$contingency`. For an item removed during purification, it holds
#'     the table from the iteration in which the item was removed; for the other
#'     items, it holds the table from the last iteration. If no item is flagged
#'     in the initial analysis, this is `no_purify$contingency`.}
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
#' @seealso [irtQ::rdif()], [irtQ::grdif()], [irtQ::est_irt()],
#'   [irtQ::est_item()], [irtQ::simdat()], [irtQ::shape_df()],
#'   [irtQ::est_score()]
#'
#' @references Bock, R. D., & Mislevy, R. J. (1982). Adaptive EAP estimation of
#'   ability in a microcomputer environment. *Applied Psychological Measurement,
#'   6*(4), 431-444. \doi{10.1177/014662168200600405}.
#'
#'   Hambleton, R. K., Swaminathan, H., & Rogers, H. J. (1991). *Fundamentals of
#'   item response theory*. Newbury Park, CA: Sage.
#'
#'   Li, H.-H., & Stout, W. (1996). A new procedure for detection of crossing
#'   DIF. *Psychometrika, 61*(4), 647-677. \doi{10.1007/BF02294041}.
#'
#'   Lim, H., Choe, E. M., & Han, K. T. (2022). A residual-based differential
#'   item functioning detection framework in item response theory. *Journal of
#'   Educational Measurement, 59*(1), 80-104. \doi{10.1111/jedm.12313}.
#'
#'   Nandakumar, R., & Roussos, L. (2004). Evaluation of the CATSIB DIF
#'   procedure in a pretest setting. *Journal of Educational and Behavioral
#'   Statistics, 29*(2), 177-199. \doi{10.3102/10769986029002177}.
#'
#'   Shealy, R., & Stout, W. (1993). A model-based standardization approach that
#'   separates true bias/DIF from group ability differences and detects test
#'   bias/DTF as well as item bias/DIF. *Psychometrika, 58*(2), 159-194.
#'   \doi{10.1007/BF02294572}.
#'
#'   Warm, T. A. (1989). Weighted likelihood estimation of ability in item
#'   response theory. *Psychometrika, 54*(3), 427-450.
#'   \doi{10.1007/BF02294627}.
#'
#'
#' @examples
#' \donttest{
#' # Load required package
#' library("dplyr")
#'
#' ## Uniform DIF Detection
#' ###############################################
#' # (1) Simulate data with true uniform DIF
#' ###############################################
#'
#' # Import the "-prm.txt" output file from flexMIRT
#' flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
#'
#' # Select 36 3PLM items that are non-DIF
#' par_nstd <-
#'   bring.flexmirt(file = flex_sam, "par")$Group1$full_df %>%
#'   dplyr::filter(.data$model == "3PLM") %>%
#'   dplyr::filter(dplyr::row_number() %in% 1:36) %>%
#'   dplyr::select(1:6)
#' par_nstd$id <- paste0("nondif", 1:36)
#'
#' # Generate four new items to contain uniform DIF
#' difpar_ref <-
#'   shape_df(
#'     par.drm = list(a = c(0.8, 1.5, 0.8, 1.5), b = c(0.0, 0.0, -0.5, -0.5), g = rep(0.15, 4)),
#'     item.id = paste0("dif", 1:4), cats = 2, model = "3PLM"
#'   )
#'
#' # Introduce uniform DIF in the focal group by shifting b-parameters
#' difpar_foc <-
#'   difpar_ref %>%
#'   dplyr::mutate_at(.vars = "par.2", .funs = function(x) x + rep(0.7, 4))
#'
#' # Combine the 4 DIF and 36 non-DIF items for both reference and focal groups
#' # Therefore, the first four items now exhibit uniform DIF
#' par_ref <- rbind(difpar_ref, par_nstd)
#' par_foc <- rbind(difpar_foc, par_nstd)
#'
#' # Generate true theta values
#' set.seed(123)
#' theta_ref <- rnorm(500, 0.0, 1.0)
#' theta_foc <- rnorm(500, 0.0, 1.0)
#'
#' # Simulate response data
#' resp_ref <- simdat(par_ref, theta = theta_ref, D = 1)
#' resp_foc <- simdat(par_foc, theta = theta_foc, D = 1)
#' data <- rbind(resp_ref, resp_foc)
#'
#' ###############################################
#' # (2) Estimate item and ability parameters
#' #     using the aggregated data
#' ###############################################
#'
#' # Estimate item parameters
#' est_mod <- est_irt(data = data, D = 1, model = "3PLM")
#' est_par <- est_mod$par.est
#'
#' # Estimate ability parameters using ML
#' theta_est <- est_score(x = est_par, data = data, method = "ML")
#' score <- theta_est$est.theta
#' se <- theta_est$se.theta
#'
#' ###############################################
#' # (3) Conduct DIF analysis
#' ###############################################
#' # Create a vector of group membership indicators
#' # where '1' indicates the focal group
#' group <- c(rep(0, 500), rep(1, 500))
#'
#' # (1) Compute the CATSIB statistic using provided scores and standard errors,
#' #     without purification (x can be NULL in this case)
#' dif_1 <- catsib(
#'   x = NULL, data = data, D = 1, score = score, se = se, group = group, focal.name = 1,
#'   weight.group = "comb", alpha = 0.05, missing = NA, purify = FALSE
#' )
#' print(dif_1)
#'
#' # (2) Compute the CATSIB statistic with purification (x is required; not
#' #     run: it takes about 7 seconds because the abilities of the examinees
#' #     are re-estimated at every iteration)
#' #
#' #   dif_2 <- catsib(
#' #     x = est_par, data = data, D = 1, score = score, se = se, group = group,
#' #     focal.name = 1, weight.group = "comb", alpha = 0.05, missing = NA,
#' #     purify = TRUE
#' #   )
#' #   print(dif_2)
#' }
#'
#' @import dplyr
#' @importFrom janitor adorn_totals
#' @export
catsib <- function(x = NULL,
                   data,
                   score = NULL,
                   se = NULL,
                   group,
                   focal.name,
                   item.skip = NULL,
                   D = 1,
                   n.bin = c(80, 10),
                   min.binsize = 3,
                   max.del = 0.075,
                   weight.group = c("comb", "foc", "ref"),
                   alpha = 0.05,
                   missing = NA,
                   purify = FALSE,
                   max.iter = 10,
                   min.resp = NULL,
                   method = "ML",
                   range = c(-5, 5),
                   norm.prior = c(0, 1),
                   nquad = 41,
                   weights = NULL,
                   ncore = 1,
                   verbose = TRUE,
                   rho.screen = TRUE,
                   ...) {
  # match.call
  cl <- match.call()

  ## ----------------------------------
  ## (1) prepare DIF analysis
  ## ----------------------------------
  # when purify = TRUE, the item metadata should be provided.
  # if not, stop.
  if (purify & is.null(x)) {
    stop("To implement a purification process, the item metadata must be provided in the argument 'x'.", call. = FALSE)
  }

  # clean the data frame of the item metadata
  if (!is.null(x)) {
    # confirm and correct all item metadata information
    x <- confirm_df(x)

    # stop when the model includes any polytomous model
    if (any(x$model %in% c("GRM", "GPCM")) | any(x$cats > 2)) {
      stop("catsib() supports dichotomous items only.", call. = FALSE)
    }
  }

  # stop when the abilities must be estimated but the item metadata are not provided
  if (is.null(score) && is.null(x)) {
    stop("To estimate the abilities, the item metadata must be provided in the argument 'x' ",
         "when 'score' is not supplied.", call. = FALSE)
  }

  # transform the response data to a matrix form
  data <- as.matrix(data)

  # re-code missing values
  if (!is.na(missing)) {
    data[data == missing] <- NA
  }

  # stop when the model includes any polytomous response data
  if (any(suppressWarnings(as.numeric(data)) > 1, na.rm = TRUE)) {
    stop("catsib() supports dichotomous responses (0 or 1) only.", call. = FALSE)
  }

  # stop when the ability estimates are given without their standard errors
  if (!is.null(score) && is.null(se)) {
    stop(
      "The standard errors of the ability estimates must be supplied in the argument 'se' ",
      "when 'score' is supplied.",
      call. = FALSE
    )
  }

  # create an item id
  if (is.null(x)) {
    item.id <- paste0("item.", 1:ncol(data))
  } else {
    item.id <- x$id
  }

  # transform the response data to a numeric matrix and check the responses
  data <- resp_to_matrix(data, rep(2, length(item.id)), item.id)

  # check the group vector and the focal group
  check_group(group, focal.name, nrow(data))

  # stop when the group vector has more than two groups
  if (length(unique(group[!is.na(group)])) > 2L) {
    stop("'group' must contain exactly two groups (the reference and focal groups).", call. = FALSE)
  }

  # stop when the significance level is not a single value between 0 and 1
  check_alpha(alpha)

  # stop when the maximum number of iterations is not a positive whole number
  check_max_iter(max.iter, purify)

  # stop when the minimum number of responses is not a single nonnegative number
  check_min_resp(min.resp)

  # stop when the abilities are estimated or updated with an unsupported scoring method
  if (is.null(score) || purify) {
    check_pattern_method(method, fn = "catsib()")
  }

  # check the positions of the items to be skipped
  item.skip <- check_item_skip(item.skip, length(item.id))

  # stop when the numbers of bins are not two positive whole numbers in decreasing order
  if (!(is.numeric(n.bin) && length(n.bin) == 2L && all(is.finite(n.bin)) && all(n.bin >= 1) &&
        all(n.bin == round(n.bin)) && n.bin[1] >= n.bin[2])) {
    stop("'n.bin' must contain two positive whole numbers, the maximum and the minimum numbers ",
         "of bins, in this order.", call. = FALSE)
  }

  # stop when the minimum bin size cannot give the variance of the responses in a bin
  if (!(is.numeric(min.binsize) && length(min.binsize) == 1L && is.finite(min.binsize) &&
        min.binsize >= 2 && min.binsize == round(min.binsize))) {
    stop("'min.binsize' must be a single whole number of 2 or more.", call. = FALSE)
  }

  # stop when the maximum proportion of excluded examinees is not a single proportion
  if (!(is.numeric(max.del) && length(max.del) == 1L && is.finite(max.del) && max.del >= 0 &&
        max.del < 1)) {
    stop("'max.del' must be a single number from 0 to less than 1.", call. = FALSE)
  }

  # stop when rho.screen is not a single logical value
  if (!is.logical(rho.screen) || length(rho.screen) != 1L || is.na(rho.screen)) {
    stop("'rho.screen' must be TRUE or FALSE.", call. = FALSE)
  }

  # record whether the ability estimates are supplied
  score_given <- !is.null(score)

  # if min.resp is not NULL, find the examinees who have the number of responses
  # less than specified value (e.g., 5). Then, replace their all responses with NA
  loc_less <- integer(0)
  if (!is.null(min.resp)) {
    n_resp <- rowSums(!is.na(data))
    loc_less <- which(n_resp < min.resp & n_resp > 0)

    # warn that the examinees with too few responses are excluded
    if (length(loc_less) > 0L) {
      warning(length(loc_less), " examinee(s) with fewer than ", min.resp,
              " responses were excluded.", call. = FALSE)
    }
    data[loc_less, ] <- NA
  }

  # compute the score if score = NULL
  if (!is.null(score)) {
    # transform the scores to a vector form and check them
    score <- check_score(score, nrow(data))

    # transform the standard errors to a vector form
    if (is.matrix(se) | is.data.frame(se)) {
      se <- as.numeric(data.matrix(se))
    }

    # stop when the standard errors are not numbers or do not match the rows of the response data
    if (!is.numeric(se) || length(se) != nrow(data)) {
      stop("'se' must be a numeric vector whose length equals the number of rows in 'data'.",
           call. = FALSE)
    }

    # treat the ability estimates of the examinees excluded by min.resp as missing
    score[loc_less] <- NA
    se[loc_less] <- NA
  } else {
    # estimate the abilities; the examinees without a response are excluded without a warning
    score_rst <- dif_rescore(
      x = x, data = data, score = NULL, loc_resp = NULL, first = TRUE, D = D, method = method,
      range = range, norm.prior = norm.prior, nquad = nquad, weights = weights, ncore = ncore,
      quiet_na = TRUE, return_se = TRUE, ...
    )
    score <- score_rst$score
    se <- score_rst$se
  }

  # stop when no examinee has an ability estimate
  if (all(is.na(score[!is.na(group)]))) {
    stop("No examinee has an ability estimate.", call. = FALSE)
  }

  # exclude the examinees with responses but without an ability estimate, with a warning
  na_score <- is.na(score) & rowSums(!is.na(data)) > 0
  if (any(na_score)) {
    warning(sum(na_score), " examinee(s) with a missing ability estimate were excluded.",
            call. = FALSE)
    data[na_score, ] <- NA
  }

  # stop when a group has fewer than two ability estimates to compute the regression correction
  if (!catsib_enough(score = score, group = group, focal.name = focal.name, range = range)) {
    stop("Each group must have at least two ability estimates within 'range' to compute ",
         "the regression correction.", call. = FALSE)
  }

  # confirm the target ability density function
  weight.group <- match.arg(weight.group)

  # set the maximum & minimum number of bins (intervals)
  max.bin <- n.bin[1]
  min.bin <- n.bin[2]

  # a) when no purification is set
  # do only one iteration of DIF analysis:
  # compute the beta statistic and its SE for all items using
  # corrected theta scores
  dif_rst <-
    catsib_one(
      data = data, group = group, focal.name = focal.name, item.skip = item.skip,
      score = score, se = se, range = range, item.id = item.id, max.bin = max.bin,
      min.bin = min.bin, min.binsize = min.binsize, max.del = max.del,
      weight.group = weight.group, alpha = alpha, rho.screen = rho.screen
    )

  # record the lowest reliabilities of the ability estimates of the reference and focal groups
  rho_min <- dif_rst$rho_raw

  # create two empty lists to contain the results
  no_purify <- list(dif_stat = NULL, dif_item = NULL, contingency = NULL)
  with_purify <- list(
    dif_stat = NULL, dif_item = NULL, n.iter = NULL,
    complete = NULL, contingency = NULL
  )

  # record the first DIF detection results into the no purification list
  no_purify$dif_stat <- dif_rst$dif_stat
  no_purify$dif_item <- dif_rst$dif_item
  no_purify$contingency <- dif_rst$contingency

  # when purification is used
  if (purify) {
    # create an empty vector and empty data frames
    # to contain the detected DIF items, statistics, and contingency tables
    dif_item <- NULL
    dif_stat <-
      data.frame(
        id = rep(NA_character_, nrow(x)), beta = NA, se = NA,
        z.beta = NA, p = NA, n.ref = NA, n.foc = NA, n.total = NA, n.iter = NA,
        stringsAsFactors = FALSE
      )
    contingency <- vector("list", nrow(x))
    names(contingency) <- item.id

    # extract the first DIF analysis results
    # and check if at least one DIF item is detected
    dif_item_tmp <- dif_rst$dif_item
    log_p_tmp <- dif_rst$log_p
    dif_stat_tmp <- dif_rst$dif_stat
    contingency_tmp <- dif_rst$contingency

    # copy the response data, item meta data, ability estimates, and standard errors
    x_puri <- x
    data_puri <- data
    score_puri <- score
    se_puri <- se

    # start the iteration if any item is detected as a DIF item
    if (!is.null(dif_item_tmp)) {
      # record unique item numbers
      item_num <- 1:nrow(x)

      # indicate whether the purification stops because every item is flagged
      all_flagged <- FALSE

      # indicate whether the purification stops because a group has too few ability estimates
      no_score <- FALSE

      # indicate whether the purification stops because no remaining item has a usable bin
      no_test <- FALSE

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

        # add the DIF statistics and the contingency table of the removed item
        dif_stat[del_item, 1:8] <- dif_stat_tmp[flag_del, ]
        dif_stat[del_item, 9] <- i - 1
        contingency[del_item] <- contingency_tmp[flag_del]

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
        se_prev <- se_puri

        # compute the updated ability estimates and their standard errors after deleting the
        # detected DIF item data; only the examinees who responded to the deleted item are
        # rescored after the first iteration
        score_rst_puri <-
          dif_rescore(
            x = x_puri, data = data_puri, score = score_puri, loc_resp = loc_resp,
            first = (i == 1L), D = D, method = method, range = range,
            norm.prior = norm.prior, nquad = nquad, weights = weights, ncore = ncore,
            quiet_na = TRUE, quiet_ncore = (i > 1L || !score_given), se_prev = se_puri,
            return_se = TRUE, ...
          )
        score_puri <- score_rst_puri$score
        se_puri <- score_rst_puri$se

        # stop the purification when fewer than two examinees of a group have an ability estimate
        # to compute the regression correction; the remaining items keep the results of the
        # previous iteration
        if (!catsib_enough(score = score_puri, group = group, focal.name = focal.name, range = range)) {
          no_score <- TRUE
          score_puri <- score_prev
          se_puri <- se_prev
          dif_stat[item_num, 1:8] <- dif_stat_tmp[-flag_del, ]
          dif_stat[item_num, 9] <- i - 1
          contingency[item_num] <- contingency_tmp[-flag_del]

          # add the items flagged at the previous iteration, except the item deleted at this one
          dif_item <- c(dif_item, item_num[match(setdiff(dif_item_tmp, flag_del),
                                                 seq_len(nrow(dif_stat_tmp))[-flag_del])])
          break
        }

        # do DIF analysis using the updated ability estimates
        item.id <- x_puri$id
        dif_rst_tmp <-
          catsib_one(
            data = data_puri, group = group, focal.name = focal.name,
            item.skip = item.skip.puri, score = score_puri, se = se_puri,
            range = range, item.id = item.id, max.bin = max.bin, min.bin = min.bin,
            min.binsize = min.binsize, max.del = max.del, weight.group = weight.group,
            alpha = alpha, rho.screen = rho.screen
          )

        # count the remaining items that are tested
        n_test <- length(setdiff(seq_len(nrow(x_puri)), item.skip.puri))

        # stop the purification when no remaining item has a usable bin; the remaining items keep the
        # results of the previous iteration. The skipped items are not counted, so the purification
        # ends normally when only skipped items remain
        if (n_test > 0L && length(dif_rst_tmp$no_bin) == n_test) {
          no_test <- TRUE
          score_puri <- score_prev
          se_puri <- se_prev
          dif_stat[item_num, 1:8] <- dif_stat_tmp[-flag_del, ]
          dif_stat[item_num, 9] <- i - 1
          contingency[item_num] <- contingency_tmp[-flag_del]

          # add the items flagged at the previous iteration, except the item deleted at this one
          dif_item <- c(dif_item, item_num[match(setdiff(dif_item_tmp, flag_del),
                                                 seq_len(nrow(dif_stat_tmp))[-flag_del])])
          break
        }

        # record the lowest reliabilities of the ability estimates in an iteration whose results are kept
        rho_min <- pmin(rho_min, dif_rst_tmp$rho_raw, na.rm = TRUE)

        # extract the DIF analysis results
        # and check if at least one DIF item is detected
        dif_item_tmp <- dif_rst_tmp$dif_item
        log_p_tmp <- dif_rst_tmp$log_p
        dif_stat_tmp <- dif_rst_tmp$dif_stat
        contingency_tmp <- dif_rst_tmp$contingency

        # check if a further DIF item is flagged
        if (is.null(dif_item_tmp)) {
          # add the DIF statistics for rest of items
          dif_stat[item_num, 1:8] <- dif_stat_tmp
          dif_stat[item_num, 9] <- i
          contingency[item_num] <- contingency_tmp

          break
        }
      }

      # print a message
      if (verbose) {
        cat("", "\n")
      }

      # record the actual number of iterations; the statistics are not recomputed in the
      # iteration that stops the purification
      n_iter <- if (all_flagged || no_score || no_test) i - 1L else i

      if (all_flagged) {
        # if every item is flagged, then, return a warning message
        warning("All items were flagged during the purification, so the purification stopped ",
                "with no item left to analyze.", call. = FALSE)
        complete <- FALSE
      } else if (no_score) {
        # if a group has too few ability estimates left, then, return a warning message
        warning("The purification stopped because fewer than two examinees of a group had an ",
                "ability estimate within 'range' for the remaining items.", call. = FALSE)
        complete <- FALSE
      } else if (no_test) {
        # if no remaining item has a usable bin, then, return a warning message
        warning("The purification stopped because no bin contained enough examinees of both groups ",
                "for any remaining item after the regression correction.", call. = FALSE)
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
        dif_stat[item_num, 1:8] <- dif_stat_tmp
        dif_stat[item_num, 9] <- i
        contingency[item_num] <- contingency_tmp
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
      with_purify$dif_stat <- dif_stat
      with_purify$dif_item <- sort(dif_item)
      with_purify$n.iter <- n_iter
      with_purify$complete <- complete
      with_purify$contingency <- contingency
    } else {
      # in case when no DIF item is detected from the first DIF analysis results
      with_purify$dif_stat <- cbind(no_purify$dif_stat, n.iter = 0)
      with_purify$n.iter <- 0L
      with_purify$complete <- TRUE
      with_purify$contingency <- no_purify$contingency
    }
  }


  # warn once when the reliability of the scores of a group was below the lower bound in any analysis
  if (any(rho_min < 0.05, na.rm = TRUE)) {
    warning(
      "The estimated reliability (rho^2) of the ability estimates was below 0.05 ",
      "for one or more groups (lowest rho^2_ref = ", round(rho_min[1], 3),
      ", lowest rho^2_foc = ", round(rho_min[2], 3), ") and was set to 0.05. ",
      "This can occur when some examinees have very large standard errors or when few items ",
      "are used to estimate the abilities. Interpret the DIF results with caution.",
      call. = FALSE
    )
  }

  # find the items whose statistics are missing in the reported tables, except the skipped items
  na_item <- which(is.na(no_purify$dif_stat$beta))
  if (purify) na_item <- union(na_item, which(is.na(with_purify$dif_stat$beta)))
  no_bin_id <- no_purify$dif_stat$id[setdiff(sort(na_item), item.skip)]

  # warn that the statistics of the items without a usable bin are missing
  if (length(no_bin_id) > 0L) {
    warning("No bin contained enough examinees of both groups after the regression correction for ",
            "item(s) ", paste(no_bin_id, collapse = ", "), ", so their statistics are missing and ",
            "they are not flagged.", call. = FALSE)
  }

  # summarize the results
  rst <- list(no_purify = no_purify, purify = purify,
              with_purify = with_purify, alpha = alpha)

  # return the DIF detection results
  class(rst) <- "catsib"
  rst$call <- cl
  rst
}


# This function checks whether each of the two groups has at least two ability estimates within
# the range used to compute the regression correction
catsib_enough <- function(score, group, focal.name, range) {
  # find the ability estimates strictly within the range
  in_range <- !is.na(score) & score > range[1] & score < range[2]

  # count the estimates of the focal group and the reference group
  sum(in_range[which(group == focal.name)]) >= 2L && sum(in_range[which(group != focal.name)]) >= 2L
}


# This function leaves out of the computation of rho^2 the scores whose squared standard error is
# equal to or greater than the variance of the scores of a group; the screening is not applied
# when fewer than two scores would be left
rho_screen <- function(score, se, use) {
  # find the scores whose squared standard error is not smaller than the variance of the scores
  out <- which(use & se^2 >= stats::var(score[use], na.rm = TRUE))

  # leave out those scores when at least two scores remain
  if (length(out) > 0L && sum(use, na.rm = TRUE) - length(out) >= 2L) {
    use[out] <- FALSE
  }
  use
}


# This function performs a regression correction for ability estimates and, then
# computes the beta statistic and its SE for all items
catsib_one <- function(data,
                       group,
                       focal.name,
                       item.skip = NULL,
                       score,
                       se,
                       range,
                       item.id,
                       max.bin,
                       min.bin,
                       min.binsize = 3,
                       max.del = 0.075,
                       weight.group,
                       alpha = 0.05,
                       rho.screen = TRUE) {

  ## ------------------------------------------
  ## prepare data sets
  ## ------------------------------------------
  # count the number of items
  nitem <- length(item.id)

  # find the location of examinees for the reference and the focal groups
  loc_ref <- which(group != focal.name)
  loc_foc <- which(group == focal.name)

  # divide the response data into the two group data
  resp_ref <- data[loc_ref, , drop = FALSE]
  resp_foc <- data[loc_foc, , drop = FALSE]

  # divide the thetas into the two group data
  score_ref <- score[loc_ref]
  score_foc <- score[loc_foc]

  # divide the SEs into the two group data
  se_ref <- se[loc_ref]
  se_foc <- se[loc_foc]

  ## ------------------------------------------
  ## a regression correction for theta scores
  ## ------------------------------------------
  # find the scores used to compute the mean, variance, and error variance of each group;
  # note that, the bound scores are excluded.
  lb_score <- range[1]
  up_score <- range[2]
  use_ref <- score_ref > lb_score & score_ref < up_score
  use_foc <- score_foc > lb_score & score_foc < up_score

  # when rho.screen = TRUE, leave out the scores whose squared standard error is equal to or
  # greater than the variance of the scores of the group
  if (rho.screen) {
    use_ref <- rho_screen(score = score_ref, se = se_ref, use = use_ref)
    use_foc <- rho_screen(score = score_foc, se = se_foc, use = use_foc)
  }

  # compute the mean and variance of scores for each group
  mu_ref <- mean(score_ref[use_ref], na.rm = TRUE)
  mu_foc <- mean(score_foc[use_foc], na.rm = TRUE)
  sigma2_ref <- stats::var(score_ref[use_ref], na.rm = TRUE)
  sigma2_foc <- stats::var(score_foc[use_foc], na.rm = TRUE)

  # compute the error variance of scores for each group
  errvar_ref <- mean((se_ref^2)[use_ref], na.rm = TRUE)
  errvar_foc <- mean((se_foc^2)[use_foc], na.rm = TRUE)

  # compute the raw (unclamped) squared correlation (reliability) between theta estimate and true theta
  rho_ref2_raw <- suppressWarnings(1 - errvar_ref / sigma2_ref)
  rho_foc2_raw <- suppressWarnings(1 - errvar_foc / sigma2_foc)

  # bound rho2 to [0.05, 1]; a negative value would reverse the order of the scores within a group,
  # and a value near zero would place almost all corrected scores of a group in one bin
  rho_ref2 <- max(0.05, min(1, rho_ref2_raw))
  rho_foc2 <- max(0.05, min(1, rho_foc2_raw))

  # apply a regression correction to the ability estimates
  crscore_ref <- mu_ref + rho_ref2 * (score_ref - mu_ref)
  crscore_foc <- mu_foc + rho_foc2 * (score_foc - mu_foc)

  ## ------------------------------------------
  ## compute CATSIB statistic
  ## ------------------------------------------
  # conduct a DIF analysis
  catsib_dif <- purrr::map(
    .x = 1:nitem,
    .f = function(i) {
      catsib_item(
        crscore_ref = crscore_ref, crscore_foc = crscore_foc,
        resp.ref = resp_ref[, i], resp.foc = resp_foc[, i],
        max.bin = max.bin, min.bin = min.bin, min.binsize = min.binsize,
        max.del = max.del, weight.group = weight.group
      )
    }
  )

  # extract the DIF analysis results
  dif_stat <-
    purrr::map(.x = catsib_dif, "dif_stat") %>%
    do.call(what = "rbind")
  contingency <- purrr::map(.x = catsib_dif, "contingency")
  names(contingency) <- item.id

  # find the items for which no bin contains enough examinees of both groups
  no_bin <- vapply(catsib_dif, function(d) isTRUE(d$no_bin), logical(1))
  no_bin[item.skip] <- FALSE

  # keep the unrounded p-values to flag the items and the log p-values to choose the item to be
  # removed in the purification
  p_beta <- dif_stat$p
  log_p <- log(2) + stats::pnorm(q = -abs(dif_stat$z.beta), mean = 0, sd = 1, log.p = TRUE)

  # round the statistics
  dif_stat <-
    dif_stat %>%
    dplyr::mutate_if("is.numeric", "round", digits = 4) %>%
    dplyr::mutate(id = item.id) %>%
    dplyr::relocate("id", .before = "beta")

  # when there are items that should be skipped for the DIF analysis
  # insert NAs to the corresponding results of the items
  if (!is.null(item.skip)) {
    dif_stat[item.skip, 2:5] <- NA
    p_beta[item.skip] <- NA
    log_p[item.skip] <- NA
    contingency[item.skip] <- NA
  }

  # find the flagged items using the unrounded p-values
  dif_item <- as.numeric(which(p_beta <= alpha))
  if (length(dif_item) == 0) dif_item <- NULL

  # return the results
  rst <- list(
    dif_stat = dif_stat, dif_item = dif_item, contingency = contingency, log_p = log_p,
    no_bin = item.id[no_bin], rho_raw = c(ref = rho_ref2_raw, foc = rho_foc2_raw)
  )
  rst
}


# This function computes the beta statistic and its SE for an item
catsib_item <- function(crscore_ref, crscore_foc, resp.ref, resp.foc,
                        max.bin, min.bin, min.binsize = 3, max.del = 0.075,
                        weight.group) {

  # check if all responses are NAs
  all.na.ref <- all(is.na(resp.ref))
  all.na.foc <- all(is.na(resp.foc))
  
  # when all responses are NAs for either group, return NAs
  if(all.na.ref | all.na.foc) {
    list(dif_stat = data.frame(
      beta = NA, se = NA, z.beta = NA, p = NA,
      n.ref = 0, n.foc = 0, n.total = 0
    ),
    contingency = data.frame(
      bin = NA, n.ref = 0, prop.ref = NA, var.ref = NA,
      n.foc = 0, prop.foc = NA, var.foc = NA,
      n.total = 0, beta = NA, var.beta = NA
    ))
    
  } else {
    
    # when not all responses are NAs, proceed to compute the beta statistic
  
    # combine all corrected theta scores
    crscore <- c(crscore_ref, crscore_foc)
    
    # set the range of the ability scale
    min.crscore <- min(crscore, na.rm = TRUE)
    max.crscore <- max(crscore, na.rm = TRUE)
    
    # decide the number of bins and create an initial frequency table
    for (num.bin in max.bin:min.bin) {
      # compute the cut-scores to divide the theta scale into the bins
      cutscore <- seq(from = min.crscore, to = max.crscore, length.out = num.bin + 1)
      
      # assign a group variable to each score
      bin_ref <- cut(crscore_ref, breaks = cutscore, include.lowest = TRUE, dig.lab = 7)
      bin_foc <- cut(crscore_foc, breaks = cutscore, include.lowest = TRUE, dig.lab = 7)
      
      # keep the bins of the examinees who responded to the item
      non.na.ref <- !is.na(resp.ref)
      non.na.foc <- !is.na(resp.foc)
      bin_ref <- bin_ref[non.na.ref]
      bin_foc <- bin_foc[non.na.foc]
      
      # create a temporary data frame of bin frequency for both groups
      bin.n.ref <- stats::xtabs(~bin_ref, drop.unused.levels = FALSE)
      bin.n.foc <- stats::xtabs(~bin_foc, drop.unused.levels = FALSE)
      
      # check if the proportion of retained examinees is at least 1 - max.del in each group
      isok_ref <- (sum(bin.n.ref[bin.n.ref >= min.binsize & bin.n.foc >= min.binsize]) /
                     sum(bin.n.ref)) >= 1 - max.del
      isok_foc <- (sum(bin.n.foc[bin.n.ref >= min.binsize & bin.n.foc >= min.binsize]) /
                     sum(bin.n.foc)) >= 1 - max.del
      
      # if the criterion is met, then break out the loop
      if (all(isok_ref, isok_foc)) {
        break
      }
    }
    
    # final data frame containing all components to compute the beta statistic
    prop.ref <- 
      as.numeric(table(bin_ref, 
                       factor(resp.ref[non.na.ref], levels = c(0, 1)))[, 2] / bin.n.ref)
    prop.foc <- 
      as.numeric(table(bin_foc, 
                       factor(resp.foc[non.na.foc], levels = c(0, 1)))[, 2] / bin.n.foc)
    var.ref <- as.numeric(by(
      data = resp.ref[non.na.ref],
      INDICES = bin_ref,
      FUN = stats::var, na.rm = TRUE
    ))
    var.foc <- as.numeric(by(
      data = resp.foc[non.na.foc],
      INDICES = bin_foc,
      FUN = stats::var, na.rm = TRUE
    ))
    ref.df <-
      as.data.frame(bin.n.ref, stringsAsFactors = FALSE) %>%
      stats::setNames(nm = c("bin", "n.ref")) %>%
      data.frame(prop.ref = prop.ref, var.ref = var.ref)
    foc.df <-
      as.data.frame(bin.n.foc, stringsAsFactors = FALSE) %>%
      stats::setNames(nm = c("bin", "n.foc")) %>%
      data.frame(prop.foc = prop.foc, var.foc = var.foc)
    n.ref <- n.foc <- weight <- NULL
    item_df <-
      merge(x = ref.df, y = foc.df, by = "bin", all = TRUE, sort = FALSE) %>%
      # keep only bins with at least 'min.binsize' examinees in BOTH groups,
      # consistent with the bin-count selection loop above
      subset(n.ref >= min.binsize & n.foc >= min.binsize) %>%
      transform(n.total = n.ref + n.foc) %>%
      dplyr::mutate(
        # weight each bin by the group whose ability distribution is the target
        weight = switch(weight.group,
          comb = .data$n.total / sum(.data$n.total),
          foc = .data$n.foc / sum(.data$n.foc),
          ref = .data$n.ref / sum(.data$n.ref)
        )
      ) %>%
      transform(
        beta = ((prop.ref - prop.foc) * weight),
        var.beta = ((var.ref / n.ref + var.foc / n.foc) * weight^2)
      )
    
    # when no bin contains enough examinees of both groups, return NAs and mark the item
    if (nrow(item_df) == 0L) {
      return(list(
        dif_stat = data.frame(
          beta = NA, se = NA, z.beta = NA, p = NA,
          n.ref = 0, n.foc = 0, n.total = 0
        ),
        contingency = data.frame(
          bin = NA, n.ref = 0, prop.ref = NA, var.ref = NA,
          n.foc = 0, prop.foc = NA, var.foc = NA,
          n.total = 0, beta = NA, var.beta = NA
        ),
        no_bin = TRUE
      ))
    }

    # compute the beta statistic and its SE
    beta <- sum(item_df$beta)
    se_beta <- sqrt(sum(item_df$var.beta))
    z_beta <- beta / se_beta
    pval_beta <- 2 * stats::pnorm(q = abs(z_beta), mean = 0, sd = 1, lower.tail = FALSE)
    n.ref <- sum(item_df$n.ref)
    n.foc <- sum(item_df$n.foc)
    stat_df <- data.frame(
      beta = beta, se = se_beta, z.beta = z_beta, p = pval_beta,
      n.ref = n.ref, n.foc = n.foc, n.total = n.ref + n.foc
    )
    
    # round the numbers of the data frame
    item_df2 <-
      janitor::adorn_totals(
        dat = item_df, where = "row", fill = NA, , , "n.ref", "n.foc", "n.total",
        "beta", "var.beta"
      ) %>%
      dplyr::mutate_at(.vars = c(3, 4, 6, 7), "round", digits = 4)
    
    # return the results
    list(dif_stat = stat_df, contingency = item_df2)
  
  }
  
}
