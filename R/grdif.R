#' Generalized IRT residual-based DIF detection framework for multiple groups
#' (GRDIF)
#'
#' This function computes three GRDIF statistics, \eqn{GRDIF_{R}},
#' \eqn{GRDIF_{S}}, and \eqn{GRDIF_{RS}}, for analyzing differential item
#' functioning (DIF) among multiple groups (Lim et al., 2024). They are
#' designed mainly to detect uniform DIF, nonuniform DIF, and mixed DIF,
#' respectively.
#'
#' @inheritParams rdif
#' @param x A data frame containing item metadata (e.g., item parameters,
#'   number of categories, IRT model types, etc.), or an object of class
#'   `est_irt` or `est_item` obtained from [irtQ::est_irt()] or
#'   [irtQ::est_item()]. All items must be dichotomous (`"1PLM"`, `"2PLM"`,
#'   `"3PLM"`, or `"DRM"`); the function stops when any item is polytomous.
#'
#'   See [irtQ::est_irt()] or [irtQ::simdat()] for more details about the item
#'   metadata. This data frame can be easily created using the
#'   [irtQ::shape_df()] function.
#' @param data A matrix of examinees' scored responses (0 or 1) to the
#'   dichotomous items specified in the `x` argument. Rows represent examinees
#'   and columns represent items.
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
#' @param group A numeric or character vector of examinees' group membership,
#'   with one element per row of `data`. It should contain the value of the
#'   reference group and the values in `focal.name`. Examinees whose value is not
#'   in `focal.name` form the reference group, and the function stops when they
#'   have more than one value. Examinees with `NA` in `group` are excluded from
#'   the analysis.
#' @param focal.name A numeric or character vector of the values in `group` that
#'   identify the focal groups, with at least one value. For example, given
#'   `group = c(0, 0, 1, 2, 2, 3, 3)`, where '0' indicates the reference group
#'   and '1', '2', and '3' indicate three focal groups, set `focal.name = c(1,
#'   2, 3)`. The order of `focal.name` determines the order of the focal groups
#'   in the results (for example, the columns `n.foc1`, `n.foc2`, and `n.foc3`);
#'   the GRDIF statistics do not depend on it. All examinees whose group value
#'   is not in `focal.name` form the reference group (see **Details**). Each
#'   value must be found in `group`, and the values must be distinct.
#' @param alpha A numeric value specifying the significance level (\eqn{\alpha})
#'   of the tests. An item is flagged by a statistic when its p-value is less
#'   than or equal to `alpha`. The comparison uses unrounded p-values, although
#'   the p-values reported in `dif_stat` are rounded to four decimal places. The
#'   same level is used in the post-hoc pairwise analyses, without an adjustment
#'   for multiple comparisons. Default is `0.05`.
#' @param purify.by A character string specifying the GRDIF statistic used for
#'   purification: "grdifrs" for \eqn{GRDIF_{RS}}, "grdifr" for
#'   \eqn{GRDIF_{R}}, or "grdifs" for \eqn{GRDIF_{S}}. Used only when `purify =
#'   TRUE`. Default is `"grdifrs"`.
#' @param post.hoc Logical. If `TRUE`, post-hoc pairwise RDIF analyses are
#'   performed, for every pair of groups, on the items flagged by each GRDIF
#'   statistic. The p-values of the post-hoc analyses are not adjusted for
#'   multiple comparisons. See **Details**. Default is `TRUE`.
#' @param range A numeric vector of length two specifying the lower and upper
#'   bounds of the ability scale used when abilities are estimated internally
#'   with the scoring methods `"ML"`, `"MLF"`, `"WL"`, and `"MAP"`. Default is
#'   `c(-5, 5)`, the default of the other DIF functions.
#'
#' @details The GRDIF framework (Lim et al., 2024) generalizes the RDIF
#'   framework (Lim & Choe, 2023; Lim et al., 2022) to more than two groups. For
#'   each item and group, the mean raw residual (MRR) and the mean squared
#'   residual (MSR) are computed, where a residual is the observed item score
#'   minus the model-expected item score evaluated at the examinee's ability
#'   estimate. Under the null hypothesis that the item has no DIF between any
#'   pair of groups, the MRRs and MSRs asymptotically follow normal
#'   distributions whose means, variances, and covariances are computed
#'   analytically from the model-predicted probabilities (Lim et al., 2024,
#'   Equations 2 to 6). \eqn{GRDIF_{R}}, \eqn{GRDIF_{S}}, and \eqn{GRDIF_{RS}}
#'   are quadratic forms of the contrasts between the focal groups and the
#'   reference group of the MRRs, of the MSRs, and of both, centered at their
#'   null means, in the inverse of their null covariance matrices (Lim et al.,
#'   2024, Equation 10). Under the null hypothesis, they asymptotically follow
#'   \eqn{\chi^{2}} distributions with G-1, G-1, and 2(G-1) degrees of freedom,
#'   respectively, where G is the number of groups. The statistics do not depend
#'   on which group is the reference group. For two groups, \eqn{GRDIF_{R}} and
#'   \eqn{GRDIF_{S}} equal the squares of the standardized \eqn{RDIF_{R}} and
#'   \eqn{RDIF_{S}} statistics of [irtQ::rdif()], \eqn{GRDIF_{RS}} equals
#'   \eqn{RDIF_{RS}}, and the p-values are the same.
#'
#'   [irtQ::grdif()] supports dichotomous items only, because the GRDIF
#'   framework has been developed and evaluated for dichotomous items (Lim et
#'   al., 2024). The function stops with an error when `x` includes any item
#'   with three or more score categories or any item of the "GRM" or "GPCM"
#'   model.
#'
#'   To compute the GRDIF statistics, [irtQ::grdif()] requires: (1) item
#'   parameter estimates obtained from aggregate data (regardless of group
#'   membership), (2) examinees' ability estimates (e.g., ML), and (3) examinees'
#'   item response data. Note that the ability estimates must be based on the
#'   aggregate-data item parameters. The item parameter estimates should be
#'   provided in the `x` argument, the ability estimates in the `score`
#'   argument, and the response data in the `data` argument. If ability
#'   estimates are not provided (i.e., `score = NULL`), [irtQ::grdif()] estimates
#'   them with the scoring method specified in the `method` argument (e.g.,
#'   `method = "ML"`). When `x` is an object of class `est_irt`, the item
#'   parameter estimates, the response data, and the scaling constant `D` are
#'   taken from the object. When `x` is an object of class `est_item`, the
#'   ability values stored in the object are also used as `score`.
#'
#'   In rare cases, the covariance matrix of the contrasts that define
#'   \eqn{GRDIF_{R}}, \eqn{GRDIF_{S}}, or \eqn{GRDIF_{RS}} is singular or nearly
#'   singular (its reciprocal condition number is 1e-10 or less), for example
#'   when all examinees of a group who responded to an item have the same
#'   probability of a correct response, so that the MSR is a linear function of
#'   the MRR. Then the statistic is computed with the Moore-Penrose generalized
#'   inverse of the covariance matrix and compared with a chi-square
#'   distribution whose degrees of freedom equal the rank of the matrix (Moore,
#'   1977), and one warning names the items. This handling is not part of the
#'   original method (Lim et al., 2024) and is added in irtQ. In most data the
#'   matrix is not singular, and the degrees of freedom are G-1, G-1, and
#'   2(G-1).
#'
#'   The `group` argument contains the group membership of each examinee as
#'   numeric or character values, and its length must equal the number of rows
#'   in `data`. The values listed in `focal.name` identify the focal groups, and
#'   all other examinees form the reference group. The `group` vector should
#'   therefore contain one value for the reference group in addition to the
#'   values in `focal.name`, and the function stops when it contains more than
#'   one such value. Examinees with `NA` in `group` are excluded from the
#'   analysis.
#'
#'   As in the RDIF framework, an iterative purification procedure can be applied
#'   (Lim et al., 2022, 2024). When `purify = TRUE`, the statistic specified in
#'   `purify.by` drives the procedure (e.g., `purify.by = "grdifrs"`). At each
#'   iteration, the flagged item with the smallest p-value of the `purify.by`
#'   statistic is removed (the p-values are compared on the log scale), the
#'   abilities are re-estimated from the remaining items with the scoring method
#'   specified in `method`, and the GRDIF statistics of the remaining items are
#'   recomputed. A supplied `score` is therefore used only in the initial
#'   analysis. The procedure stops when no remaining item is flagged or when
#'   `max.iter` iterations have been performed. If no item is flagged in the
#'   initial analysis, no iteration is performed. If every item is flagged, the
#'   procedure stops with a warning, because no item is left to estimate the
#'   abilities, and `complete` is `FALSE`. The procedure also stops with a
#'   warning and `complete = FALSE` when, after an item is removed, no examinee
#'   of the reference group or of a focal group has an ability estimate left
#'   (see `min.resp` below); the remaining items then keep the results and the
#'   scores of the previous iteration, and the items flagged in that iteration
#'   are added to the flagged items.
#'
#'   Scoring based on a small number of item responses can lead to large
#'   standard errors, potentially reducing the accuracy of DIF detection in the
#'   GRDIF framework. The `min.resp` argument excludes such examinees. For
#'   example, if `min.resp = 5`, examinees who responded to fewer than five
#'   items (but at least one) have all their responses treated as missing (i.e.,
#'   `NA`) before the abilities are estimated, in the initial analysis and at
#'   every purification iteration, and a warning reports their number in the
#'   initial analysis. As a result, their ability estimates are missing, and they
#'   are not used in the computation of the GRDIF statistics. When `score` is
#'   supplied, the scores of these examinees are also treated as missing. During
#'   purification, the count is based on the items that remain after the flagged
#'   items are removed, and an examinee excluded once stays excluded. If
#'   `min.resp = NULL`, a score is computed for any examinee with at least one
#'   valid item response.
#'
#'   Examinees whose ability estimate is `NA`, such as examinees without any item
#'   response or examinees with `NA` in a supplied `score`, are excluded from the
#'   computation of the GRDIF statistics and from the post-hoc analyses. A
#'   warning reports the number of examinees with item responses who are
#'   excluded because of a missing value in `score`; these examinees stay
#'   excluded during purification. During purification, an examinee whose
#'   remaining responses are all removed with the flagged items, or who has fewer
#'   than `min.resp` responses to the remaining items, has no ability estimate
#'   and is excluded, and one warning at the end reports the number of such
#'   examinees. The sample sizes `n.ref`, `n.foc1`, `n.foc2`, and so on count
#'   only the examinees who are used.
#'
#'   If no examinee of one of the groups responded to an item, the GRDIF
#'   statistics and p-values of the item are `NaN`, and the item is not flagged.
#'
#'   The GRDIF statistics test only whether DIF is present between at least one
#'   pair of groups (Lim et al., 2024). When `post.hoc = TRUE`, the three RDIF
#'   statistics of [irtQ::rdif()] are computed for every pair of distinct values
#'   in `group`, for each item flagged by \eqn{GRDIF_{R}}, \eqn{GRDIF_{S}}, or
#'   \eqn{GRDIF_{RS}}. For example, with four groups A, B, C, and D, the pairs
#'   A-B, A-C, A-D, B-C, B-D, and C-D are compared. Each pairwise analysis uses
#'   the responses of the examinees in the two groups, their ability estimates
#'   from the GRDIF analysis, and the item parameters in `x`; the abilities are
#'   not re-estimated for the pair. In a pair labeled "a & b" in the column
#'   `group.pair`, a is the reference group and b is the focal group of the
#'   pair, so \eqn{RDIF_{R}} and \eqn{RDIF_{S}} are the residual means of group
#'   b minus those of group a. The reference group of the GRDIF analysis is the
#'   reference group of every pair that includes it (Lim et al., 2024, Equation
#'   14), and for a pair of focal groups, the group listed earlier in
#'   `focal.name` is the reference group. The p-values are not adjusted for
#'   multiple comparisons; to control the familywise error rate, compare them
#'   with a corrected level, for example `alpha` divided by the number of pairs
#'   (a Bonferroni correction). When `purify = TRUE`, the post-hoc results in
#'   `with_purify` are computed for each item when it is removed during
#'   purification, with the ability estimates of that iteration; the items added
#'   to the flagged items when the purification ends before it is complete
#'   (`max.iter` reached or no ability estimate left in a group) are not
#'   analyzed.
#'
#'   The print method of this class has an argument `digits` that is not used;
#'   the statistics are printed rounded to three decimal places.
#'
#' @return This function returns an object of class `"grdif"`, which is a list
#' with the following five components:
#'
#' \item{no_purify}{A list of sub-objects presenting the results of DIF analysis
#' without a purification procedure. These include:
#'   \describe{
#'     \item{dif_stat}{A data frame with one row per item and the columns `id`,
#'     `grdifr`, `grdifs`, `grdifrs`, `p.grdifr`, `p.grdifs`, `p.grdifrs`,
#'     `n.ref`, `n.foc1`, `n.foc2`, ..., and `n.total`: the item ID, the three
#'     GRDIF statistics, their p-values, the numbers of examinees in the
#'     reference group and in each focal group (in the order of `focal.name`)
#'     who responded to the item, and their sum. Statistics and p-values are
#'     rounded to four decimal places.}
#'     \item{moments}{A list of three data frames, `mu`, `var`, and
#'     `covariance`, with the null moments of the MRR and MSR of each group,
#'     computed analytically under the null hypothesis of no DIF (Lim et al.,
#'     2024, Equations 3 to 5). `mu` and `var` have the columns `id`, `mrr.ref`,
#'     `msr.ref`, `mrr.foc1`, `msr.foc1`, and so on, with the null means and null
#'     variances of the MRR and MSR of the reference group and of each focal
#'     group. `covariance` has the columns `id`, `ref`, `foc1`, and so on, with
#'     the null covariance between the MRR and MSR of each group.}
#'     \item{dif_item}{A list with the elements `grdifr`, `grdifs`, and
#'     `grdifrs`, each giving the positions (rows of `x`) of the items flagged by
#'     the corresponding statistic, or `NULL` if no item is flagged.}
#'     \item{score}{A numeric vector of the ability estimates used in the initial
#'     analysis: the supplied `score`, or the internal estimates when
#'     `score = NULL`.}
#'     \item{post.hoc}{A list with the elements `by.grdifr`, `by.grdifs`, and
#'     `by.grdifrs`. Each element is a data frame with the post-hoc results for
#'     the items flagged by the corresponding GRDIF statistic, with one row per
#'     item and pair of groups and the columns `id`, `group.pair`, `rdifr`,
#'     `z.rdifr`, `rdifs`, `z.rdifs`, `rdifrs`, `p.rdifr`, `p.rdifs`, `p.rdifrs`,
#'     `n.ref`, `n.foc`, and `n.total`, defined as in the `dif_stat` element of
#'     [irtQ::rdif()]; `n.ref` and `n.foc` refer to the first and second groups
#'     in `group.pair`. An element is `NULL` if no item is flagged by the
#'     statistic or if `post.hoc = FALSE`.}
#'   }
#' }
#'
#' \item{purify}{A logical value indicating whether purification was requested
#' (the `purify` argument).}
#'
#' \item{with_purify}{A list with the results of the purification procedure. All
#' elements are `NULL` when `purify = FALSE`. The elements are:
#'   \describe{
#'     \item{purify.by}{A character string indicating which GRDIF statistic was
#'     used for purification: "grdifr", "grdifs", or "grdifrs", corresponding to
#'     \eqn{GRDIF_{R}}, \eqn{GRDIF_{S}}, and \eqn{GRDIF_{RS}}, respectively.}
#'     \item{dif_stat}{A data frame with the same columns as
#'     `no_purify$dif_stat` and an additional column `n.iter`. For an item
#'     removed during purification, the row reports the statistics from the
#'     iteration in which the item was removed, and `n.iter` is that iteration
#'     minus 1 (0 for an item removed on the basis of the initial analysis). For
#'     the other items, the row reports the statistics from the last iteration,
#'     and `n.iter` is the number of that iteration. If no item is flagged in the
#'     initial analysis, this is `no_purify$dif_stat` with `n.iter = 0`.}
#'     \item{moments}{A list of three data frames with the same columns as in
#'     `no_purify$moments` and an additional column `n.iter`, defined in the
#'     same way as in `dif_stat`.}
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
#'     iteration counted in `n.iter`. `NULL` if no item is flagged in the
#'     initial analysis, because purification is not carried out; the initial
#'     estimates are then in `no_purify$score`.}
#'     \item{post.hoc}{A data frame with the same columns as the elements of
#'     `no_purify$post.hoc` and an additional column `n_iter`. It contains the
#'     post-hoc results of each item removed during purification, computed with
#'     the ability estimates of the iteration in which the item was removed;
#'     `n_iter` is defined as `n.iter` in `dif_stat`. `NULL` if
#'     `post.hoc = FALSE` or if no item is flagged in the initial analysis.}
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
#'   [irtQ::simdat()], [irtQ::shape_df()], [irtQ::est_score()]
#'
#' @references Bock, R. D., & Mislevy, R. J. (1982). Adaptive EAP estimation of
#'   ability in a microcomputer environment. *Applied Psychological Measurement,
#'   6*(4), 431-444. \doi{10.1177/014662168200600405}.
#'
#'   Hambleton, R. K., Swaminathan, H., & Rogers, H. J. (1991). *Fundamentals of
#'   item response theory*. Newbury Park, CA: Sage.
#'
#'   Lim, H., & Choe, E. M. (2023). Detecting differential item functioning in
#'   CAT using IRT residual DIF approach. *Journal of Educational Measurement,
#'   60*(4), 626-650. \doi{10.1111/jedm.12366}.
#'
#'   Lim, H., Choe, E. M., & Han, K. T. (2022). A residual-based differential
#'   item functioning detection framework in item response theory. *Journal of
#'   Educational Measurement, 59*(1), 80-104. \doi{10.1111/jedm.12313}.
#'
#'   Lim, H., Zhu, D., Choe, E. M., & Han, K. T. (2024). Detecting differential
#'   item functioning among multiple groups using IRT residual DIF framework.
#'   *Journal of Educational Measurement, 61*(4), 656-681.
#'   \doi{10.1111/jedm.12415}.
#'
#'   Moore, D. S. (1977). Generalized inverses, Wald's method, and the
#'   construction of chi-squared tests of fit. *Journal of the American
#'   Statistical Association, 72*(357), 131-137.
#'   \doi{10.1080/01621459.1977.10479921}.
#'
#'   Warm, T. A. (1989). Weighted likelihood estimation of ability in item
#'   response theory. *Psychometrika, 54*(3), 427-450.
#'   \doi{10.1007/BF02294627}.
#'
#' @examples
#' \donttest{
#' # Load required library
#' library("dplyr")
#'
#' ## Uniform DIF detection for four groups (1 reference, 3 focal)
#' ########################################################
#' # (1) Manipulate uniform DIF for all three focal groups
#' ########################################################
#'
#' # Import the "-prm.txt" output file from flexMIRT
#' flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
#'
#' # Select 36 non-DIF items modeled under 3PLM
#' par_nstd <-
#'   bring.flexmirt(file = flex_sam, "par")$Group1$full_df %>%
#'   dplyr::filter(.data$model == "3PLM") %>%
#'   dplyr::filter(dplyr::row_number() %in% 1:36) %>%
#'   dplyr::select(1:6)
#' par_nstd$id <- paste0("nondif", 1:36)
#'
#' # Generate four new items on which uniform DIF will be imposed
#' difpar_ref <-
#'   shape_df(
#'     par.drm = list(a = c(0.8, 1.5, 0.8, 1.5), b = c(0.0, 0.0, -0.5, -0.5), g = rep(0.15, 4)),
#'     item.id = paste0("dif", 1:4), cats = 2, model = "3PLM"
#'   )
#'
#' # Introduce DIF by shifting b-parameters differently for each focal group
#' difpar_foc1 <-
#'   difpar_ref %>%
#'   dplyr::mutate_at(.vars = "par.2", .funs = function(x) x + c(0.7, 0.7, 0, 0))
#' difpar_foc2 <-
#'   difpar_ref %>%
#'   dplyr::mutate_at(.vars = "par.2", .funs = function(x) x + c(0, 0, 0.7, 0.7))
#' difpar_foc3 <-
#'   difpar_ref %>%
#'   dplyr::mutate_at(.vars = "par.2", .funs = function(x) x + c(-0.4, -0.4, -0.5, -0.5))
#'
#' # Combine the 4 DIF and 36 non-DIF items for all four groups
#' # Therefore, each of the first four items has uniform DIF in two of the three
#' # focal groups
#' par_ref <- rbind(difpar_ref, par_nstd)
#' par_foc1 <- rbind(difpar_foc1, par_nstd)
#' par_foc2 <- rbind(difpar_foc2, par_nstd)
#' par_foc3 <- rbind(difpar_foc3, par_nstd)
#'
#' # Generate true abilities from different distributions
#' set.seed(128)
#' theta_ref <- rnorm(500, 0.0, 1.0)
#' theta_foc1 <- rnorm(500, -1.0, 1.0)
#' theta_foc2 <- rnorm(500, 1.0, 1.0)
#' theta_foc3 <- rnorm(500, 0.5, 1.0)
#'
#' # Simulate response data for each group
#' resp_ref <- irtQ::simdat(par_ref, theta = theta_ref, D = 1)
#' resp_foc1 <- irtQ::simdat(par_foc1, theta = theta_foc1, D = 1)
#' resp_foc2 <- irtQ::simdat(par_foc2, theta = theta_foc2, D = 1)
#' resp_foc3 <- irtQ::simdat(par_foc3, theta = theta_foc3, D = 1)
#' data <- rbind(resp_ref, resp_foc1, resp_foc2, resp_foc3)
#'
#' ########################################################
#' # (2) Estimate item and ability parameters
#' #     using aggregated data
#' ########################################################
#'
#' # Estimate item parameters
#' est_mod <- irtQ::est_irt(data = data, D = 1, model = "3PLM")
#' est_par <- est_mod$par.est
#'
#' # Estimate ability parameters using MLE
#' score <- irtQ::est_score(x = est_par, data = data, method = "ML")$est.theta
#'
#' ########################################################
#' # (3) Conduct DIF analysis
#' ########################################################
#'
#' # Create a group membership vector:
#' # 0 = reference group; 1, 2, 3 = focal groups
#' group <- c(rep(0, 500), rep(1, 500), rep(2, 500), rep(3, 500))
#'
#' # (a) Compute GRDIF statistics without purification,
#' #     and perform post-hoc pairwise comparisons for flagged items
#' dif_nopuri <- grdif(
#'   x = est_par, data = data, score = score, group = group,
#'   focal.name = c(1, 2, 3), D = 1, alpha = 0.05,
#'   purify = FALSE, post.hoc = TRUE
#' )
#' print(dif_nopuri)
#'
#' # Display post-hoc pairwise comparison results
#' print(dif_nopuri$no_purify$post.hoc)
#'
#' # (b) Compute GRDIF statistics with purification based on GRDIF_R,
#' #     including post-hoc comparisons (not run: it takes about 10 seconds
#' #     because the abilities of the examinees are re-estimated at every
#' #     iteration)
#' #
#' #   dif_puri_r <- grdif(
#' #     x = est_par, data = data, score = score, group = group,
#' #     focal.name = c(1, 2, 3), D = 1, alpha = 0.05,
#' #     purify = TRUE, purify.by = "grdifr", post.hoc = TRUE
#' #   )
#' #   print(dif_puri_r)
#' #
#' #   # Display post-hoc results after purification
#' #   print(dif_puri_r$with_purify$post.hoc)
#' #
#' #   purify.by = "grdifs" or "grdifrs" can be used in the same way.
#' }
#'
#' @export
grdif <- function(x, ...) UseMethod("grdif")

#' @describeIn grdif Default method for a data frame `x` containing item metadata.
#'
#' @import dplyr
#' @export
grdif.default <- function(x,
                          data,
                          score = NULL,
                          group,
                          focal.name,
                          D = 1,
                          alpha = 0.05,
                          missing = NA,
                          purify = FALSE,
                          purify.by = c("grdifrs", "grdifr", "grdifs"),
                          max.iter = 10,
                          min.resp = NULL,
                          post.hoc = TRUE,
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
  rst <- grdif_main(
    x = x, data = data, score = score, group = group, focal.name = focal.name, D = D,
    alpha = alpha, missing = missing, purify = purify, purify.by = purify.by,
    max.iter = max.iter, min.resp = min.resp, post.hoc = post.hoc, method = method,
    range = range, norm.prior = norm.prior, nquad = nquad, weights = weights,
    ncore = ncore, verbose = verbose, ...
  )

  # return the DIF detection results
  rst$call <- cl
  rst
}

#' @describeIn grdif Method for an object of class `est_irt` created by [irtQ::est_irt()].
#' The response data and the scaling constant `D` are taken from the object.
#' @import dplyr
#' @export
#'
grdif.est_irt <- function(x, score = NULL, group, focal.name, alpha = 0.05, missing = NA, purify = FALSE,
                          purify.by = c("grdifrs", "grdifr", "grdifs"), max.iter = 10, min.resp = NULL, post.hoc = TRUE,
                          method = "ML", range = c(-5, 5), norm.prior = c(0, 1), nquad = 41, weights = NULL, ncore = 1,
                          verbose = TRUE, ...) {
  # match.call
  cl <- match.call()

  # stop when the response data or the scaling constant is passed, since they are taken from the object
  check_obj_dots(list(...))

  # conduct the DIF analysis with the data, scaling constant, and items of the object
  rst <- grdif_main(
    x = x$par.est, data = x$data, score = score, group = group, focal.name = focal.name,
    D = x$scale.D, alpha = alpha, missing = missing, purify = purify, purify.by = purify.by,
    max.iter = max.iter, min.resp = min.resp, post.hoc = post.hoc, method = method,
    range = range, norm.prior = norm.prior, nquad = nquad, weights = weights,
    ncore = ncore, verbose = verbose, ...
  )

  # return the DIF detection results
  rst$call <- cl
  rst
}


#' @describeIn grdif Method for an object of class `est_item` created by [irtQ::est_item()].
#' The response data, the ability estimates, and the scaling constant `D` are taken
#' from the object.
#' @import dplyr
#' @export
#'
grdif.est_item <- function(x, group, focal.name, alpha = 0.05, missing = NA, purify = FALSE,
                           purify.by = c("grdifrs", "grdifr", "grdifs"), max.iter = 10, min.resp = NULL, post.hoc = TRUE,
                           method = "ML", range = c(-5, 5), norm.prior = c(0, 1), nquad = 41, weights = NULL, ncore = 1,
                           verbose = TRUE, ...) {
  # match.call
  cl <- match.call()

  # stop when the response data, the ability estimates, or the scaling constant is passed, since
  # they are taken from the object
  check_obj_dots(list(...), obj_score = TRUE)

  # conduct the DIF analysis with the data, ability estimates, scaling constant, and items of the object
  rst <- grdif_main(
    x = x$par.est, data = x$data, score = x$score, group = group, focal.name = focal.name,
    D = x$scale.D, alpha = alpha, missing = missing, purify = purify, purify.by = purify.by,
    max.iter = max.iter, min.resp = min.resp, post.hoc = post.hoc, method = method,
    range = range, norm.prior = norm.prior, nquad = nquad, weights = weights,
    ncore = ncore, verbose = verbose, ...
  )

  # return the DIF detection results
  rst$call <- cl
  rst
}


# This function conducts the GRDIF analysis shared by the methods of grdif(): it prepares the
# data and ability estimates, runs the analysis once, applies the purification procedure, and
# conducts the post-hoc pairwise RDIF analyses
grdif_main <- function(x, data, score, group, focal.name, D, alpha, missing, purify, purify.by,
                       max.iter, min.resp, post.hoc, method, range, norm.prior, nquad, weights,
                       ncore, verbose, ...) {

  ## ----------------------------------
  ## (1) prepare DIF analysis
  ## ----------------------------------
  # confirm and correct all item metadata information
  x <- confirm_df(x)

  # stop when any item is polytomous because the GRDIF framework is for dichotomous items
  if (any(x$cats != 2L) || any(x$model %in% c("GRM", "GPCM"))) {
    stop("grdif() supports dichotomous items only. The GRDIF framework has been evaluated for ",
         "dichotomous items (Lim et al., 2024).", call. = FALSE)
  }

  # record whether the ability estimates are supplied
  score_given <- !is.null(score)

  # check the inputs and prepare the response data and ability estimates
  prep <- dif_prepare(
    x = x, data = data, score = score, group = group, focal.name = focal.name,
    item.skip = NULL, D = D, alpha = alpha, missing = missing, purify = purify,
    max.iter = max.iter, min.resp = min.resp, method = method, range = range,
    norm.prior = norm.prior, nquad = nquad, weights = weights, ncore = ncore,
    fn = "grdif()", multi_focal = TRUE, ...
  )
  data <- prep$data
  score <- prep$score

  # a) when no purification is set
  # do only one iteration of DIF analysis
  dif_rst <- grdif_one(x = x, data = data, score = score, group = group, focal.name = focal.name, D = D, alpha = alpha)

  # record the items whose covariance matrix is singular
  singular_id <- dif_rst$singular

  # record the items whose covariance matrix is singular in the post-hoc analyses
  singular_post <- NULL

  # create two empty lists to contain the results
  no_purify <- list(
    dif_stat = NULL, moments = NULL, dif_item = NULL, score = NULL,
    post.hoc = list(by.grdifr = NULL, by.grdifs = NULL, by.grdifrs = NULL)
  )
  with_purify <- list(
    purify.by = NULL, dif_stat = NULL, moments = NULL,
    dif_item = NULL, n.iter = NULL, score = NULL, post.hoc = NULL, complete = NULL
  )

  # record the first DIF detection results into the no purification list
  no_purify$dif_stat <- dif_rst$dif_stat
  no_purify$dif_item <- dif_rst$dif_item
  no_purify$moments <-
    list(
      mu = data.frame(id = x$id, dif_rst$moments$mu, stringsAsFactors = FALSE),
      var = data.frame(id = x$id, dif_rst$moments$var, stringsAsFactors = FALSE),
      covariance = data.frame(id = x$id, dif_rst$moments$covariance, stringsAsFactors = FALSE)
    )
  no_purify$score <- score

  # when purification is used
  if (purify) {
    # verify the criterion for purification
    purify.by <- match.arg(purify.by, c("grdifrs", "grdifr", "grdifs"))

    # create an empty vector and empty data frames
    # to contain the detected DIF items, statistics, and moments
    dif_item <- NULL
    n_df <-
      no_purify$dif_stat %>%
      dplyr::select(dplyr::contains("n.")) %>%
      dplyr::mutate_all(~ {
        NA_integer_
      })
    dif_stat <-
      data.frame(
        id = rep(NA_character_, nrow(x)),
        grdifr = NA,
        grdifs = NA,
        grdifrs = NA,
        p.grdifr = NA,
        p.grdifs = NA,
        p.grdifrs = NA,
        n_df,
        n.iter = NA,
        stringsAsFactors = FALSE
      )
    mu_empty <-
      dif_rst$moments$mu %>%
      dplyr::mutate_all(~ {
        NA_integer_
      })
    var_empty <-
      dif_rst$moments$var %>%
      dplyr::mutate_all(~ {
        NA_integer_
      })
    covar_empty <-
      dif_rst$moments$covariance %>%
      dplyr::mutate_all(~ {
        NA_integer_
      })
    mmt_df <-
      list(
        mu = data.frame(id = rep(NA_character_, nrow(x)), mu_empty, n.iter = NA, stringsAsFactors = FALSE),
        var = data.frame(id = x$id, var_empty, n.iter = NA, stringsAsFactors = FALSE),
        covariance = data.frame(id = x$id, covar_empty, n.iter = NA, stringsAsFactors = FALSE)
      )
    ncol.stat <- ncol(dif_stat)
    ncol.mu <- ncol.var <- ncol(mmt_df$mu)
    ncol.covar <- ncol(mmt_df$covariance)

    # extract the first DIF analysis results
    # and check if at least one DIF item is detected
    dif_item_tmp <- dif_rst$dif_item[[purify.by]]
    log_p_tmp <- dif_rst$log_p[[purify.by]]
    dif_stat_tmp <- dif_rst$dif_stat
    mmt_df_tmp <- no_purify$moments

    # copy the response data, item meta data, and ability estimates
    x_puri <- x
    data_puri <- data
    score_puri <- score

    # create an empty list to contain the post-hoc list
    if (post.hoc) {
      post_df_pury <- list()
    }

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
        dif_stat[del_item, 1:(ncol.stat - 1)] <- dif_stat_tmp[flag_del, ]
        dif_stat[del_item, ncol.stat] <- i - 1
        mmt_df$mu[del_item, 1:(ncol.mu - 1)] <- mmt_df_tmp$mu[flag_del, ]
        mmt_df$var[del_item, 1:(ncol.var - 1)] <- mmt_df_tmp$var[flag_del, ]
        mmt_df$covariance[del_item, 1:(ncol.covar - 1)] <- mmt_df_tmp$covariance[flag_del, ]
        mmt_df$mu[del_item, ncol.mu] <- i - 1
        mmt_df$var[del_item, ncol.var] <- i - 1
        mmt_df$covariance[del_item, ncol.covar] <- i - 1

        # find the examinees who responded to the item to be deleted
        loc_resp <- which(!is.na(data_puri[, flag_del]))

        # refine the leftover items
        item_num <- item_num[-flag_del]

        # post-hoc DIF analysis of the item to be deleted
        if (post.hoc) {
          # post-hoc DIF analyses for all pairs of the groups
          post_rst <- grdif_posthoc(
            x = x_puri[flag_del, , drop = FALSE], data = data_puri[, flag_del, drop = FALSE],
            score = score_puri, group = group, focal.name = focal.name, D = D, alpha = alpha
          )
          singular_post <- union(singular_post, post_rst$singular)

          # re-organize the stat list
          post_df_pury[[i]] <-
            purrr::map_df(.x = post_rst$stats, ~ {
              .x
            }) %>%
            dplyr::mutate(n_iter = i - 1)
        }

        # stop the purification when no item is left to estimate the abilities
        if (length(item_num) == 0L) {
          all_flagged <- TRUE
          break
        }

        # remove the detected DIF item data which has the smallest p-value from the item metadata
        x_puri <- x_puri[-flag_del, , drop = FALSE]

        # remove the detected DIF item data which has the smallest p-value from the response data
        data_puri <- data_puri[, -flag_del, drop = FALSE]

        # if min.resp is not NULL, find the examinees who have the number of responses
        # less than specified value (e.g., 5). Then, replace their all responses with NA
        if (!is.null(min.resp)) {
          n_resp <- Rfast::rowsums(!is.na(data_puri))
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
        if (any(vapply(unique(group[!is.na(group)]), function(g) all(is.na(score_puri[which(group == g)])),
                       logical(1)))) {
          no_score <- TRUE
          score_puri <- score_prev
          dif_stat[item_num, 1:(ncol.stat - 1)] <- dif_stat_tmp[-flag_del, ]
          dif_stat[item_num, ncol.stat] <- i - 1
          mmt_df$mu[item_num, 1:(ncol.mu - 1)] <- mmt_df_tmp$mu[-flag_del, ]
          mmt_df$var[item_num, 1:(ncol.var - 1)] <- mmt_df_tmp$var[-flag_del, ]
          mmt_df$covariance[item_num, 1:(ncol.covar - 1)] <- mmt_df_tmp$covariance[-flag_del, ]
          mmt_df$mu[item_num, ncol.mu] <- i - 1
          mmt_df$var[item_num, ncol.var] <- i - 1
          mmt_df$covariance[item_num, ncol.covar] <- i - 1

          # add the items flagged at the previous iteration, except the item deleted at this one
          dif_item <- c(dif_item, item_num[match(setdiff(dif_item_tmp, flag_del),
                                                 seq_along(log_p_tmp)[-flag_del])])
          break
        }

        # do DIF analysis using the updated ability estimates
        dif_rst_tmp <- grdif_one(
          x = x_puri, data = data_puri, score = score_puri, group = group,
          focal.name = focal.name, D = D, alpha = alpha
        )

        # record the items whose covariance matrix is singular
        singular_id <- union(singular_id, dif_rst_tmp$singular)

        # extract the DIF analysis results
        # and check if at least one DIF item is detected
        dif_item_tmp <- dif_rst_tmp$dif_item[[purify.by]]
        log_p_tmp <- dif_rst_tmp$log_p[[purify.by]]
        dif_stat_tmp <- dif_rst_tmp$dif_stat
        mmt_df_tmp <-
          list(
            mu = data.frame(id = dif_rst_tmp$dif_stat$id, dif_rst_tmp$moments$mu, stringsAsFactors = FALSE),
            var = data.frame(id = dif_rst_tmp$dif_stat$id, dif_rst_tmp$moments$var, stringsAsFactors = FALSE),
            covariance = data.frame(id = dif_rst_tmp$dif_stat$id, dif_rst_tmp$moments$covariance, stringsAsFactors = FALSE)
          )

        # check if a further DIF item is flagged
        if (is.null(dif_item_tmp)) {
          # add the DIF statistics for rest of items
          dif_stat[item_num, 1:(ncol.stat - 1)] <- dif_stat_tmp
          dif_stat[item_num, ncol.stat] <- i
          mmt_df$mu[item_num, 1:(ncol.mu - 1)] <- mmt_df_tmp$mu
          mmt_df$var[item_num, 1:(ncol.var - 1)] <- mmt_df_tmp$var
          mmt_df$covariance[item_num, 1:(ncol.covar - 1)] <- mmt_df_tmp$covariance
          mmt_df$mu[item_num, ncol.mu] <- i
          mmt_df$var[item_num, ncol.var] <- i
          mmt_df$covariance[item_num, ncol.covar] <- i

          break
        }
      }

      # print a message
      if (verbose) {
        cat("", "\n")
      }

      # record the actual number of iterations; the statistics are not recomputed in the
      # iteration that stops the purification
      n_iter <- if (all_flagged || no_score) i - 1L else i

      if (all_flagged) {
        # if every item is flagged, then, return a warning message
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
        dif_stat[item_num, 1:(ncol.stat - 1)] <- dif_stat_tmp
        dif_stat[item_num, ncol.stat] <- i
        mmt_df$mu[item_num, 1:(ncol.mu - 1)] <- mmt_df_tmp$mu
        mmt_df$var[item_num, 1:(ncol.var - 1)] <- mmt_df_tmp$var
        mmt_df$covariance[item_num, 1:(ncol.covar - 1)] <- mmt_df_tmp$covariance
        mmt_df$mu[item_num, ncol.mu] <- i
        mmt_df$var[item_num, ncol.var] <- i
        mmt_df$covariance[item_num, ncol.covar] <- i
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
      with_purify$moments <- mmt_df
      with_purify$dif_item <- sort(dif_item)
      with_purify$n.iter <- n_iter
      with_purify$score <- score_puri
      with_purify$complete <- complete
      if (post.hoc) {
        with_purify$post.hoc <-
          post_df_pury %>%
          dplyr::bind_rows()
      }
    } else {
      # in case when no DIF item is detected from the first DIF analysis results
      with_purify$purify.by <- purify.by
      with_purify$dif_stat <- cbind(no_purify$dif_stat, n.iter = 0)
      with_purify$moments <-
        purrr::map(.x = no_purify$moments, ~ {
          cbind(.x, n.iter = 0)
        })
      with_purify$n.iter <- 0L
      with_purify$complete <- TRUE
    }
  }


  # post-hoc DIF analysis for the non-purified DIF analysis results
  if (post.hoc) {
    # select the flagged items for the post-hoc DIF analyses
    item4post <- no_purify$dif_item

    for (i in 1:3) {
      if (length(item4post[[i]]) != 0) {
        # post-hoc DIF analyses for all pairs of the groups
        post_rst <- grdif_posthoc(
          x = x[item4post[[i]], , drop = FALSE], data = data[, item4post[[i]], drop = FALSE],
          score = score, group = group, focal.name = focal.name, D = D, alpha = alpha
        )
        post.stats <- post_rst$stats
        singular_post <- union(singular_post, post_rst$singular)

        # re-organize the stat list
        post_df_nopury <-
          purrr::map(
            .x = 1:length(item4post[[i]]),
            .f = function(j) {
              purrr::map_df(.x = post.stats, ~ {
                .x[j, ]
              })
            }
          ) %>%
          dplyr::bind_rows()
        no_purify$post.hoc[[i]] <- post_df_nopury
      }
    }
  }

  # warn that the generalized inverse is used for the items with a singular covariance matrix
  if (length(singular_id) > 0L) {
    warning("The covariance matrix of GRDIF_R, GRDIF_S, or GRDIF_RS is singular or nearly singular ",
            "for item(s) ", paste(singular_id, collapse = ", "), ". The statistic of these items is ",
            "computed with a generalized inverse and fewer degrees of freedom.", call. = FALSE)
  }
  if (length(singular_post) > 0L) {
    warning("In the post-hoc analyses, the covariance matrix of RDIF_R and RDIF_S is singular or ",
            "nearly singular for item(s) ", paste(singular_post, collapse = ", "), ". RDIF_RS of ",
            "these items is computed with a generalized inverse and fewer degrees of freedom.",
            call. = FALSE)
  }

  # summarize the results
  rst <- list(no_purify = no_purify, purify = purify, with_purify = with_purify, alpha = alpha)

  # return the DIF detection results
  class(rst) <- "grdif"
  rst
}


# This function conducts the post-hoc RDIF analyses of the given items for all pairs of the groups;
# the reference group is the reference group of its pairs, and the focal group listed earlier in
# focal.name is the reference group of a pair of focal groups
grdif_posthoc <- function(x, data, score, group, focal.name, D, alpha) {
  # list the groups in the order of the reference group followed by the focal groups, as text so
  # that the group values are matched by their labels also when group is a factor
  group_ord <- c(as.character(unique(group[!is.na(group) & !group %in% focal.name])),
                 as.character(focal.name))

  # find all pairs of the groups, and put the group listed earlier first in each pair
  pair.g <- utils::combn(x = sort(unique(group)), m = 2, simplify = FALSE)
  pair.g <- lapply(pair.g, function(pair) pair[order(match(as.character(pair), group_ord))])

  # post-hoc DIF analyses
  post.stats <- list()
  singular <- NULL
  for (p in 1:length(pair.g)) {
    # select a pair of two groups
    pair.tmp <- pair.g[[p]]
    group.pair <- paste(pair.tmp, collapse = " & ")

    # find the location of examinees for the two compared groups
    loc_g1 <- which(group == pair.tmp[1])
    loc_g2 <- which(group == pair.tmp[2])
    loc.tmp <- c(loc_g1, loc_g2)

    # DIF analysis of the two groups
    dif_post.tmp <-
      rdif_one(
        x = x, data = data[loc.tmp, , drop = FALSE], score = score[loc.tmp],
        group = group[loc.tmp], focal.name = pair.tmp[2], D = D, alpha = alpha
      )

    # record the items whose covariance matrix is singular
    singular <- union(singular, dif_post.tmp$singular)

    # add the stats to the list
    post.stats[[p]] <-
      dif_post.tmp$dif_stat %>%
      dplyr::mutate(group.pair = group.pair) %>%
      dplyr::relocate("group.pair", .after = "id")
  }

  # return the results
  list(stats = post.stats, singular = singular)
}


# This function conducts one iteration of the GRDIF analysis
grdif_one <- function(x, data, score, group, focal.name, D = 1, alpha = 0.05) {
  # break down the item metadata into several elements
  elm_item <- breakdown(x)

  # check the unique score categories
  cats <- elm_item$cats

  # check the number of items
  nitem <- length(cats)

  # count N of the focal groups
  n.fg <- length(unique(focal.name))

  ## ---------------------------------
  # compute the GRDIF statistics
  ## ---------------------------------
  # treat the responses of examinees without an ability estimate as missing
  na_score <- is.na(score)
  if (any(na_score)) {
    data[na_score, ] <- NA
    score[na_score] <- 0
  }

  # find the location of examinees for the reference and the focal groups;
  # the examinees without a group are not used
  loc_ref <- which(!is.na(group) & !group %in% focal.name)
  loc_foc <-
    purrr::map(.x = focal.name, ~ {
      which(group == .x)
    })

  # divide the response data into the multiple groups data
  resp_ref <- data[loc_ref, , drop = FALSE]
  resp_foc <- purrr::map(.x = loc_foc, ~ {
    data[.x, , drop = FALSE]
  })

  # count sample sizes
  n_ref <- Rfast::colsums(!is.na(resp_ref))
  n_foc <- purrr::map(.x = resp_foc, ~ {
    Rfast::colsums(!is.na(.x))
  })

  # check if an item has all missing data for either of multiple groups
  all_miss <- sort(unique(c(
    which(n_ref == 0),
    unlist(purrr::map(.x = n_foc, ~ {
      which(.x == 0)
    }))
  )))

  # divide the thetas into multiple groups data
  score_ref <- score[loc_ref]
  score_foc <- purrr::map(.x = loc_foc, ~ {
    score[.x]
  })

  # compute the model-expected item scores (for a dichotomous item, the probability of a correct
  # response) and the model probabilities of score categories
  trace_ref <- trace(elm_item = elm_item, theta = score_ref, D = D, tcc = TRUE)
  trace_foc <-
    purrr::map(.x = score_foc, ~ {
      trace(elm_item = elm_item, theta = .x, D = D, tcc = TRUE)
    })
  extscore_ref <- trace_ref$icc
  extscore_foc <- purrr::map(.x = trace_foc, "icc")
  prob_ref <- trace_ref$prob.cats
  prob_foc <- purrr::map(.x = trace_foc, "prob.cats")

  # replace NA values into the missing data location
  extscore_ref[is.na(resp_ref)] <- NA
  for (g in 1:n.fg) {
    extscore_foc[[g]][is.na(resp_foc[[g]])] <- NA
  }

  # compute the raw residuals across all groups
  resid_ref <- resp_ref - extscore_ref
  resid_foc <-
    purrr::map2(.x = resp_foc, .y = extscore_foc, ~ {
      .x - .y
    })

  # compute the mean raw residuals (mrr) and mean squared residuals (msr) across all groups
  mrr_r <- colMeans(resid_ref, na.rm = TRUE)
  mrr_f <- purrr::map(.x = resid_foc, ~ {
    colMeans(.x, na.rm = TRUE)
  })
  msr_r <- colMeans(resid_ref^2, na.rm = TRUE)
  msr_f <- purrr::map(.x = resid_foc, ~ {
    colMeans(.x^2, na.rm = TRUE)
  })

  # compute the means, variances, and covariances of mrr and msr across all groups
  moments <- resid_moments_mg(
    p_ref = prob_ref, p_foc = prob_foc, n_ref = n_ref, n_foc = n_foc,
    resp_ref = resp_ref, resp_foc = resp_foc, cats = cats, n.fg = n.fg
  )
  mu_rdif <- moments$mu
  var_rdif <- moments$var
  covar_rdif <- moments$covar

  # create a design matrix for grdif_r and grdif_s
  mat.tmp1 <- cbind(rep(-1, n.fg))
  mat.tmp2 <- diag(n.fg)
  dmat_r <- dmat_s <- cbind(mat.tmp1, mat.tmp2)

  # create a design matrix for grdif_rs
  mat.tmp1 <- cbind(rep(c(-1, 0), n.fg))
  mat.tmp2 <- cbind(rep(c(0, -1), n.fg))
  mat.tmp3 <- diag(n.fg * 2)
  dmat_rs <- cbind(mat.tmp1, mat.tmp2, mat.tmp3)

  # compute the chi-square statistics, which are missing for the items that a group did not answer
  chisq_r <- chisq_s <- chisq_rs <- rep(NA_real_, nitem)
  df_r <- df_s <- rep(n.fg, nitem)
  df_rs <- rep(n.fg * 2, nitem)
  for (i in 1:nitem) {
    if (i %in% all_miss) {
      chisq_r[i] <- chisq_s[i] <- chisq_rs[i] <- NaN
    } else {
      # vectors of mrr and msr estimates
      mu_mrr_est <- c(mrr_r[i], purrr::map_dbl(.x = mrr_f, ~ {
        .x[i]
      }))
      mu_msr_est <- c(msr_r[i], purrr::map_dbl(.x = msr_f, ~ {
        .x[i]
      }))

      # vectors of mrr and msr parameters
      mu_mrr_par <- c(
        mu_rdif$ref[[1]][i],
        purrr::map_dbl(.x = purrr::map(.x = mu_rdif$foc, ~ {
          .x[[1]]
        }), ~ {
          .x[i]
        })
      )
      mu_msr_par <- c(
        mu_rdif$ref[[2]][i],
        purrr::map_dbl(.x = purrr::map(.x = mu_rdif$foc, ~ {
          .x[[2]]
        }), ~ {
          .x[i]
        })
      )

      # a vector including both mrr and msr estimates & parameters
      mu_est <- mu_par <- rep(NA, (n.fg + 1) * 2)
      mu_est[2 * (1:(n.fg + 1)) - 1] <- mu_mrr_est
      mu_est[2 * (1:(n.fg + 1))] <- mu_msr_est
      mu_par[2 * (1:(n.fg + 1)) - 1] <- mu_mrr_par
      mu_par[2 * (1:(n.fg + 1))] <- mu_msr_par

      # covariance matrices for grdif_r, grdif_s, and grdif_rs
      var_mrr <- c(
        var_rdif$ref[[1]][i],
        purrr::map_dbl(.x = purrr::map(.x = var_rdif$foc, ~ {
          .x[[1]]
        }), ~ {
          .x[i]
        })
      )
      var_msr <- c(
        var_rdif$ref[[2]][i],
        purrr::map_dbl(.x = purrr::map(.x = var_rdif$foc, ~ {
          .x[[2]]
        }), ~ {
          .x[i]
        })
      )
      covar <- c(covar_rdif$ref[i], purrr::map_dbl(.x = covar_rdif$foc, ~ {
        .x[i]
      }))
      cov_rdifr <- diag(var_mrr)
      cov_rdifs <- diag(var_msr)
      covlist <- list()
      for (g in 1:(n.fg + 1)) {
        # create a var-covariance matrix between mrr and msr
        cov_mat <- array(NA, c(2, 2))

        # replace NAs with the analytically computed covariance
        cov_mat[col(cov_mat) != row(cov_mat)] <- covar[g]

        # replace NAs with the analytically computed variances
        diag(cov_mat) <- c(var_mrr[g], var_msr[g])

        # add the cov-matrix
        covlist[[g]] <- cov_mat
      }
      cov_rdifrs <-
        Matrix::bdiag(covlist) %>%
        data.matrix()

      # multiply the design matrix with statistic vectors
      vec_r <- dmat_r %*% (cbind(mu_mrr_est) - cbind(mu_mrr_par))
      vec_s <- dmat_s %*% (cbind(mu_msr_est) - cbind(mu_msr_par))
      vec_rs <- dmat_rs %*% (cbind(mu_est) - cbind(mu_par))

      # multiply the design matrix with the covariance matrix
      dcov_r <- dmat_r %*% cov_rdifr %*% t(dmat_r)
      dcov_s <- dmat_s %*% cov_rdifs %*% t(dmat_s)
      dcov_rs <- dmat_rs %*% cov_rdifrs %*% t(dmat_rs)

      # compute the chi-square statistics with the inverse or the generalized inverse of the
      # covariance matrices, and the degrees of freedom
      qf_r <- quad_form(cov_mat = dcov_r, dev_vec = vec_r, df = n.fg)
      qf_s <- quad_form(cov_mat = dcov_s, dev_vec = vec_s, df = n.fg)
      qf_rs <- quad_form(cov_mat = dcov_rs, dev_vec = vec_rs, df = n.fg * 2)
      chisq_r[i] <- qf_r$stat
      chisq_s[i] <- qf_s$stat
      chisq_rs[i] <- qf_rs$stat
      df_r[i] <- qf_r$df
      df_s[i] <- qf_s$df
      df_rs[i] <- qf_rs$df
    }
  }

  # calculate p-values for all three statistics
  p_grdifr <- stats::pchisq(chisq_r, df = df_r, lower.tail = FALSE)
  p_grdifs <- stats::pchisq(chisq_s, df = df_s, lower.tail = FALSE)
  p_grdifrs <- stats::pchisq(chisq_rs, df = df_rs, lower.tail = FALSE)

  # find the items whose covariance matrix is singular
  singular_id <- x$id[which(df_r < n.fg | df_s < n.fg | df_rs < n.fg * 2)]

  # create a data.frame of sample sizes
  names(n_foc) <- paste0("n.foc", 1:n.fg)
  n_df <-
    data.frame(n.ref = n_ref, n_foc) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(n.total = sum(dplyr::c_across(cols = dplyr::everything())))

  # create a data frame to contain the results
  stat_df <-
    data.frame(
      id = x$id,
      grdifr = round(chisq_r, 4),
      grdifs = round(chisq_s, 4),
      grdifrs = round(chisq_rs, 4),
      p.grdifr = round(p_grdifr, 4),
      p.grdifs = round(p_grdifs, 4),
      p.grdifrs = round(p_grdifrs, 4),
      n_df,
      stringsAsFactors = FALSE
    )
  rownames(stat_df) <- NULL

  # compute the log p-values to choose the item to be removed in the purification
  log_p <- list(
    grdifr = stats::pchisq(chisq_r, df = df_r, lower.tail = FALSE, log.p = TRUE),
    grdifs = stats::pchisq(chisq_s, df = df_s, lower.tail = FALSE, log.p = TRUE),
    grdifrs = stats::pchisq(chisq_rs, df = df_rs, lower.tail = FALSE, log.p = TRUE)
  )

  # find the flagged items
  dif_item_grdifr <- which(p_grdifr <= alpha)
  dif_item_grdifs <- which(p_grdifs <= alpha)
  dif_item_grdifrs <- which(p_grdifrs <= alpha)
  if (length(dif_item_grdifr) == 0) dif_item_grdifr <- NULL
  if (length(dif_item_grdifs) == 0) dif_item_grdifs <- NULL
  if (length(dif_item_grdifrs) == 0) dif_item_grdifrs <- NULL

  # data frame of the moments
  mu_df <- data.frame(mu_rdif)
  var_df <- data.frame(var_rdif)
  covar_df <- data.frame(covar_rdif)
  names(mu_df) <- c(
    paste0(c("mrr.", "msr."), "ref"),
    paste0(
      rep(c("mrr.", "msr."), n.fg),
      rep(paste0("foc", 1:n.fg), each = 2)
    )
  )
  names(var_df) <- c(
    paste0(c("mrr.", "msr."), "ref"),
    paste0(
      rep(c("mrr.", "msr."), n.fg),
      rep(paste0("foc", 1:n.fg), each = 2)
    )
  )
  names(covar_df) <- c("ref", paste0("foc", 1:n.fg))
  mmt_df <- list(mu = mu_df, var = var_df, covariance = covar_df)

  # summarize the results
  rst <- list(
    dif_stat = stat_df,
    dif_item = list(
      grdifr = dif_item_grdifr,
      grdifs = dif_item_grdifs,
      grdifrs = dif_item_grdifrs
    ),
    moments = mmt_df, alpha = alpha, log_p = log_p, singular = singular_id
  )

  # return the results
  rst
}

# This function computes the mean and variance of the IRT based residual statistics across all groups
# Also, it computes the covariance between mrr and msr across all groups
resid_moments_mg <- function(p_ref, p_foc, n_ref, n_foc, resp_ref, resp_foc, cats, n.fg) {
  # check the number of items
  nitem <- length(cats)

  # count the number of rows for the multiple group response data
  nrow_ref <- nrow(resp_ref)
  nrow_foc <- purrr::map(.x = resp_foc, ~ {
    nrow(.x)
  })

  # create an empty list of two matrices to contain the first moments of raw residuals and squared residuals
  mu_ref <- purrr::map(.x = 1:2, .f = function(x) matrix(NA, nrow = nrow_ref, ncol = nitem))
  mu_foc <-
    purrr::map(
      .x = nrow_foc,
      ~ {
        purrr::map(.x = 1:2, .f = function(x) matrix(NA, nrow = .x, ncol = nitem))
      }
    )

  # create an empty list of two matrices to contain the second moments of raw residuals and squared residuals
  mu2_ref <- mu_ref
  mu2_foc <- mu_foc

  # create an empty matrix to contain the expectation of cube of raw residuals at each theta for all items
  mu_resid3_ref <- matrix(NA, nrow = nrow_ref, ncol = nitem)
  mu_resid3_foc <- purrr::map(.x = nrow_foc, ~ {
    matrix(NA, nrow = .x, ncol = nitem)
  })

  # compute the first and second moments of raw residual and squared residuals for all items
  for (i in 1:nitem) {
    # compute the expected residuals for each score category
    Emat_ref <- matrix(0:(cats[i] - 1), nrow = nrow_ref, ncol = cats[i], byrow = TRUE)
    Emat_foc <-
      purrr::map(
        .x = nrow_foc,
        ~ {
          matrix(0:(cats[i] - 1), nrow = .x, ncol = cats[i], byrow = TRUE)
        }
      )
    exp_resid_ref <-
      Emat_ref - matrix(rowSums(Emat_ref * p_ref[[i]]), nrow = nrow_ref, ncol = cats[i], byrow = FALSE)
    exp_resid_foc <-
      purrr::map(.x = 1:n.fg, ~ {
        Emat_foc[[.x]] - matrix(rowSums(Emat_foc[[.x]] * p_foc[[.x]][[i]]),
          nrow = nrow_foc[[.x]], ncol = cats[i], byrow = FALSE
        )
      })

    # replace NA values into the missing data location
    exp_resid_ref[is.na(resp_ref[, i]), ] <- NA
    for (g in 1:n.fg) {
      exp_resid_foc[[g]][is.na(resp_foc[[g]][, i]), ] <- NA
    }

    # compute the expected values of raw and squared residuals to be used to compute the mean and variance parameters
    value_ref <- vector("list", 2)
    value_ref[[1]] <- exp_resid_ref
    value_ref[[2]] <- exp_resid_ref^2
    resid3_ref <- exp_resid_ref^3
    value_foc <- purrr::map(.x = 1:n.fg, ~ {
      vector("list", 2)
    })
    resid3_foc <- vector("list", n.fg)
    for (g in 1:n.fg) {
      value_foc[[g]][[1]] <- exp_resid_foc[[g]]
      value_foc[[g]][[2]] <- exp_resid_foc[[g]]^2
      resid3_foc[[g]] <- exp_resid_foc[[g]]^3
    }

    # compute the first and second moments across all thetas for the reference and focal groups
    for (j in 1:2) {
      mu_ref[[j]][, i] <- rowSums(value_ref[[j]] * p_ref[[i]], na.rm = FALSE)
      mu2_ref[[j]][, i] <- rowSums((value_ref[[j]])^2 * p_ref[[i]], na.rm = FALSE)
      for (g in 1:n.fg) {
        mu_foc[[g]][[j]][, i] <- rowSums(value_foc[[g]][[j]] * p_foc[[g]][[i]], na.rm = FALSE)
        mu2_foc[[g]][[j]][, i] <- rowSums((value_foc[[g]][[j]])^2 * p_foc[[g]][[i]], na.rm = FALSE)
      }
    }

    # compute the expectation (first moment) of the cube of raw residuals across all thetas for the reference and focal groups
    mu_resid3_ref[, i] <- rowSums(resid3_ref * p_ref[[i]], na.rm = FALSE)
    for (g in 1:n.fg) {
      mu_resid3_foc[[g]][, i] <- rowSums(resid3_foc[[g]] * p_foc[[g]][[i]], na.rm = FALSE)
    }
  }

  # compute the variances across all thetas for the reference and focal groups
  var_ref <- purrr::map2(.x = mu2_ref, .y = mu_ref, .f = function(x, y) x - y^2)
  var_foc <-
    purrr::map(
      .x = 1:n.fg,
      ~ {
        purrr::map2(.x = mu2_foc[[.x]], .y = mu_foc[[.x]], .f = function(x, y) x - y^2)
      }
    )

  # compute the variances of MRR and MSR, respectively, across all groups
  # the variance of a group mean is the sum of the conditional variances divided by the squared sample size
  const_ref <- purrr::map(
    .x = var_ref,
    ~ {
      matrix((1 / n_ref^2), nrow = nrow(.x), ncol = ncol(.x), byrow = TRUE)
    }
  )
  const_foc <-
    purrr::map(
      .x = 1:n.fg,
      ~ {
        purrr::map(
          .x = var_foc[[.x]],
          .f = function(x) matrix((1 / n_foc[[.x]]^2), nrow = nrow(x), ncol = ncol(x), byrow = TRUE)
        )
      }
    )
  var_rdif_ref <- purrr::map2(.x = const_ref, .y = var_ref, ~ {
    colSums(.x * .y, na.rm = TRUE)
  })
  var_rdif_foc <-
    purrr::map(
      .x = 1:n.fg,
      ~ {
        purrr::map2(.x = const_foc[[.x]], .y = var_foc[[.x]], .f = function(x, y) colSums(x * y, na.rm = TRUE))
      }
    )
  names(var_rdif_ref) <- c("mrr", "msr")
  names(var_rdif_foc) <- paste0("foc", 1:n.fg)
  var_rdif_foc <-
    purrr::map(.x = var_rdif_foc, .f = function(x) {
      stats::setNames(x, c("mrr", "msr"))
    })

  # compute the means of MRR and MSR, respectively, across all groups
  mu_rdif_ref <- purrr::map(.x = mu_ref, ~ {
    colMeans(x = .x, na.rm = TRUE)
  })
  mu_rdif_foc <- purrr::map(.x = 1:n.fg, .f = function(x) {
    purrr::map(.x = mu_foc[[x]], ~ {
      colMeans(x = .x, na.rm = TRUE)
    })
  })
  names(mu_rdif_ref) <- c("mrr", "msr")
  names(mu_rdif_foc) <- paste0("foc", 1:n.fg)
  mu_rdif_foc <-
    purrr::map(
      .x = mu_rdif_foc,
      .f = function(x) {
        list(mrr = round(x[[1]], 10), msr = x[[2]])
      }
    )
  mu_rdif_ref[[1]] <- round(mu_rdif_ref[[1]], 10)

  # compute the covariance between mrr and msr across all groups
  covar_ref <- colSums(mu_resid3_ref, na.rm = TRUE) * (1 / n_ref^2)
  covar_foc <-
    purrr::map2(.x = mu_resid3_foc, .y = n_foc, ~ {
      colSums(.x, na.rm = TRUE) * (1 / .y^2)
    })
  names(covar_foc) <- paste0("foc", 1:n.fg)

  # return the results
  mu_rdif <- list(ref = mu_rdif_ref, foc = mu_rdif_foc)
  var_rdif <- list(ref = var_rdif_ref, foc = var_rdif_foc)
  covar_rdif <- list(ref = covar_ref, foc = covar_foc)
  rst <- list(mu = mu_rdif, var = var_rdif, covar = covar_rdif)
  rst
}
