#' IRT Residual-Based Differential Item Functioning (RDIF) Detection Framework
#'
#' This function computes three RDIF statistics for each item: \eqn{RDIF_{R}},
#' \eqn{RDIF_{S}}, and \eqn{RDIF_{RS}} (Lim & Choe, 2023; Lim et al., 2022).
#' \eqn{RDIF_{R}} primarily captures differences in raw residuals between two
#' groups, which are typically associated with uniform DIF. \eqn{RDIF_{S}}
#' primarily captures differences in squared residuals, which are typically
#' associated with nonuniform DIF. \eqn{RDIF_{RS}} jointly considers both types
#' of differences and is capable of detecting both uniform and nonuniform DIF.
#'
#' @inheritParams irtfit
#' @inheritParams est_score
#' @param score A numeric vector of examinees' ability estimates (theta values),
#'   in the same row order as `data`. The estimates should be based on the item
#'   parameters in `x` (see **Details**). If `NULL`, abilities are estimated
#'   internally with [irtQ::est_score()] using `method`, `range`, `norm.prior`,
#'   `nquad`, `weights`, and `ncore`. When `purify = TRUE`, abilities are
#'   re-estimated internally at every purification iteration, so a supplied
#'   `score` is used only in the initial analysis. A missing value in `score`
#'   excludes the examinee from the analysis, also during purification, with a
#'   warning. The `est_item`
#'   method has no `score` argument; it uses the abilities stored in the object.
#'   Default is `NULL`.
#' @param group A numeric or character vector indicating examinees' group
#'   membership. The length of the vector must match the number of rows in the
#'   response data matrix.
#' @param focal.name A single numeric or character value specifying the focal
#'   group. For instance, given `group = c(0, 1, 0, 1, 1)` and '1' indicating
#'   the focal group, set `focal.name = 1`.
#' @param item.skip A numeric vector of item positions (row numbers of `x`), or a
#'   logical vector with one value per item, identifying the items to exclude
#'   from the analysis. Skipped items still contribute to the ability
#'   estimates, but their results (statistics, p-values, and moments) are set to
#'   `NA`, and they are never flagged or removed during purification. If `NULL`,
#'   all items are analyzed. Default is `NULL`.
#' @param alpha A numeric value specifying the significance level (\eqn{\alpha})
#'   of the tests. An item is flagged by a statistic when its p-value is less
#'   than or equal to `alpha`. The comparison uses unrounded p-values, although
#'   the p-values reported in `dif_stat` are rounded to four decimal places.
#'   Default is `0.05`.
#' @param missing  A value indicating missing responses in the data set. Default
#'   is `NA`.
#' @param purify Logical. Indicates whether to apply a purification procedure.
#'   Default is `FALSE`.
#' @param purify.by A character string specifying the RDIF statistic used for
#'   purification: "rdifrs" for \eqn{RDIF_{RS}}, "rdifr" for \eqn{RDIF_{R}}, or
#'   "rdifs" for \eqn{RDIF_{S}}. Used only when `purify = TRUE`. Default is
#'   `"rdifrs"`.
#' @param max.iter A positive integer specifying the maximum number of
#'   iterations allowed for the purification process. Default is `10`. Each
#'   iteration removes one item, so `max.iter` should be at least as large as
#'   the number of items expected to be flagged. If the limit is reached while
#'   flagged items remain, a warning is issued, the `complete` element of the
#'   purification results is `FALSE`, and the items flagged in the last
#'   iteration are added to the flagged items. This argument is not passed to
#'   [irtQ::est_score()].
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
#'   pattern. Default is `"ML"`.
#' @param range A numeric vector of length two specifying the lower and upper
#'   bounds of the ability scale. This is used for the following scoring
#'   methods: `"ML"`, `"WL"`, and `"MAP"`. Default is `c(-5, 5)`.
#' @param norm.prior A numeric vector of length two specifying the mean and
#'   standard deviation of the normal prior distribution. These values are used
#'   to generate the Gaussian quadrature points and weights. Ignored if `method`
#'   is `"ML"` or `"WL"`. Default is `c(0, 1)`.
#' @param nquad An integer indicating the number of Gaussian quadrature points
#'   to be generated from the normal prior distribution. Used only when `method`
#'   is `"EAP"`. Ignored for `"ML"`, `"WL"`, and `"MAP"`. Default is 41.
#' @param weights A two-column matrix or data frame containing the quadrature
#'   points (in the first column) and their corresponding weights (in the second
#'   column) for the latent variable prior distribution. The weights and points
#'   can be conveniently generated using the function [irtQ::gen.weight()].
#'
#'   If `NULL` and `method = "EAP"`, default quadrature values are generated
#'   based on the `norm.prior` and `nquad` arguments. Ignored if `method` is
#'   `"ML"`, `"WL"`, or `"MAP"`.
#' @param ncore An integer specifying the number of logical CPU cores used by
#'   [irtQ::est_score()] when abilities are estimated internally. The test
#'   statistics themselves are not computed in parallel. Default is `1`.
#' @param verbose Logical. If `TRUE`, the progress of the purification procedure
#'   is printed to the console. Default is `TRUE`.
#' @param ... Additional arguments passed to [irtQ::est_score()] when abilities
#'   are estimated internally, for example `tol` or `se`. The `est_score()`
#'   arguments `max.iter` and `missing` cannot be passed this way, because the
#'   function uses arguments with these names for its own purposes.
#'
#' @details The RDIF framework (Lim & Choe, 2023; Lim et al., 2022) consists of
#'   three IRT residual-based statistics: \eqn{RDIF_{R}}, \eqn{RDIF_{S}}, and
#'   \eqn{RDIF_{RS}}. For each item, a residual is the observed item score minus
#'   the model-expected item score evaluated at the examinee's ability estimate.
#'   \eqn{RDIF_{R}} is the difference between the focal and reference groups in
#'   the mean raw residual, and \eqn{RDIF_{S}} is the difference in the mean
#'   squared residual. Under the null hypothesis that the item has no DIF,
#'   \eqn{RDIF_{R}} and \eqn{RDIF_{S}} asymptotically follow normal
#'   distributions whose means and variances are computed analytically from the
#'   model-predicted probabilities (Lim et al., 2022, Equations 9, 10, 12, and
#'   13). The null mean of \eqn{RDIF_{R}} is zero, whereas the null mean of
#'   \eqn{RDIF_{S}} is generally not zero. Each statistic is standardized by
#'   subtracting its null mean and dividing by its null standard deviation
#'   (`z.rdifr` and `z.rdifs`), and a two-sided p-value is obtained from the
#'   standard normal distribution. \eqn{RDIF_{RS}} is the quadratic form of the
#'   vector of \eqn{RDIF_{R}} and \eqn{RDIF_{S}}, centered at their null means,
#'   in the inverse of their null covariance matrix. Under the null hypothesis,
#'   \eqn{RDIF_{RS}} asymptotically follows a \eqn{\chi^{2}} distribution with 2
#'   degrees of freedom (Lim et al., 2022, Equation 20).
#'
#'   In rare cases, the covariance matrix of \eqn{RDIF_{R}} and \eqn{RDIF_{S}}
#'   is singular or nearly singular (its reciprocal condition number is 1e-10 or
#'   less), for example when all examinees who responded to an item have the
#'   same probability of a correct response, so that \eqn{RDIF_{S}} is a linear
#'   function of \eqn{RDIF_{R}}. Then \eqn{RDIF_{RS}} is computed with the
#'   Moore-Penrose generalized inverse of the covariance matrix and compared with
#'   a chi-square distribution whose degrees of freedom equal the rank of the
#'   matrix (Moore, 1977), and a warning names the item. This handling is not
#'   part of the original method (Lim et al., 2022) and is added in irtQ. In
#'   most data the matrix is not singular, and \eqn{RDIF_{RS}} has two degrees of
#'   freedom.
#'
#'   [irtQ::rdif()] accepts both dichotomous and polytomous items. For a
#'   polytomous item, the residual is the observed item score minus the
#'   model-expected item score, and the null means, variances, and covariance are
#'   computed from the model-predicted category probabilities in the same way,
#'   so the three statistics assess net DIF (Jung & Lim, 2026; Lim, Malatesta, &
#'   Lee, 2024). To evaluate global DIF in polytomous items, use
#'   [irtQ::crdif()].
#'
#'   To compute the RDIF statistics, the [irtQ::rdif()] function requires:
#'   (1) item parameter estimates obtained from aggregate data (regardless
#'   of group membership), (2) examinees' ability estimates (e.g., ML), and
#'   (3) examinees' item response data. Note that the ability estimates must
#'   be based on the aggregate-data item parameters. The item parameter estimates
#'   should be provided in the `x` argument, the ability estimates in the `score`
#'   argument, and the response data in the `data` argument. If ability
#'   estimates are not provided (i.e., `score = NULL`), [irtQ::rdif()] will
#'   estimate them automatically using the scoring method specified via the
#'   `method` argument (e.g., `method = "ML"`). When `x` is an object of class
#'   `est_irt`, the item parameter estimates, the response data, and the scaling
#'   constant `D` are taken from the object. When `x` is an object of class
#'   `est_item`, the ability values stored in the object are also used as
#'   `score`.
#'
#'   The `group` argument should be a vector containing exactly two distinct
#'   values (either numeric or character), representing the reference and focal
#'   groups, and the function stops when it contains more than two. Its length
#'   must match the number of rows in the response data, where each element
#'   corresponds to an examinee. Once `group` is specified, a single numeric or
#'   character value must be provided in the `focal.name` argument to indicate
#'   which level in `group` represents the focal group. Use [irtQ::grdif()] to
#'   compare more than two groups. Examinees with `NA` in `group` are excluded
#'   from both groups.
#'
#'   Examinees whose ability estimate is `NA`, such as examinees without any item
#'   response or examinees with `NA` in a supplied `score`, are excluded from the
#'   computation of the RDIF statistics. A warning reports the number of
#'   examinees with item responses who are excluded because of a missing value in
#'   `score`; these examinees stay excluded during purification. During
#'   purification, an examinee whose remaining responses are all removed with
#'   the flagged items, or who has fewer than `min.resp` responses to the
#'   remaining items, has no ability estimate and is excluded, and one warning
#'   at the end reports the number of such examinees. The sample sizes `n.ref` and
#'   `n.foc` count only the examinees who are used.
#'
#'   If no examinee of one of the two groups responded to an item, the
#'   statistics and p-values of the item are `NaN`, and the item is not flagged.
#'
#'   Similar to other DIF detection approaches, the RDIF framework supports an
#'   iterative purification procedure (Lim et al., 2022). When `purify = TRUE`,
#'   the statistic specified in `purify.by` drives the procedure (e.g.,
#'   `purify.by = "rdifrs"`). At each iteration, the flagged item with the
#'   smallest p-value of the `purify.by` statistic is removed (the p-values are
#'   compared on the log scale), the abilities are re-estimated from the
#'   remaining items with the scoring method specified in `method`, and the RDIF
#'   statistics of the remaining items are recomputed. A supplied `score` is
#'   therefore used only in the initial analysis. The
#'   procedure stops when no remaining item is flagged or when `max.iter`
#'   iterations have been performed, where `max.iter` must be a single whole
#'   number of at least 1. If no item is flagged in the initial analysis, no
#'   iteration is performed. If every item is flagged, the procedure stops with
#'   a warning, because no item is left to estimate the abilities, and
#'   `complete` is `FALSE`. The procedure also stops with a warning and
#'   `complete = FALSE` when, after an item is removed, no examinee of the
#'   reference group or of the focal group has an ability estimate left (see
#'   `min.resp` below); the remaining items then keep the results of the
#'   previous iteration, and the items flagged in that iteration are added to
#'   the flagged items.
#'
#'   Scoring based on a small number of item responses can lead to large
#'   standard errors, potentially reducing the accuracy of DIF detection in the
#'   RDIF framework. The `min.resp` argument excludes such examinees. For
#'   example, if `min.resp = 5`, examinees who responded to fewer than five items
#'   (but at least one) have all their responses treated as missing (i.e., `NA`)
#'   before the abilities are estimated, in the initial analysis and at every
#'   purification iteration. As a result, their ability estimates are missing,
#'   and they are not used in the computation of the RDIF statistics. When
#'   `score` is supplied, the scores of these examinees are also treated as
#'   missing in the initial analysis. During purification, the count is based on
#'   the items that remain after the flagged items are removed, and an examinee
#'   excluded once stays excluded. Because the count is based on the remaining
#'   items, a large `min.resp` can leave a group without any ability estimate in
#'   a later iteration, and the purification then stops with a warning. If
#'   `min.resp = NULL`, a score is computed for any examinee with at least one
#'   valid item response.
#'
#'   The print method of this class has an argument `digits` that is not used;
#'   the statistics are printed rounded to three decimal places.
#'
#' @return This function returns an object of class `"rdif"`, which is a list
#' with the following five components:
#'
#' \item{no_purify}{A list of sub-objects containing the results of DIF analysis
#' without applying a purification procedure. The sub-objects include:
#'   \describe{
#'     \item{dif_stat}{A data frame with one row per item and the columns `id`,
#'     `rdifr`, `z.rdifr`, `rdifs`, `z.rdifs`, `rdifrs`, `p.rdifr`, `p.rdifs`,
#'     `p.rdifrs`, `n.ref`, `n.foc`, and `n.total`: the item ID, the
#'     \eqn{RDIF_{R}} statistic and its standardized value, the \eqn{RDIF_{S}}
#'     statistic and its standardized value, the \eqn{RDIF_{RS}} statistic, the
#'     p-values of the three statistics (two-sided for \eqn{RDIF_{R}} and
#'     \eqn{RDIF_{S}}), the numbers of examinees in the reference and focal
#'     groups who responded to the item, and their sum. Statistics and p-values
#'     are rounded to four decimal places. \eqn{RDIF_{RS}} has no standardized
#'     value because it is a \eqn{\chi^{2}}-based statistic.}
#'     \item{moments}{A data frame with the columns `id`, `mu.rdifr`,
#'     `sigma.rdifr`, `mu.rdifs`, `sigma.rdifs`, and `covariance`: the item ID,
#'     the null mean and standard deviation of \eqn{RDIF_{R}}, the null mean and
#'     standard deviation of \eqn{RDIF_{S}}, and the null covariance between
#'     \eqn{RDIF_{R}} and \eqn{RDIF_{S}}. These moments are computed
#'     analytically under the null hypothesis of no DIF and are used to
#'     standardize the statistics.}
#'     \item{dif_item}{A list with the elements `rdifr`, `rdifs`, and `rdifrs`,
#'     each giving the positions (rows of `x`) of the items flagged by the
#'     corresponding statistic, or `NULL` if no item is flagged.}
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
#'     \item{purify.by}{A character string indicating the RDIF statistic used for
#'     purification. Possible values are "rdifr", "rdifs", and "rdifrs",
#'     corresponding to \eqn{RDIF_{R}}, \eqn{RDIF_{S}}, and \eqn{RDIF_{RS}},
#'     respectively.}
#'     \item{dif_stat}{A data frame with the same columns as
#'     `no_purify$dif_stat` and an additional column `n.iter`. For an item
#'     removed during purification, the row reports the statistics from the
#'     iteration in which the item was removed, and `n.iter` is that iteration
#'     minus 1 (0 for an item removed on the basis of the initial analysis). For
#'     the other items, the row reports the statistics from the last iteration,
#'     and `n.iter` is the number of that iteration. If no item is flagged in the
#'     initial analysis, this is `no_purify$dif_stat` with `n.iter = 0`.}
#'     \item{moments}{A data frame with the same columns as
#'     `no_purify$moments` and an additional column `n.iter`, defined in the same
#'     way as in `dif_stat`.}
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
#' @seealso [irtQ::est_irt()], [irtQ::est_item()], [irtQ::simdat()],
#'   [irtQ::shape_df()], [irtQ::est_score()]
#'
#' @references Bock, R. D., & Mislevy, R. J. (1982). Adaptive EAP estimation of
#'   ability in a microcomputer environment. *Applied Psychological Measurement,
#'   6*(4), 431-444. \doi{10.1177/014662168200600405}.
#'
#'   Hambleton, R. K., Swaminathan, H., & Rogers, H. J. (1991). *Fundamentals of
#'   item response theory*. Newbury Park, CA: Sage.
#'
#'   Jung, H., & Lim, H. (2026, April). Detecting global and net DIF in
#'   polytomous items using RDIF. Paper presented at the annual meeting of the
#'   National Council on Measurement in Education, Los Angeles, CA.
#'
#'   Lim, H., & Choe, E. M. (2023). Detecting differential item functioning in
#'   CAT using IRT residual DIF approach. *Journal of Educational Measurement,
#'   60*(4), 626-650. \doi{10.1111/jedm.12366}.
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
#'   Warm, T. A. (1989). Weighted likelihood estimation of ability in item
#'   response theory. *Psychometrika, 54*(3), 427-450.
#'   \doi{10.1007/BF02294627}.
#'
#' @examples
#' \donttest{
#' # Load required package
#' library("dplyr")
#'
#' ## Uniform DIF detection
#' ###############################################
#' # (1) Generate data with known uniform DIF
#' ###############################################
#'
#' # Import the "-prm.txt" output file from flexMIRT
#' flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
#'
#' # Select 36 non-DIF items using the 3PLM model
#' par_nstd <-
#'   bring.flexmirt(file = flex_sam, "par")$Group1$full_df %>%
#'   dplyr::filter(.data$model == "3PLM") %>%
#'   dplyr::filter(dplyr::row_number() %in% 1:36) %>%
#'   dplyr::select(1:6)
#' par_nstd$id <- paste0("nondif", 1:36)
#'
#' # Generate 4 new DIF items for the reference group
#' difpar_ref <-
#'   shape_df(
#'     par.drm = list(a = c(0.8, 1.5, 0.8, 1.5), b = c(0.0, 0.0, -0.5, -0.5), g = rep(0.15, 4)),
#'     item.id = paste0("dif", 1:4), cats = 2, model = "3PLM"
#'   )
#'
#' # Add uniform DIF by shifting the b-parameters for the focal group
#' difpar_foc <-
#'   difpar_ref %>%
#'   dplyr::mutate_at(.vars = "par.2", .funs = function(x) x + rep(0.7, 4))
#'
#' # Combine the DIF and non-DIF items for both reference and focal groups
#' # Therefore, the first 4 items exhibit uniform DIF
#' par_ref <- rbind(difpar_ref, par_nstd)
#' par_foc <- rbind(difpar_foc, par_nstd)
#'
#' # Generate true ability values
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
#' #     from the combined response data
#' ###############################################
#'
#' # Estimate item parameters
#' est_mod <- est_irt(data = data, D = 1, model = "3PLM")
#' est_par <- est_mod$par.est
#'
#' # Estimate ability parameters using ML
#' score <- est_score(x = est_par, data = data, method = "ML")$est.theta
#'
#' ###############################################
#' # (3) Perform DIF analysis
#' ###############################################
#'
#' # Define group membership: 1 = focal group
#' group <- c(rep(0, 500), rep(1, 500))
#'
#' # (a)-1 Compute RDIF statistics with provided ability scores
#' #       (no purification)
#' dif_nopuri_1 <- rdif(
#'   x = est_par, data = data, score = score,
#'   group = group, focal.name = 1, D = 1, alpha = 0.05
#' )
#' print(dif_nopuri_1)
#'
#' # (a)-2 Compute RDIF statistics without providing ability scores
#' #       (no purification)
#' dif_nopuri_2 <- rdif(
#'   x = est_par, data = data, score = NULL,
#'   group = group, focal.name = 1, D = 1, alpha = 0.05,
#'   method = "ML"
#' )
#' print(dif_nopuri_2)
#'
#' # (b) Compute RDIF statistics with purification based on RDIF(RS)
#' dif_puri_rs <- rdif(
#'   x = est_par, data = data, score = score,
#'   group = group, focal.name = 1, D = 1, alpha = 0.05,
#'   purify = TRUE, purify.by = "rdifrs"
#' )
#' print(dif_puri_rs)
#'
#' # purify.by = "rdifr" or "rdifs" can be used in the same way.
#'
#' }
#'
#' @order 1
#' @export
rdif <- function(x, ...) UseMethod("rdif")

#' @describeIn rdif Default method for a data frame `x` containing item metadata.
#'
#' @order 2
#' @export
rdif.default <- function(x,
                         data,
                         score = NULL,
                         group,
                         focal.name,
                         item.skip = NULL,
                         D = 1,
                         alpha = 0.05,
                         missing = NA,
                         purify = FALSE,
                         purify.by = c("rdifrs", "rdifr", "rdifs"),
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
  rst <- rdif_main(
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

#' @describeIn rdif Method for an object of class `est_irt` created by
#' [irtQ::est_irt()]. The response data and the scaling factor `D` are taken from
#' the object.
#'
#' @order 3
#' @export
#'
rdif.est_irt <- function(x,
                         score = NULL,
                         group,
                         focal.name,
                         item.skip = NULL,
                         alpha = 0.05,
                         missing = NA,
                         purify = FALSE,
                         purify.by = c("rdifrs", "rdifr", "rdifs"),
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
  rst <- rdif_main(
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

#' @describeIn rdif Method for an object of class `est_item` created by
#' [irtQ::est_item()]. The response data, ability estimates, and the scaling
#' factor `D` are taken from the object.
#'
#' @order 4
#' @export
#'
rdif.est_item <- function(x,
                          group,
                          focal.name,
                          item.skip = NULL,
                          alpha = 0.05,
                          missing = NA,
                          purify = FALSE,
                          purify.by = c("rdifrs", "rdifr", "rdifs"),
                          max.iter = 10,
                          min.resp = NULL,
                          method = "ML",
                          range = c(-5, 5),
                          norm.prior =
                            c(0, 1),
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
  rst <- rdif_main(
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


# This function conducts the DIF analysis shared by the methods of rdif(): it prepares the
# data and ability estimates, runs the analysis once, and applies the purification procedure
rdif_main <- function(x, data, score, group, focal.name, item.skip, D, alpha, missing,
                      purify, purify.by, max.iter, min.resp, method, range, norm.prior,
                      nquad, weights, ncore, verbose, ...) {

  ## ----------------------------------
  ## (1) prepare DIF analysis
  ## ----------------------------------
  # confirm and correct all item metadata information
  x <- confirm_df(x)

  # record whether the ability estimates are supplied
  score_given <- !is.null(score)

  # check the inputs and prepare the response data and ability estimates
  prep <- dif_prepare(
    x = x, data = data, score = score, group = group, focal.name = focal.name,
    item.skip = item.skip, D = D, alpha = alpha, missing = missing, purify = purify,
    max.iter = max.iter, min.resp = min.resp, method = method, range = range,
    norm.prior = norm.prior, nquad = nquad, weights = weights, ncore = ncore, fn = "rdif()", ...
  )
  data <- prep$data
  score <- prep$score
  item.skip <- prep$item.skip

  # a) when no purification is set
  # do only one iteration of DIF analysis
  dif_rst <-
    rdif_one(x = x, data = data, score = score, group = group,
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
  no_purify$moments <- rdif_moments(dif_rst)
  no_purify$score <- score

  # when purification is used
  if (purify) {
    # verify the criterion for purification
    purify.by <- match.arg(purify.by, c("rdifrs", "rdifr", "rdifs"))

    # create an empty vector and empty data frames
    # to contain the detected DIF items, statistics, and moments
    dif_item <- NULL
    dif_stat <-
      data.frame(
        id = rep(NA_character_, nrow(x)),
        rdifr = NA, z.rdifr = NA,
        rdifs = NA, z.rdifs = NA,
        rdifrs = NA, p.rdifr = NA,
        p.rdifs = NA, p.rdifrs = NA,
        n.ref = NA, n.foc = NA, n.total = NA, n.iter = NA,
        stringsAsFactors = FALSE
      )
    mmt_df <-
      data.frame(
        id = rep(NA_character_, nrow(x)), mu.rdifr = NA, sigma.rdifr = NA,
        mu.rdifs = NA, sigma.rdifs = NA, covariance = NA, n.iter = NA,
        stringsAsFactors = FALSE
      )

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
        dif_stat[del_item, 1:12] <- dif_stat_tmp[flag_del, ]
        dif_stat[del_item, 13] <- i - 1
        mmt_df[del_item, 1:6] <- mmt_df_tmp[flag_del, ]
        mmt_df[del_item, 7] <- i - 1

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
        # and drop the item parameter columns that no remaining item uses
        x_puri <- trim_par_cols(x_puri[-flag_del, , drop = FALSE])

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
          dif_stat[item_num, 1:12] <- dif_stat_tmp[-flag_del, ]
          dif_stat[item_num, 13] <- i - 1
          mmt_df[item_num, 1:6] <- mmt_df_tmp[-flag_del, ]
          mmt_df[item_num, 7] <- i - 1

          # add the items flagged at the previous iteration, except the item deleted at this one
          dif_item <- c(dif_item, item_num[match(setdiff(dif_item_tmp, flag_del),
                                                 seq_along(log_p_tmp)[-flag_del])])
          break
        }

        # do DIF analysis using the updated ability estimates
        dif_rst_tmp <- rdif_one(
          x = x_puri, data = data_puri, score = score_puri, group = group,
          focal.name = focal.name, item.skip = item.skip.puri, D = D,
          alpha = alpha
        )

        # record the items whose covariance matrix is singular
        singular_id <- union(singular_id, dif_rst_tmp$singular)

        # extract the DIF analysis results
        # and check if at least one DIF item is detected
        dif_item_tmp <- dif_rst_tmp$dif_item[[purify.by]]
        log_p_tmp <- dif_rst_tmp$log_p[[purify.by]]
        dif_stat_tmp <- dif_rst_tmp$dif_stat
        mmt_df_tmp <- rdif_moments(dif_rst_tmp)

        # check if a further DIF item is flagged
        if (is.null(dif_item_tmp)) {
          # add the DIF statistics for rest of items
          dif_stat[item_num, 1:12] <- dif_stat_tmp
          dif_stat[item_num, 13] <- i
          mmt_df[item_num, 1:6] <- mmt_df_tmp
          mmt_df[item_num, 7] <- i

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
        dif_stat[item_num, 1:12] <- dif_stat_tmp
        dif_stat[item_num, 13] <- i
        mmt_df[item_num, 1:6] <- mmt_df_tmp
        mmt_df[item_num, 7] <- i
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
    } else {
      # in case when no DIF item is detected from the first DIF analysis results
      with_purify$purify.by <- purify.by
      with_purify$dif_stat <- cbind(no_purify$dif_stat, n.iter = 0)
      with_purify$moments <- cbind(no_purify$moments, n.iter = 0)
      with_purify$n.iter <- 0L
      with_purify$complete <- TRUE
    }
  }

  # warn that the generalized inverse is used for the items with a singular covariance matrix
  if (length(singular_id) > 0L) {
    warning("The covariance matrix of RDIF_R and RDIF_S is singular or nearly singular for item(s) ",
            paste(singular_id, collapse = ", "), ". RDIF_RS of these items is computed ",
            "with a generalized inverse and fewer degrees of freedom.", call. = FALSE)
  }

  # summarize the results
  rst <- list(no_purify = no_purify, purify = purify, with_purify = with_purify, alpha = alpha)

  # return the DIF detection results
  class(rst) <- "rdif"
  rst
}


# This function creates the data frame of the null moments of the two statistics
rdif_moments <- function(dif_rst) {
  # combine the means, standard deviations, and covariance
  mmt_df <- data.frame(
    id = dif_rst$dif_stat$id,
    dif_rst$moments$rdifr[, c(1, 3)],
    dif_rst$moments$rdifs[, c(1, 3)],
    dif_rst$covariance, stringsAsFactors = FALSE
  )
  names(mmt_df) <- c("id", "mu.rdifr", "sigma.rdifr", "mu.rdifs", "sigma.rdifs", "covariance")
  mmt_df
}


# This function conducts one iteration of DIF analysis using the IRT residual based statistics
#' @import dplyr
rdif_one <- function(x,
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
  cats <- elm_item$cats

  # check the number of items
  nitem <- length(cats)

  ## ---------------------------------
  # compute the two statistics
  ## ---------------------------------
  # treat the responses of examinees without an ability estimate as missing
  na_score <- is.na(score)
  if (any(na_score)) {
    data[na_score, ] <- NA
    score[na_score] <- 0
  }

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

  # compute the model-expected item scores (for a dichotomous item, the probability of a correct
  # response) and the model probabilities of score categories
  trace_ref <- trace(elm_item = elm_item, theta = score_ref, D = D, tcc = TRUE)
  trace_foc <- trace(elm_item = elm_item, theta = score_foc, D = D, tcc = TRUE)
  extscore_ref <- trace_ref$icc
  extscore_foc <- trace_foc$icc
  prob_ref <- trace_ref$prob.cats
  prob_foc <- trace_foc$prob.cats

  # replace NA values into the missing data location
  extscore_ref[is.na(resp_ref)] <- NA
  extscore_foc[is.na(resp_foc)] <- NA

  # compute the raw residuals
  resid_ref <- resp_ref - extscore_ref
  resid_foc <- resp_foc - extscore_foc

  # compute the residual-based DIF statistic
  rdifr <- colMeans(resid_foc, na.rm = TRUE) - colMeans(resid_ref, na.rm = TRUE) # the difference of mean raw residuals (DMRR)
  rdifs <- colMeans(resid_foc^2, na.rm = TRUE) - colMeans(resid_ref^2, na.rm = TRUE) # the difference of mean squared residuals (DMSR)

  # compute the means and variances of the two statistics for the hypothesis testing
  moments <- resid_moments(
    p_ref = prob_ref, p_foc = prob_foc, n_ref = n_ref, n_foc = n_foc,
    resp_ref = resp_ref, resp_foc = resp_foc, cats = cats
  )
  moments_rdifr <- moments$rdifr
  moments_rdifs <- moments$rdifs
  covar <- moments$covariance

  # compute the chi-square statistics and their degrees of freedom
  chisq <- rep(NA_real_, nitem)
  df_chisq <- rep(2, nitem)
  for (i in 1:nitem) {
    if (i %in% all_miss) {
      chisq[i] <- NaN
    } else {
      # create a var-covariance matrix between rdifr and rdifs
      cov_mat <- array(NA, c(2, 2))

      # replace NAs with the analytically computed covariance
      cov_mat[col(cov_mat) != row(cov_mat)] <- covar[i]

      # replace NAs with the analytically computed variances
      diag(cov_mat) <- c(moments_rdifr[i, 2], moments_rdifs[i, 2])

      # create a vector of mean rdifr and mean rdifs
      mu_vec <- cbind(moments_rdifr[i, 1], moments_rdifs[i, 1])

      # create a vector of rdifr and rdifs
      est_mu_vec <- cbind(rdifr[i], rdifs[i])

      # compute the chi-square statistic with the inverse or the generalized inverse of the
      # covariance matrix, and the degrees of freedom
      qf <- quad_form(cov_mat = cov_mat, dev_vec = t(est_mu_vec - mu_vec), df = 2)
      chisq[i] <- qf$stat
      df_chisq[i] <- qf$df
    }
  }

  # find the items whose covariance matrix is singular
  singular_id <- x$id[which(df_chisq < 2 & !(seq_len(nitem) %in% item.skip))]

  # standardize the two statistics of rdifr and rdifs
  z_stat_rdifr <- (rdifr - moments_rdifr$mu) / moments_rdifr$sigma
  z_stat_rdifs <- (rdifs - moments_rdifs$mu) / moments_rdifs$sigma

  # calculate p-values for all three statistics
  p_rdifr <- 2 * stats::pnorm(q = abs(z_stat_rdifr), mean = 0, sd = 1, lower.tail = FALSE)
  p_rdifs <- 2 * stats::pnorm(q = abs(z_stat_rdifs), mean = 0, sd = 1, lower.tail = FALSE)
  p_rdifrs <- stats::pchisq(chisq, df = df_chisq, lower.tail = FALSE)

  # compute total sample size
  n_total <- n_foc + n_ref

  # create a data frame to contain the results
  stat_df <-
    data.frame(
      id = x$id,
      rdifr = round(rdifr, 4), z.rdifr = round(z_stat_rdifr, 4),
      rdifs = round(rdifs, 4), z.rdifs = round(z_stat_rdifs, 4),
      rdifrs = round(chisq, 4),
      p.rdifr = round(p_rdifr, 4), p.rdifs = round(p_rdifs, 4), p.rdifrs = round(p_rdifrs, 4),
      n.ref = n_ref, n.foc = n_foc, n.total = n_total, stringsAsFactors = FALSE
    )
  rownames(stat_df) <- NULL

  # when there are items that should be skipped for the DIF analysis
  # insert NAs to the corresponding results of the items
  if (!is.null(item.skip)) {
    stat_df[item.skip, 2:9] <- NA
    p_rdifr[item.skip] <- NA
    p_rdifs[item.skip] <- NA
    p_rdifrs[item.skip] <- NA
    moments_rdifr[item.skip, ] <- NA
    moments_rdifs[item.skip, ] <- NA
    covar[item.skip] <- NA
  }

  # compute the log p-values to choose the item to be removed in the purification
  log_p <- list(
    rdifr = log(2) + stats::pnorm(q = -abs(z_stat_rdifr), mean = 0, sd = 1, log.p = TRUE),
    rdifs = log(2) + stats::pnorm(q = -abs(z_stat_rdifs), mean = 0, sd = 1, log.p = TRUE),
    rdifrs = stats::pchisq(chisq, df = df_chisq, lower.tail = FALSE, log.p = TRUE)
  )

  # find the flagged items using the unrounded p-values
  dif_item_rdifr <- as.integer(which(p_rdifr <= alpha))
  dif_item_rdifs <- as.integer(which(p_rdifs <= alpha))
  dif_item_rdifrs <- which(p_rdifrs <= alpha)
  if(length(dif_item_rdifr) == 0) dif_item_rdifr <- NULL
  if(length(dif_item_rdifs) == 0) dif_item_rdifs <- NULL
  if(length(dif_item_rdifrs) == 0) dif_item_rdifrs <- NULL

  # summarize the results
  rst <- list(
    dif_stat = stat_df,
    dif_item = list(rdifr = dif_item_rdifr, rdifs = dif_item_rdifs, rdifrs = dif_item_rdifrs),
    moments = list(rdifr = moments_rdifr, rdifs = moments_rdifs), covariance = covar, alpha = alpha,
    log_p = log_p, singular = singular_id
  )

  # return the results
  rst
}


# This function computes the mean and variance of the IRT based residual statistics
# and it computes the covariance between rdifr and rdifs
resid_moments <- function(p_ref, p_foc, n_ref, n_foc, resp_ref, resp_foc, cats) {
  # check the number of items
  nitem <- length(cats)

  # check the number of rows and columns in the two groups response data
  nrow_ref <- nrow(resp_ref)
  nrow_foc <- nrow(resp_foc)

  # create an empty list of two matrices to contain for the first and second moments of raw residual and squared residuals
  mu_ref <- purrr::map(.x = 1:2, .f = function(x) matrix(NA, nrow = nrow_ref, ncol = nitem))
  mu_foc <- purrr::map(.x = 1:2, .f = function(x) matrix(NA, nrow = nrow_foc, ncol = nitem))
  mu2_ref <- mu_ref
  mu2_foc <- mu_foc

  # create an empty matrix to contain the expectation of cube of raw residuals at each theta for all items
  mu_resid3_ref <- matrix(NA, nrow = nrow_ref, ncol = nitem)
  mu_resid3_foc <- matrix(NA, nrow = nrow_foc, ncol = nitem)

  # compute the first and second moments of raw residual and squared residuals for all items
  for (i in 1:nitem) {
    # compute the expected residuals for each score category
    Emat_ref <- matrix(0:(cats[i] - 1), nrow = nrow_ref, ncol = cats[i], byrow = TRUE)
    Emat_foc <- matrix(0:(cats[i] - 1), nrow = nrow_foc, ncol = cats[i], byrow = TRUE)
    exp_resid_ref <- Emat_ref - matrix(rowSums(Emat_ref * p_ref[[i]]), nrow = nrow_ref, ncol = cats[i], byrow = FALSE)
    exp_resid_foc <- Emat_foc - matrix(rowSums(Emat_foc * p_foc[[i]]), nrow = nrow_foc, ncol = cats[i], byrow = FALSE)

    # replace NA values into the missing data location
    exp_resid_ref[is.na(resp_ref[, i]), ] <- NA
    exp_resid_foc[is.na(resp_foc[, i]), ] <- NA

    # compute the expected values of raw and squared residuals to be used to compute the mean and variance parameters
    value_foc <- value_ref <- vector("list", 2)
    value_ref[[1]] <- exp_resid_ref
    value_foc[[1]] <- exp_resid_foc
    value_ref[[2]] <- exp_resid_ref^2
    value_foc[[2]] <- exp_resid_foc^2
    resid3_ref <- exp_resid_ref^3
    resid3_foc <- exp_resid_foc^3

    # compute the first and second moments across all thetas for the reference and focal groups
    for (j in 1:2) {
      mu_ref[[j]][, i] <- Rfast::rowsums(value_ref[[j]] * p_ref[[i]], na.rm = FALSE)
      mu2_ref[[j]][, i] <- Rfast::rowsums((value_ref[[j]])^2 * p_ref[[i]], na.rm = FALSE)
      mu_foc[[j]][, i] <- Rfast::rowsums(value_foc[[j]] * p_foc[[i]], na.rm = FALSE)
      mu2_foc[[j]][, i] <- Rfast::rowsums((value_foc[[j]])^2 * p_foc[[i]], na.rm = FALSE)
    }

    # compute the expectation (first moment) of the cube of raw residuals across all thetas for the reference and focal groups
    mu_resid3_ref[, i] <- Rfast::rowsums(resid3_ref * p_ref[[i]], na.rm = FALSE)
    mu_resid3_foc[, i] <- Rfast::rowsums(resid3_foc * p_foc[[i]], na.rm = FALSE)
  }

  # compute the variances across all thetas for the reference and focal groups
  var_ref <- purrr::map2(
    .x = mu2_ref, .y = mu_ref,
    .f = function(x, y) {
      x - y^2
    }
  )
  var_foc <- purrr::map2(
    .x = mu2_foc, .y = mu_foc,
    .f = function(x, y) {
      x - y^2
    }
  )

  # compute the sum of variances for each group
  # use V(aX - bY) = a^2 * V(X) + b^2 * V(Y)
  const_ref <- purrr::map(
    .x = var_ref,
    .f = function(x) {
      matrix((1 / n_ref^2), nrow = nrow(x), ncol = ncol(x), byrow = TRUE)
    }
  )
  const_foc <- purrr::map(
    .x = var_foc,
    .f = function(x) {
      matrix((1 / n_foc^2), nrow = nrow(x), ncol = ncol(x), byrow = TRUE)
    }
  )
  var_mu_ref <- purrr::map2(
    .x = const_ref, .y = var_ref,
    .f = function(x, y) {
      Rfast::colsums(x * y, na.rm = TRUE)
    }
  )
  var_mu_foc <- purrr::map2(
    .x = const_foc, .y = var_foc,
    .f = function(x, y) {
      Rfast::colsums(x * y, na.rm = TRUE)
    }
  )

  # compute the final mean and variance
  # use E(X - Y) = E(X) - E(Y)
  # use V(X - Y) = V(X) + V(Y)
  mu <- purrr::map2(
    .x = mu_foc, .y = mu_ref,
    .f = function(x, y) {
      colMeans(x, na.rm = TRUE) - colMeans(y, na.rm = TRUE)
    }
  )
  mu[[1]] <- round(mu[[1]], digits = 10)
  sigma2 <- purrr::map2(.x = var_mu_foc, .y = var_mu_ref, .f = function(x, y) x + y)
  sigma <- purrr::map(.x = sigma2, .f = function(x) sqrt(x))

  # compute the covariance between rdifr and rdifs
  covar <- Rfast::colsums(mu_resid3_foc, na.rm = TRUE) * (1 / n_foc^2) + Rfast::colsums(mu_resid3_ref, na.rm = TRUE) * (1 / n_ref^2)

  # return the results
  rst <- purrr::pmap(
    .l = list(x = mu, y = sigma2, z = sigma),
    .f = function(x, y, z) data.frame(mu = x, sigma2 = y, sigma = z)
  )
  names(rst) <- c("rdifr", "rdifs")
  rst$covariance <- covar
  rst
}
