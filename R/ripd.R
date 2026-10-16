#' Residual-based Item Parameter Drift (RIPD) Detection Framework
#'
#' This function computes three RIPD statistics, \eqn{RIPD_{R}}, \eqn{RIPD_{S}},
#' and \eqn{RIPD_{RS}}, for each item. \eqn{RIPD_{R}} captures differences in
#' mean raw residuals between groups, which is typically indicative of uniform
#' item parameter drift (IPD). \eqn{RIPD_{S}} captures differences in mean
#' squared residuals between groups, reflecting nonuniform IPD. \eqn{RIPD_{RS}},
#' a combined chi-square-based index, is sensitive to both uniform and
#' nonuniform IPD.
#'
#' @inheritParams rdif
#' @param score A numeric vector of examinees' ability estimates (theta values)
#'   in the same row order as `data`. If `NULL`, the ability estimates are
#'   computed internally with [irtQ::est_score()] using `method`, `range`,
#'   `norm.prior`, `nquad`, `weights`, and `ncore`. When `purify = TRUE`, the
#'   ability estimates of all examinees in both groups are re-estimated
#'   internally at every purification iteration with these settings, so a
#'   supplied `score` is used only for the initial (non-purified) analysis. A
#'   missing value in `score` excludes the examinee from the analysis, also
#'   during purification, with a warning. Default is `NULL`.
#' @param group A numeric or character vector indicating the group membership of
#'   examinees. Its length must equal the number of rows in `data`. Examinees
#'   whose value equals `focal.name` form the focal group; all other examinees
#'   are pooled into the reference group (for example, the synthetic reference
#'   group described in Details). Examinees with `NA` in `group` are excluded
#'   from both groups.
#' @param item.skip A numeric vector of item positions (row numbers of `x`), or a
#'   logical vector with one value per item, identifying the items to exclude
#'   from the IPD analysis. If `NULL`, all items are analyzed. Skipped
#'   items still contribute to the ability estimates; only their RIPD
#'   statistics and moments are set to `NA`, and they are never flagged. In CAT,
#'   this is typically used to skip low-exposure items (see Details). Default is
#'   `NULL`.
#' @param alpha A numeric value specifying the significance level (\eqn{\alpha})
#'   for hypothesis testing using the RIPD statistics. Default is `0.05`.
#' @param missing A value indicating missing values in the response data set.
#'   Default is NA.
#' @param method A character string indicating the scoring method used to
#'   estimate the abilities when `score = NULL` and in every purification
#'   iteration. Available options are:
#'   - `"ML"`: Maximum likelihood estimation
#'   - `"MLF"`: Maximum likelihood estimation with fences
#'   - `"WL"`: Weighted likelihood estimation
#'   - `"MAP"`: Maximum a posteriori estimation
#'   - `"EAP"`: Expected a posteriori estimation
#'
#'   The methods `"EAP.SUM"` and `"INV.TCC"` of [irtQ::est_score()] are not
#'   supported, because they do not give an ability estimate for each response
#'   pattern. Default is `"ML"`.
#' @param purify.by A character string specifying which RIPD statistic is used
#'   to perform the purification. Available options are "ripdrs" for
#'   \eqn{RIPD_{RS}}, "ripdr" for \eqn{RIPD_{R}}, and "ripds" for
#'   \eqn{RIPD_{S}}. Used only when `purify = TRUE`. Default is `"ripdrs"`.
#' @param max.iter A positive integer specifying the maximum number of
#'   purification iterations. Default is 10. If the limit is reached while
#'   flagged items remain, a warning is issued, the items flagged at the last
#'   iteration are added to the flagged set, and the `complete` element of
#'   `with_purify` is `FALSE`. Because one item is removed per iteration,
#'   `max.iter` should be at least the number of items expected to drift (Lim &
#'   Han, 2026, used 80).
#' @param min.resp A single nonnegative number, usually a positive integer,
#'   specifying the minimum number of valid item responses required from a
#'   focal group examinee for the examinee's responses to be used. All
#'   responses of focal group examinees with fewer than `min.resp` (but at least
#'   one) responses are set to `NA` before the ability estimation, in the
#'   initial analysis and at every purification iteration, and these examinees
#'   are excluded from the analysis. A warning reports the number of examinees
#'   excluded in the initial analysis, and one warning at the end of the
#'   purification reports the examinees excluded during it. A large `min.resp`
#'   can leave the focal group without any ability estimate in a later
#'   purification iteration, which stops the purification with a warning.
#'   Reference group examinees are not affected. If `NULL`, no minimum is
#'   applied. Default is `NULL`.
#'
#' @return This function returns an object of class `"ripd"`, a list with the
#' following elements:
#'
#' \item{no_purify}{A list of sub-objects containing the results of IPD analysis
#' without applying a purification procedure. The sub-objects include:
#'   \describe{
#'     \item{ipd_stat}{A data frame of RIPD results for all items with the
#'     columns `id`, `ripdr`, `z.ripdr` (standardized \eqn{RIPD_{R}}), `ripds`,
#'     `z.ripds` (standardized \eqn{RIPD_{S}}), `ripdrs`, `p.ripdr`, `p.ripds`,
#'     `p.ripdrs`, `n.ref`, `n.foc`, and `n.total`. `n.ref` and `n.foc` are the
#'     numbers of examinees in each group who responded to the item and have an
#'     ability estimate. Statistics and p-values are rounded to four decimal
#'     places, and items are flagged using the unrounded p-values; they are `NA`
#'     for items in `item.skip`. \eqn{RIPD_{RS}} has no standardized value
#'     because it is a \eqn{\chi^{2}}-based statistic.}
#'     \item{moments}{A data frame of the null means and standard deviations of
#'     the RIPD statistics with the columns `id`, `mu.ripdr`, `sigma.ripdr`,
#'     `mu.ripds`, `sigma.ripds`, and `covariance` (the covariance between
#'     \eqn{RIPD_{R}} and \eqn{RIPD_{S}}). Values are `NA` for items in
#'     `item.skip`.}
#'     \item{ipd_item}{A list with the elements `ripdr`, `ripds`, and `ripdrs`,
#'     each giving the indices (row positions in `x`) of the items flagged by the
#'     corresponding statistic at level `alpha`. An element is `NULL` when no
#'     item is flagged.}
#'     \item{score}{A numeric vector of ability estimates used to compute the RIPD
#'     statistics.}
#'   }
#' }
#'
#' \item{purify}{A logical value indicating whether the purification procedure
#' was applied.}
#'
#' \item{with_purify}{A list of sub-objects containing the results of IPD analysis
#' with a purification procedure. All elements are `NULL` when `purify = FALSE`.
#' The sub-objects include:
#'   \describe{
#'     \item{purify.by}{A character string indicating the RIPD statistic used for
#'     purification. Possible values are "ripdr", "ripds", and "ripdrs",
#'     corresponding to \eqn{RIPD_{R}}, \eqn{RIPD_{S}}, and \eqn{RIPD_{RS}},
#'     respectively.}
#'     \item{ipd_stat}{A data frame of the RIPD results for all items with the
#'     same columns as in \code{no_purify} plus `n.iter`, the analysis from
#'     which each row is taken. For a flagged item, `n.iter` is the iteration
#'     at which it was flagged (0 is the initial analysis); for the other items,
#'     it is the final iteration.}
#'     \item{moments}{A data frame of the null moments for all items with the
#'     same columns as in \code{no_purify} plus `n.iter`, taken from the same
#'     analysis as the corresponding row of `ipd_stat`.}
#'     \item{ipd_item}{A numeric vector of item indices (row positions in
#'     \code{x}) identified as IPD items across all purification iterations,
#'     sorted in ascending order. If \code{max.iter} is reached before
#'     convergence or the purification stops because a group has no ability
#'     estimate left, the items flagged at the last iteration (iteration
#'     \code{n.iter}) are also included. It contains all items if every item is
#'     flagged. \code{NULL} if no item is flagged in the initial analysis.}
#'     \item{n.iter}{An integer indicating the total number of iterations
#'     performed during the purification process. The iteration that removes the
#'     last item, or after which a group has no ability estimate left, is not
#'     counted, because the statistics are not recomputed.}
#'     \item{score}{A numeric vector of purified ability estimates used to
#'     compute the final RIPD statistics, from the last iteration counted in
#'     \code{n.iter}. \code{NULL} when no item is flagged in the first analysis,
#'     because purification is not carried out.}
#'     \item{complete}{A logical value indicating whether the purification
#'     process converged. \code{FALSE} if the procedure stopped because the
#'     maximum number of iterations was reached, because every item was flagged,
#'     or because a group had no examinee with an ability estimate left.}
#'   }
#' }
#'
#' \item{alpha}{A numeric value indicating the significance level (\eqn{\alpha})
#' used in hypothesis testing for RIPD statistics.}
#'
#' \item{call}{The matched function call.}
#'
#' @details
#' The current version supports dichotomous items only. The function stops if
#' `x` contains GRM or GPCM items or if `data` contain response values other
#' than 0 and 1.
#'
#' \strong{Theoretical Background: From RDIF to RIPD}
#'
#' The RIPD framework directly adapts the residual-based differential item
#' functioning (RDIF) detection framework (Lim et al., 2022; Lim & Choe, 2023)
#' to detect item parameter drift (IPD) in computerized adaptive testing (CAT).
#' Each RIPD statistic (\eqn{RIPD_R}, \eqn{RIPD_S}, \eqn{RIPD_{RS}}) mirrors
#' its RDIF counterpart (\eqn{RDIF_R}, \eqn{RDIF_S}, \eqn{RDIF_{RS}}) in both
#' computation and asymptotic theory, with two key adaptations. First, whereas
#' the original RDIF framework estimates item parameters once from the pooled
#' data of both groups, the RIPD framework uses the \emph{original,
#' pre-calibrated} item parameters from the CAT pool directly, eliminating the
#' need for item recalibration. Second, the reference group is a synthetic,
#' drift-free group simulated from the focal group's ability estimates (see
#' CAT-Specific Workflow). This makes RIPD practically scalable for operational
#' CAT programs, where the response data are sparse and recalibration is often
#' impractical.
#'
#' IPD in CAT can be viewed as a special case of DIF across time (Lord, 1980;
#' Veerkamp & Glas, 2000): the focal group represents current test-takers whose
#' responses may reflect drift, while the reference group represents a baseline
#' under invariant item parameters.
#'
#' \strong{Residual Definition}
#'
#' For examinee \eqn{i} and item \eqn{j}, the raw residual is defined as
#' \deqn{e_{ij} = u_{ij} - P_j(\hat{\theta}_i),}
#' where \eqn{u_{ij}} is the observed binary response (0 or 1) and
#' \eqn{P_j(\hat{\theta}_i)} is the model-predicted probability of a correct
#' response under the \emph{original} item parameters evaluated at the
#' examinee's ability estimate \eqn{\hat{\theta}_i}. Because item responses
#' are independent but not identically distributed (due to varying predicted
#' probabilities across examinees with different ability levels), the asymptotic
#' distributions of the RIPD statistics are established via Lyapunov's central
#' limit theorem. Below, \eqn{P_{ij} = P_j(\hat{\theta}_i)}, and \eqn{N_F} and
#' \eqn{N_R} denote the numbers of focal and reference group examinees who
#' responded to item \eqn{j}; all sums run over these examinees only.
#'
#' \strong{Three RIPD Statistics and Their Asymptotic Distributions}
#'
#' \emph{RIPD\eqn{_R}} (uniform drift): Measures the difference in mean raw
#' residuals between the focal (\eqn{F}) and reference (\eqn{R}) groups,
#' \deqn{RIPD_{R,j} = \frac{\sum_{i \in F} e_{ij}}{N_F} -
#'   \frac{\sum_{i \in R} e_{ij}}{N_R}.}
#' \eqn{RIPD_{R,j}} is most effective when the item response functions (IRFs)
#' of the two groups differ primarily in location (e.g., drift in the difficulty
#' \eqn{b} or guessing \eqn{c} parameter). Under the null hypothesis \eqn{H_0}
#' of no IPD, \eqn{RIPD_{R,j}} asymptotically follows a normal distribution
#' with mean zero and analytically derived variance
#' \eqn{\sigma^2_{RIPD_R} = \frac{\sum_{i \in F} P_{ij}(1 - P_{ij})}{N_F^2} +
#'   \frac{\sum_{i \in R} P_{ij}(1 - P_{ij})}{N_R^2}}.
#' A standard Z-test is then applied as \eqn{Z_R = RIPD_{R,j} / \sigma_{RIPD_R}}.
#'
#' \emph{RIPD\eqn{_S}} (nonuniform drift): Measures the difference in mean
#' squared residuals between groups,
#' \deqn{RIPD_{S,j} = \frac{\sum_{i \in F} e_{ij}^2}{N_F} -
#'   \frac{\sum_{i \in R} e_{ij}^2}{N_R}.}
#' \eqn{RIPD_{S,j}} captures drift that varies across the ability continuum,
#' such as a change in the discrimination parameter \eqn{a}, which
#' \eqn{RIPD_{R,j}} alone may not detect. Importantly, under \eqn{H_0}, the
#' null mean of \eqn{RIPD_{S,j}} is \strong{not zero} in general; it equals
#' \deqn{\mu_{RIPD_S} = \frac{\sum_{i \in F} P_{ij}(1 - P_{ij})}{N_F} -
#'   \frac{\sum_{i \in R} P_{ij}(1 - P_{ij})}{N_R},}
#' which depends on the predicted probabilities of the examinees in each group
#' who responded to item \eqn{j}. Even with a synthetic reference group, the
#' ability estimates and the administered items of the reference group differ
#' from those of the focal group, so this null mean is generally nonzero and is
#' computed item by item. The asymptotic variance is also analytically
#' derived. A standard Z-test is applied as
#' \eqn{Z_S = (RIPD_{S,j} - \mu_{RIPD_S}) / \sigma_{RIPD_S}}.
#'
#' \emph{RIPD\eqn{_{RS}}} (combined): A Wald-type statistic that jointly tests
#' for uniform and nonuniform drift within a single framework,
#' \deqn{RIPD_{RS,j} = (\boldsymbol{\nu}_j - \boldsymbol{\mu}_j)^\top
#'   \hat{\boldsymbol{\Sigma}}_j^{-1}
#'   (\boldsymbol{\nu}_j - \boldsymbol{\mu}_j),}
#' where \eqn{\boldsymbol{\nu}_j = (RIPD_{R,j},\ RIPD_{S,j})^\top},
#' \eqn{\boldsymbol{\mu}_j = (\mu_{RIPD_R},\ \mu_{RIPD_S})^\top} is the
#' vector of null means, and \eqn{\hat{\boldsymbol{\Sigma}}_j} is the
#' analytically derived \eqn{2 \times 2} covariance matrix of
#' \eqn{(RIPD_{R,j}, RIPD_{S,j})}. Under \eqn{H_0}, \eqn{RIPD_{RS,j}}
#' asymptotically follows a \eqn{\chi^2} distribution with 2 degrees of
#' freedom. \eqn{RIPD_{RS}} is the \strong{recommended} statistic in practice
#' because it is sensitive to both uniform and nonuniform drift and addresses
#' potential inflation of the family-wise Type I error that arises from
#' applying \eqn{RIPD_R} and \eqn{RIPD_S} separately to the same item.
#'
#' In rare cases, the covariance matrix of \eqn{RIPD_R} and \eqn{RIPD_S} is
#' singular or nearly singular (its reciprocal condition number is 1e-10 or
#' less), for example when all examinees who responded to an item have the same
#' probability of a correct response, so that \eqn{RIPD_S} is a linear
#' function of \eqn{RIPD_R}. Then \eqn{RIPD_{RS}} is computed with the
#' Moore-Penrose generalized inverse of the covariance matrix and compared with
#' a chi-square distribution whose degrees of freedom equal the rank of the
#' matrix (Moore, 1977), and a warning names the item. This handling is not part
#' of the original method (Lim et al., 2022; Lim & Han, 2026) and is added in
#' irtQ. In most data the matrix is not singular, and \eqn{RIPD_{RS}} has two
#' degrees of freedom.
#'
#' \strong{Diagnosing the Nature of Drift}
#'
#' The pattern of flagging across the three statistics can help diagnose the
#' type of drift. If an item is flagged by \eqn{RIPD_R} (and possibly
#' \eqn{RIPD_{RS}}) but not \eqn{RIPD_S}, uniform drift is indicated (e.g.,
#' difficulty shift). If flagged by \eqn{RIPD_S} (and possibly
#' \eqn{RIPD_{RS}}) but not \eqn{RIPD_R}, nonuniform drift is suggested (e.g.,
#' discrimination change). Items flagged by all statistics or by
#' \eqn{RIPD_{RS}} alone likely exhibit mixed drift.
#'
#' \strong{CAT-Specific Workflow}
#'
#' The RIPD procedure for CAT consists of three steps (Lim & Han, 2026):
#' \enumerate{
#'   \item \strong{Focal group CAT}: The current cohort of examinees takes the
#'     CAT using the operational (potentially drifted) item pool. Each examinee
#'     receives only a subset of items; the resulting response matrix is sparse.
#'     The final ability estimates \eqn{\hat{\theta}} (e.g., ML estimates) from
#'     this step are retained.
#'   \item \strong{Synthetic reference group construction}: The focal group's
#'     \eqn{\hat{\theta}} values are treated as true abilities. A CAT
#'     simulation with the same item pool and algorithm settings as in step 1 is
#'     run, generating item responses from the \emph{original} (pre-calibrated,
#'     drift-free) item parameters. This creates a synthetic reference group
#'     that matches the focal group's ability distribution, so that any
#'     systematic differences in residuals can be attributed to IPD rather than
#'     ability confounds. The reference group size is controlled by a
#'     replication factor: e.g., 1F reuses the focal \eqn{\hat{\theta}} values
#'     once (N_ref = N_foc), while kF replicates them k times. Larger reference
#'     groups reduce sampling variability in the reference residuals and improve
#'     detection power. In the simulation study of Lim and Han (2026), power
#'     gains diminished beyond 5F, and a reference group of at least 8 to 10
#'     times the focal group size (8F to 10F) was suggested as a practical
#'     guideline. Note that \code{ripd()} does not simulate the reference group;
#'     this step must be carried out with a CAT simulation program outside
#'     irtQ, and the simulated responses and ability estimates are then
#'     combined with the focal group data.
#'   \item \strong{RIPD computation with purification}: Residuals are computed
#'     for both groups using the original item parameters, and the three RIPD
#'     statistics are evaluated item-by-item. Items whose Z-test (for
#'     \eqn{RIPD_R} or \eqn{RIPD_S}) or chi-square test (for
#'     \eqn{RIPD_{RS}}) exceeds the critical value at level \code{alpha} are
#'     flagged as drifting.
#' }
#'
#' \strong{Purification Procedure}
#'
#' When \code{purify = TRUE}, an iterative purification procedure adapted from
#' Lim et al. (2022) is applied to mitigate the bias in ability estimates caused
#' by drifted items. This is analogous to the contaminating effect of DIF items
#' on matching variables. At each iteration:
#' \enumerate{
#'   \item Among the remaining items flagged by the \code{purify.by} statistic,
#'     the item with the smallest p-value (compared on the log scale) is
#'     removed and added to the flagged set. For \eqn{RIPD_R} and
#'     \eqn{RIPD_S}, this is the item with the largest absolute z value.
#'   \item Ability estimates of all examinees in both groups are recomputed with
#'     \code{\link{est_score}}, excluding all currently flagged items and using
#'     the \code{method}, \code{range}, \code{norm.prior}, \code{nquad},
#'     \code{weights}, and \code{ncore} arguments.
#'   \item RIPD statistics are recalculated using the updated ability estimates
#'     for the remaining items.
#' }
#' The process continues until no additional items are flagged (convergence) or
#' \code{max.iter} is reached. In the latter case, a warning is issued, the
#' items flagged at the last iteration are added to the flagged set, and
#' \code{complete} is set to \code{FALSE}. If every item is flagged, the
#' procedure stops with a warning, because no item is left to estimate the
#' abilities, all items are returned as flagged, and \code{complete} is
#' \code{FALSE}. The procedure also stops with a warning and \code{complete =
#' FALSE} when, after an item is removed, no examinee of the reference group or
#' of the focal group has an ability estimate left (see \code{min.resp}); the
#' remaining items then keep the results of the previous iteration, and the
#' items flagged in that iteration are added to the flagged items.
#' Examinees without any response have no ability estimate and are excluded.
#' During purification, an examinee whose remaining responses are all removed
#' with the flagged items, or who has fewer than \code{min.resp} responses to
#' the remaining items, is also excluded, and one warning at the end reports the
#' number of such examinees.
#' The statistic used to drive purification is specified by \code{purify.by};
#' \code{"ripdrs"} (\eqn{RIPD_{RS}}) is recommended as it is sensitive to both
#' types of drift.
#'
#' \strong{Key Items and item.skip}
#'
#' In CAT contexts, many items in the pool receive few responses and thus lack
#' sufficient data for reliable IPD analysis. The \code{item.skip} argument
#' excludes such items (e.g., low-exposure non-key items) from the IPD tests.
#' Skipped items are still used to compute the ability estimates. In practice,
#' key items for RIPD evaluation are typically identified through a preliminary
#' CAT simulation as those with the highest average exposure frequencies under
#' the intended operational settings.
#'
#' \strong{Printing}
#'
#' The \code{print()} method of a \code{"ripd"} object shows the flagged items
#' and the RIPD statistics of the analyses without and with purification. Its
#' argument \code{what} selects \code{"no_purify"} or \code{"with_purify"}
#' (default \code{"all"}). Items without statistics (skipped items and items
#' without responses in a group) are not shown in the tables.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @examples
#' ## --- RIPD Example: Detecting IPD in CAT ---------------------------------
#' ##
#' ## Background (Lim & Han, 2026):
#' ##   In CAT-based IPD detection using RIPD, the reference group is
#' ##   "synthetic": it is created by re-administering a CAT to examinees whose
#' ##   true abilities are set equal to the focal group's ML theta estimates,
#' ##   using the ORIGINAL (pre-drift) item parameters. This eliminates the
#' ##   need for recalibration and makes RIPD directly applicable to
#' ##   operational CAT settings.
#' ##
#' ## The simIPD dataset contains:
#' ##   - foc_resp / foc_score : focal group CAT responses + ML theta estimates
#' ##                            (IPD items: a and b each drifted by -0.5)
#' ##   - ref_resp / ref_score : synthetic reference group CAT responses + ML
#' ##                            theta estimates (original parameters, 1F size)
#' ##   - item_par  : original (non-drifted) 360-item pool in irtQ format
#' ##   - key_item  : indices of 90 key items (highly exposed in CAT)
#' ##   - item.skip : indices of 270 non-key items (excluded from RIPD)
#' ##   - ipd_item  : indices of 18 truly drifted items (ground truth)
#' ## -------------------------------------------------------------------------
#'
#' data(simIPD)
#'
#' ## Step 1. Combine focal and synthetic reference group data
#' ##         (the row order does not matter as long as group matches the rows)
#' data  <- rbind(simIPD$ref_resp,  simIPD$foc_resp)
#' score <- c(simIPD$ref_score,     simIPD$foc_score)
#' group <- c(rep(0, nrow(simIPD$ref_resp)),   # 0 = reference
#'            rep(1, nrow(simIPD$foc_resp)))   # 1 = focal
#'
#' ## Step 2. Run RIPD without purification
#' ##         item.skip excludes the 270 non-key items from analysis
#' ripd_result <- ripd(
#'   x          = simIPD$item_par,
#'   data       = data,
#'   score      = score,
#'   group      = group,
#'   focal.name = 1,           # focal group is coded as 1
#'   item.skip  = simIPD$item.skip,
#'   D          = 1.7,
#'   alpha      = 0.05,
#'   purify     = FALSE,
#'   method     = "ML",
#'   range      = c(-5, 5)
#' )
#'
#' ## Step 3. Compare the items flagged by RIPD_RS to the ground truth
#' detected <- ripd_result$no_purify$ipd_item$ripdrs
#' cat("Truly drifted items (ground truth):", simIPD$ipd_item, "\n")
#' cat("RIPD_RS-flagged items:             ", detected, "\n")
#' cat("True positives:", sum(detected %in% simIPD$ipd_item), "of",
#'     length(simIPD$ipd_item), "\n")
#'
#' ## Step 4. Purification (not run: it takes several minutes because the
#' ##         abilities of the examinees are re-estimated at every iteration)
#' ##
#' ##   ripd_purified <- ripd(
#' ##     x = simIPD$item_par, data = data, score = score, group = group,
#' ##     focal.name = 1, item.skip = simIPD$item.skip, D = 1.7, alpha = 0.05,
#' ##     purify = TRUE, purify.by = "ripdrs", max.iter = 30,
#' ##     method = "ML", range = c(-5, 5)
#' ##   )
#' ##   print(ripd_purified, what = "with_purify")
#'
#' ## -- Note on reference group size -----------------------------------------
#' ## This example uses a 1F reference group (n_ref = n_foc = 3,000).
#' ## In practice, a larger reference group improves detection power;
#' ## Lim and Han (2026) suggest at least 8F to 10F. To create, for example,
#' ## an 8F reference group, replicate the focal ability estimates and run
#' ## the same CAT simulation with the original item parameters using a CAT
#' ## simulation program (irtQ does not provide one):
#' ##
#' ##   theta_8F <- rep(simIPD$foc_score, times = 8)   # 24,000 examinees
#' ##   # simulate a CAT for theta_8F with simIPD$item_par, then call ripd()
#' ##   # with the simulated responses and final ability estimates
#' ## -------------------------------------------------------------------------
#'
#' @seealso [irtQ::rdif()], [irtQ::est_irt()], [irtQ::est_item()],
#'   [irtQ::simdat()], [irtQ::shape_df()], [irtQ::est_score()],
#'   [irtQ::pcd2()], [irtQ::simIPD]
#'
#' @references Lim, H., & Choe, E. M. (2023). Detecting differential item
#'   functioning in CAT using IRT residual DIF approach. *Journal of Educational
#'   Measurement, 60*(4), 626-650. \doi{10.1111/jedm.12366}.
#'
#'   Lim, H., Choe, E. M., & Han, K. T. (2022). A residual-based differential
#'   item functioning detection framework in item response theory. *Journal of
#'   Educational Measurement, 59*(1), 80-104. \doi{10.1111/jedm.12313}.
#'
#'   Lim, H., & Han, K. T. (2026). IRT residual-based approach to detecting item
#'   parameter drift in CAT. *Journal of Educational and Behavioral Statistics*.
#'   \doi{10.3102/10769986261460852}.
#'
#'   Lord, F. M. (1980). Applications of item response theory to practical
#'   testing problems. Lawrence Erlbaum Associates.
#'
#'   Moore, D. S. (1977). Generalized inverses, Wald's method, and the
#'   construction of chi-squared tests of fit. *Journal of the American
#'   Statistical Association, 72*(357), 131-137.
#'   \doi{10.1080/01621459.1977.10479921}.
#'
#'   Veerkamp, W. J. J., & Glas, C. A. W. (2000). Detection of known items in
#'   adaptive testing with a statistical quality control method. *Journal of
#'   Educational and Behavioral Statistics, 25*(4), 373-389.
#'   \doi{10.3102/10769986025004373}.
#'
#' @order 1
#'@export
ripd <- function(x, ...) UseMethod("ripd")

#' @describeIn ripd Default method for a data frame `x` containing item metadata.
#'
#' @order 2
#'@export
ripd.default <- function(x,
                         data,
                         score = NULL,
                         group,
                         focal.name,
                         item.skip = NULL,
                         D = 1,
                         alpha = 0.05,
                         missing = NA,
                         purify = FALSE,
                         purify.by = c("ripdrs", "ripdr", "ripds"),
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

  # conduct the IPD analysis
  rst <- ripd_main(
    x = x, data = data, score = score, group = group, focal.name = focal.name,
    item.skip = item.skip, D = D, alpha = alpha, missing = missing, purify = purify,
    purify.by = purify.by, max.iter = max.iter, min.resp = min.resp, method = method,
    range = range, norm.prior = norm.prior, nquad = nquad, weights = weights,
    ncore = ncore, verbose = verbose, ...
  )

  # return the IPD detection results
  rst$call <- cl
  rst
}

#' @describeIn ripd Method for an object of class `est_irt` created by
#' [irtQ::est_irt()]. The response data and the scaling factor `D` are taken from
#' the object.
#'
#' @order 3
#'@export
ripd.est_irt <- function(x,
                         score = NULL,
                         group,
                         focal.name,
                         item.skip = NULL,
                         alpha = 0.05,
                         missing = NA,
                         purify = FALSE,
                         purify.by = c("ripdrs", "ripdr", "ripds"),
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

  # conduct the IPD analysis with the data, scaling factor, and items of the object
  rst <- ripd_main(
    x = x$par.est, data = x$data, score = score, group = group, focal.name = focal.name,
    item.skip = item.skip, D = x$scale.D, alpha = alpha, missing = missing,
    purify = purify, purify.by = purify.by, max.iter = max.iter, min.resp = min.resp,
    method = method, range = range, norm.prior = norm.prior, nquad = nquad,
    weights = weights, ncore = ncore, verbose = verbose, ...
  )

  # return the IPD detection results
  rst$call <- cl
  rst
}

#' @describeIn ripd Method for an object of class `est_item` created by
#' [irtQ::est_item()]. The response data, ability estimates, and the scaling
#' factor `D` are taken from the object.
#'
#' @order 4
#'@export
ripd.est_item <- function(x,
                          group,
                          focal.name,
                          item.skip = NULL,
                          alpha = 0.05,
                          missing = NA,
                          purify = FALSE,
                          purify.by = c("ripdrs", "ripdr", "ripds"),
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

  # conduct the IPD analysis with the data, ability estimates, scaling factor, and items of the object
  rst <- ripd_main(
    x = x$par.est, data = x$data, score = x$score, group = group, focal.name = focal.name,
    item.skip = item.skip, D = x$scale.D, alpha = alpha, missing = missing,
    purify = purify, purify.by = purify.by, max.iter = max.iter, min.resp = min.resp,
    method = method, range = range, norm.prior = norm.prior, nquad = nquad,
    weights = weights, ncore = ncore, verbose = verbose, ...
  )

  # return the IPD detection results
  rst$call <- cl
  rst
}


# This function conducts the IPD analysis shared by the methods of ripd(): it prepares the
# data and ability estimates, runs the analysis once, and applies the purification procedure
ripd_main <- function(x, data, score, group, focal.name, item.skip, D, alpha, missing,
                      purify, purify.by, max.iter, min.resp, method, range, norm.prior,
                      nquad, weights, ncore, verbose, ...) {

  ## ----------------------------------
  ## (1) prepare IPD analysis
  ## ----------------------------------
  # confirm and correct all item metadata information
  x <- confirm_df(x)

  # stop when the model includes any polytomous model
  if (any(x$model %in% c("GRM", "GPCM")) | any(x$cats > 2)) {
    stop("The current version only supports dichotomous response data.", call. = FALSE)
  }

  # transform the response data to a matrix form
  data <- as.matrix(data)

  # re-code missing values
  if (!is.na(missing)) {
    data[data == missing] <- NA
  }

  # transform the response data to a numeric matrix and check the responses
  data <- resp_to_matrix(data, x$cats, x$id)

  # check the group vector and the focal group
  check_group(group, focal.name, nrow(data))

  # check the positions of the items to be skipped
  item.skip <- check_item_skip(item.skip, nrow(x))

  # stop when the significance level is not a single value between 0 and 1
  check_alpha(alpha)

  # stop when the maximum number of iterations is not a positive whole number
  check_max_iter(max.iter, purify)

  # stop when the minimum number of responses is not a single nonnegative number
  check_min_resp(min.resp)

  # stop when the abilities are estimated or updated with an unsupported scoring method
  if (is.null(score) || purify) {
    check_pattern_method(method, fn = "ripd()")
  }

  # record whether the ability estimates are supplied
  score_given <- !is.null(score)

  # find the focal group examinees with fewer than min.resp responses
  loc_less <- NULL
  if (!is.null(min.resp)) {
    loc_less <- ripd_min_resp(data = data, group = group, focal.name = focal.name,
                              min.resp = min.resp)
  }

  # transform the scores to a vector form when they are provided
  if (!is.null(score)) {
    # transform the scores to a vector form and check them
    score <- check_score(score, nrow(data))

    # warn that the examinees with responses but without an ability estimate are excluded
    n_na <- sum(is.na(score) & rowSums(!is.na(data)) > 0)
    if (n_na > 0L) {
      warning(n_na, " examinee(s) with a missing ability estimate were excluded.", call. = FALSE)
    }

    # set all responses of the examinees without an ability estimate to NA, so that they stay
    # excluded during the purification
    data[is.na(score), ] <- NA
  }

  # warn that the focal group examinees with too few responses are excluded
  if (length(loc_less) > 0L) {
    warning(length(loc_less), " focal group examinee(s) with fewer than ", min.resp,
            " responses were excluded.", call. = FALSE)
  }

  # set all responses of those examinees to NA
  if (length(loc_less) > 0L) {
    data[loc_less, ] <- NA
  }

  # compute the score if score = NULL
  if (is.null(score)) {
    # estimate the abilities; the examinees without a response are excluded without a warning
    score <- dif_rescore(
      x = x, data = data, score = NULL, loc_resp = NULL, first = TRUE, D = D, method = method,
      range = range, norm.prior = norm.prior, nquad = nquad, weights = weights, ncore = ncore,
      quiet_na = TRUE, ...
    )
  } else if (length(loc_less) > 0L) {
    # treat the ability estimates of those examinees as missing
    score[loc_less] <- NA
  }

  # stop when no examinee has an ability estimate
  if (all(is.na(score))) {
    stop("No examinee has an ability estimate.", call. = FALSE)
  }

  # a) when no purification is set
  # do only one iteration of IPD analysis
  ipd_rst <- ripd_one(
    x = x, data = data, score = score, group = group, focal.name = focal.name,
    item.skip = item.skip, D = D, alpha = alpha
  )

  # record the items whose covariance matrix is singular
  singular_id <- ipd_rst$singular

  # create two empty lists to contain the results
  no_purify <- list(ipd_stat = NULL, moments = NULL, ipd_item = NULL, score = NULL)
  with_purify <- list(
    purify.by = NULL, ipd_stat = NULL, moments = NULL,
    ipd_item = NULL, n.iter = NULL, score = NULL, complete = NULL
  )

  # record the first IPD detection results into the no purification list
  no_purify$ipd_stat <- ipd_rst$ipd_stat
  no_purify$ipd_item <- ipd_rst$ipd_item
  no_purify$moments <- ripd_moments(ipd_rst)
  no_purify$score <- score

  # when purification is used
  if (purify) {
    # verify the criterion for purification
    purify.by <- match.arg(purify.by, c("ripdrs", "ripdr", "ripds"))

    # create an empty vector and empty data frames
    # to contain the detected IPD items, statistics, and moments
    ipd_item <- NULL
    ipd_stat <-
      data.frame(
        id = rep(NA_character_, nrow(x)), ripdr = NA, z.ripdr = NA,
        ripds = NA, z.ripds = NA, ripdrs = NA, p.ripdr = NA, p.ripds = NA, p.ripdrs = NA,
        n.ref = NA, n.foc = NA, n.total = NA, n.iter = NA, stringsAsFactors = FALSE
      )
    mmt_df <-
      data.frame(
        id = rep(NA_character_, nrow(x)), mu.ripdr = NA, sigma.ripdr = NA,
        mu.ripds = NA, sigma.ripds = NA, covariance = NA, n.iter = NA,
        stringsAsFactors = FALSE
      )

    # extract the first IPD analysis results
    # and check if at least one IPD item is detected
    ipd_item_tmp <- ipd_rst$ipd_item[[purify.by]]
    log_p_tmp <- ipd_rst$log_p[[purify.by]]
    ipd_stat_tmp <- ipd_rst$ipd_stat
    mmt_df_tmp <- no_purify$moments

    # copy the response data and item meta data
    x_puri <- x
    data_puri <- data

    # start the iteration if any item is detected as an IPD item
    if (!is.null(ipd_item_tmp)) {
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
        flag_max <- pick_flagged_item(flag_loc = ipd_item_tmp, log_p = log_p_tmp)

        # check an item that is deleted
        del_item <- item_num[flag_max]

        # add the deleted item as the IPD item
        ipd_item <- c(ipd_item, del_item)

        # add the IPD statistics and moments for the detected IPD item
        ipd_stat[del_item, 1:12] <- ipd_stat_tmp[flag_max, ]
        ipd_stat[del_item, 13] <- i - 1
        mmt_df[del_item, 1:6] <- mmt_df_tmp[flag_max, ]
        mmt_df[del_item, 7] <- i - 1

        # find the examinees who responded to the item to be deleted
        loc_resp <- which(!is.na(data_puri[, flag_max]))

        # refine the leftover items
        item_num <- item_num[-flag_max]

        # stop the purification when no item is left to estimate the abilities
        if (length(item_num) == 0L) {
          all_flagged <- TRUE
          break
        }

        # remove the detected IPD item data which has the smallest p-value from the item metadata
        x_puri <- x_puri[-flag_max, , drop = FALSE]

        # remove the detected IPD item data which has the smallest p-value from the response data
        data_puri <- data_puri[, -flag_max, drop = FALSE]

        # update the locations of the items that should be skipped in the purified data
        if (!is.null(item.skip)) {
          item.skip.puri <- c(1:length(item_num))[item_num %in% item.skip]
        } else {
          item.skip.puri <- NULL
        }

        # set all responses of focal group examinees with fewer than min.resp responses to NA
        if (!is.null(min.resp)) {
          loc_less <- ripd_min_resp(data = data_puri, group = group, focal.name = focal.name,
                                    min.resp = min.resp)
          data_puri[loc_less, ] <- NA
        }

        # keep the ability estimates of the previous iteration in case the purification stops
        score_prev <- score

        # re-estimate the abilities without the flagged items; only the examinees who responded
        # to the deleted item are rescored after the first iteration
        score <-
          dif_rescore(
            x = x_puri, data = data_puri, score = score, loc_resp = loc_resp,
            first = (i == 1L), D = D, method = method, range = range,
            norm.prior = norm.prior, nquad = nquad, weights = weights, ncore = ncore,
            quiet_na = TRUE, quiet_ncore = (i > 1L || !score_given), ...
          )

        # stop the purification when the examinees of a group have all lost their ability estimates;
        # the remaining items keep the results of the previous iteration
        if (all(is.na(score[which(group == focal.name)])) ||
            all(is.na(score[which(group != focal.name)]))) {
          no_score <- TRUE
          score <- score_prev
          ipd_stat[item_num, 1:12] <- ipd_stat_tmp[-flag_max, ]
          ipd_stat[item_num, 13] <- i - 1
          mmt_df[item_num, 1:6] <- mmt_df_tmp[-flag_max, ]
          mmt_df[item_num, 7] <- i - 1

          # add the items flagged at the previous iteration, except the item deleted at this one
          ipd_item <- c(ipd_item, item_num[match(setdiff(ipd_item_tmp, flag_max),
                                                 seq_along(log_p_tmp)[-flag_max])])
          break
        }

        # do IPD analysis using the updated ability estimates
        ipd_rst_tmp <- ripd_one(
          x = x_puri, data = data_puri, score = score, group = group,
          focal.name = focal.name, item.skip = item.skip.puri, D = D,
          alpha = alpha
        )

        # record the items whose covariance matrix is singular
        singular_id <- union(singular_id, ipd_rst_tmp$singular)

        # extract the IPD analysis results
        # and check if at least one IPD item is detected
        ipd_item_tmp <- ipd_rst_tmp$ipd_item[[purify.by]]
        log_p_tmp <- ipd_rst_tmp$log_p[[purify.by]]
        ipd_stat_tmp <- ipd_rst_tmp$ipd_stat
        mmt_df_tmp <- ripd_moments(ipd_rst_tmp)

        # check if a further IPD item is flagged
        if (is.null(ipd_item_tmp)) {
          # add the IPD statistics for rest of items
          ipd_stat[item_num, 1:12] <- ipd_stat_tmp
          ipd_stat[item_num, 13] <- i
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
      # iteration that removes the last item or finds no ability estimate in a group
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
      } else if (max.iter == n_iter & !is.null(ipd_item_tmp)) {
        # if the maximum number of iterations is reached before purification is complete,
        # then, return a warning message
        warning("The maximum number of iterations was reached before purification was completed.",
                call. = FALSE)
        complete <- FALSE

        # add flagged IPD item at the last iteration
        ipd_item <- c(ipd_item, item_num[ipd_item_tmp])

        # add the IPD statistics for rest of items
        ipd_stat[item_num, 1:12] <- ipd_stat_tmp
        ipd_stat[item_num, 13] <- i
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
      warn_purify_excluded(sum(is.na(score) & rowSums(!is.na(data)) > 0))

      # record the final IPD detection results with the purification procedure
      with_purify$purify.by <- purify.by
      with_purify$ipd_stat <- ipd_stat
      with_purify$moments <- mmt_df
      with_purify$ipd_item <- sort(ipd_item)
      with_purify$n.iter <- n_iter
      with_purify$score <- score
      with_purify$complete <- complete
    } else {
      # in case when no IPD item is detected from the first IPD analysis results
      with_purify$purify.by <- purify.by
      with_purify$ipd_stat <- cbind(no_purify$ipd_stat, n.iter = 0)
      with_purify$moments <- cbind(no_purify$moments, n.iter = 0)
      with_purify$n.iter <- 0
      with_purify$complete <- TRUE
    }
  }

  # warn that the generalized inverse is used for the items with a singular covariance matrix
  if (length(singular_id) > 0L) {
    warning("The covariance matrix of RIPD_R and RIPD_S is singular or nearly singular for item(s) ",
            paste(singular_id, collapse = ", "), ". RIPD_RS of these items is computed ",
            "with a generalized inverse and fewer degrees of freedom.", call. = FALSE)
  }

  # summarize the results
  rst <- list(no_purify = no_purify, purify = purify, with_purify = with_purify, alpha = alpha)

  # return the IPD detection results
  class(rst) <- "ripd"
  rst
}


# This function finds the focal group examinees with fewer than min.resp (but at least one) responses
ripd_min_resp <- function(data, group, focal.name, min.resp) {
  # count the responses of the focal group examinees
  loc_foc <- which(group == focal.name)
  n_resp <- rowSums(!is.na(data[loc_foc, , drop = FALSE]))

  # return the row numbers of the examinees whose number of responses is less than min.resp
  loc_foc[which(n_resp < min.resp & n_resp > 0)]
}


# This function creates the data frame of the null moments of the two statistics
ripd_moments <- function(ipd_rst) {
  # combine the means, standard deviations, and covariance
  mmt_df <- data.frame(
    id = ipd_rst$ipd_stat$id,
    ipd_rst$moments$ripdr[, c(1, 3)],
    ipd_rst$moments$ripds[, c(1, 3)],
    ipd_rst$covariance, stringsAsFactors = FALSE
  )
  names(mmt_df) <- c("id", "mu.ripdr", "sigma.ripdr", "mu.ripds", "sigma.ripds", "covariance")
  mmt_df
}


# This function conducts one iteration of IPD analysis using the IRT residual-based statistics
ripd_one <- function(x, data, score, group, focal.name, item.skip = NULL, D = 1, alpha = 0.05) {
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

  # compute the model-predicted probabilities of answering correctly (a.k.a. model-expected
  # item scores) and the model probabilities of score categories
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

  # compute the residual-based IPD statistic
  ripdr <- colMeans(resid_foc, na.rm = TRUE) - colMeans(resid_ref, na.rm = TRUE) # the difference of mean raw residuals (DMRR)
  ripds <- colMeans(resid_foc^2, na.rm = TRUE) - colMeans(resid_ref^2, na.rm = TRUE) # the difference of mean squared residuals (DMSR)

  # compute the means and variances of the two statistics for the hypothesis testing
  moments <- resid_moments(
    p_ref = prob_ref, p_foc = prob_foc, n_ref = n_ref, n_foc = n_foc,
    resp_ref = resp_ref, resp_foc = resp_foc, cats = cats
  )
  moments_ripdr <- moments$rdifr
  moments_ripds <- moments$rdifs
  covar <- moments$covariance

  # compute the chi-square statistics and their degrees of freedom
  chisq <- rep(NA_real_, nitem)
  df_chisq <- rep(2, nitem)
  for (i in 1:nitem) {
    if (i %in% all_miss) {
      chisq[i] <- NaN
    } else {
      # create the variance-covariance matrix of ripdr and ripds
      cov_mat <- array(NA, c(2, 2))

      # replace NAs with the analytically computed covariance
      cov_mat[col(cov_mat) != row(cov_mat)] <- covar[i]

      # replace NAs with the analytically computed variances
      diag(cov_mat) <- c(moments_ripdr[i, 2], moments_ripds[i, 2])

      # create a vector of mean ripdr and mean ripds
      mu_vec <- cbind(moments_ripdr[i, 1], moments_ripds[i, 1])

      # create a vector of ripdr and ripds
      est_mu_vec <- cbind(ripdr[i], ripds[i])

      # compute the chi-square statistic with the inverse or the generalized inverse of the
      # covariance matrix, and the degrees of freedom
      qf <- quad_form(cov_mat = cov_mat, dev_vec = t(est_mu_vec - mu_vec), df = 2)
      chisq[i] <- qf$stat
      df_chisq[i] <- qf$df
    }
  }

  # find the items whose covariance matrix is singular
  singular_id <- x$id[which(df_chisq < 2 & !(seq_len(nitem) %in% item.skip))]

  # standardize the two statistics of ripdr and ripds
  z_stat_ripdr <- (ripdr - moments_ripdr$mu) / moments_ripdr$sigma
  z_stat_ripds <- (ripds - moments_ripds$mu) / moments_ripds$sigma

  # calculate p-values for all three statistics
  p_ripdr <- 2 * stats::pnorm(q = abs(z_stat_ripdr), mean = 0, sd = 1, lower.tail = FALSE)
  p_ripds <- 2 * stats::pnorm(q = abs(z_stat_ripds), mean = 0, sd = 1, lower.tail = FALSE)
  p_ripdrs <- stats::pchisq(chisq, df = df_chisq, lower.tail = FALSE)

  # compute total sample size
  n_total <- n_foc + n_ref

  # create a data frame to contain the results
  stat_df <-
    data.frame(
      id = x$id,
      ripdr = round(ripdr, 4), z.ripdr = round(z_stat_ripdr, 4),
      ripds = round(ripds, 4), z.ripds = round(z_stat_ripds, 4),
      ripdrs = round(chisq, 4),
      p.ripdr = round(p_ripdr, 4), p.ripds = round(p_ripds, 4), p.ripdrs = round(p_ripdrs, 4),
      n.ref = n_ref, n.foc = n_foc, n.total = n_total, stringsAsFactors = FALSE
    )
  rownames(stat_df) <- NULL

  # when there are items that should be skipped for the IPD analysis
  # insert NAs to the corresponding results of the items
  if (!is.null(item.skip)) {
    stat_df[item.skip, 2:9] <- NA
    p_ripdr[item.skip] <- NA
    p_ripds[item.skip] <- NA
    p_ripdrs[item.skip] <- NA
    moments_ripdr[item.skip, ] <- NA
    moments_ripds[item.skip, ] <- NA
    covar[item.skip] <- NA
  }

  # compute the log p-values to choose the item to be removed in the purification
  log_p <- list(
    ripdr = log(2) + stats::pnorm(q = -abs(z_stat_ripdr), mean = 0, sd = 1, log.p = TRUE),
    ripds = log(2) + stats::pnorm(q = -abs(z_stat_ripds), mean = 0, sd = 1, log.p = TRUE),
    ripdrs = stats::pchisq(chisq, df = df_chisq, lower.tail = FALSE, log.p = TRUE)
  )

  # find the flagged items using the unrounded p-values
  ipd_item_ripdr <- as.numeric(which(p_ripdr <= alpha))
  ipd_item_ripds <- as.numeric(which(p_ripds <= alpha))
  ipd_item_ripdrs <- which(p_ripdrs <= alpha)
  if (length(ipd_item_ripdr) == 0) ipd_item_ripdr <- NULL
  if (length(ipd_item_ripds) == 0) ipd_item_ripds <- NULL
  if (length(ipd_item_ripdrs) == 0) ipd_item_ripdrs <- NULL

  # summarize the results
  rst <- list(
    ipd_stat = stat_df,
    ipd_item = list(ripdr = ipd_item_ripdr, ripds = ipd_item_ripds, ripdrs = ipd_item_ripdrs),
    moments = list(ripdr = moments_ripdr, ripds = moments_ripds), covariance = covar, alpha = alpha,
    log_p = log_p, singular = singular_id
  )

  # return the results
  rst
}
