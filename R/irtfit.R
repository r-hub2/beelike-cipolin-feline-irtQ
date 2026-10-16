#' Traditional IRT Item Fit Statistics
#'
#' This function computes traditional IRT item fit statistics, including the
#' \eqn{\chi^{2}} fit statistic (e.g., Bock, 1960; Yen, 1981),
#' the log-likelihood ratio \eqn{\chi^{2}} fit statistic (\eqn{G^{2}}; McKinley
#' & Mills, 1985), and the infit and outfit statistics (Ames & Penfield, 2015). It
#' also returns contingency tables used to compute the \eqn{\chi^{2}} and
#' \eqn{G^{2}} statistics.
#'
#' @inheritParams est_score
#' @param x A data frame containing item metadata (e.g., item parameters,
#'   number of categories, IRT model types, etc.), or an object of class
#'   `est_irt` or `est_item` obtained from [irtQ::est_irt()] or
#'   [irtQ::est_item()].
#'
#'   See [irtQ::est_irt()] or [irtQ::simdat()] for more details about the item
#'   metadata. This data frame can be easily created using the
#'   [irtQ::shape_df()] function.
#' @param score A numeric vector of examinees' ability estimates (theta
#'   values). Examinees with a missing ability estimate are excluded. Not used
#'   by the `est_item` method, which takes the ability estimates stored in `x`.
#' @param group.method A character string specifying the method used to group
#'   examinees along the ability scale when computing the \eqn{\chi^{2}} and
#'   \eqn{G^{2}} fit statistics. Available options are:
#'    - `"equal.width"`: Divides the ability range into intervals of equal width.
#'    - `"equal.freq"`: Divides the examinees into groups of approximately
#'    equal size.
#'
#'   Default is `"equal.width"`. For each item, the intervals span the minimum
#'   and maximum ability estimates (after truncation by `range.score`) of the
#'   examinees who responded to the item. The number of intervals is set by
#'   `n.width`.
#' @param n.width An integer specifying the number of intervals (groups) into
#'   which the ability scale is divided. Duplicate cut points and empty
#'   intervals are dropped, so an item can have fewer groups. Default is 10.
#' @param loc.theta A character string indicating the point on the ability scale
#'  at which the expected category probabilities are calculated for each group.
#'  Available options are:
#'  - `"average"`: Uses the average ability estimate of examinees within each group.
#'  - `"middle"`: Uses the midpoint of each group's ability interval.
#'
#'  Default is `"average"`.
#' @param range.score A numeric vector of length two giving the lower and upper
#'   bounds of the ability scale. Ability estimates outside these bounds are set
#'   to the nearest bound before grouping and before computing infit and outfit.
#'   The bounds also set the theta range of the plots drawn by
#'   [irtQ::plot.irtfit()]. If `NULL`, no truncation is applied, and the plot
#'   range is \eqn{\pm} the largest absolute ability estimate rounded up to an
#'   integer. Default is `NULL`.
#' @param alpha A numeric value specifying the significance level (\eqn{\alpha})
#'   for the tests of the \eqn{\chi^{2}} and \eqn{G^{2}} item fit statistics. It
#'   also sets the confidence level, \eqn{1 - \alpha}, of the intervals drawn by
#'   [irtQ::plot.irtfit()]. Default is `0.05`.
#' @param missing A value indicating missing responses in `data`. Default is
#'   `NA`. Missing responses are excluded item by item.
#' @param overSR A positive numeric threshold for the absolute standardized
#'   residuals. The proportion of cells (ability groups by score categories),
#'   before collapsing, whose absolute standardized residuals exceed this value
#'   is reported in `overSR.prop`. Default is 2.
#' @param min.collapse A numeric value giving the minimum expected frequency per
#'   cell. When computing \eqn{\chi^{2}} and \eqn{G^{2}}, adjacent ability groups
#'   are merged until every expected cell frequency is at least this value.
#'   Default is 1.
#' @param pcm.loc An integer vector giving the positions (rows of `x`) of
#'   partial credit model (PCM) items whose slope parameters are fixed. It is
#'   used only to count the item parameters for the degrees of freedom of
#'   \eqn{\chi^{2}}. Default is `NULL`.
#' @param ... Further arguments passed to or from other methods.
#'
#' @details
#' To compute the \eqn{\chi^2} and \eqn{G^2} item fit statistics, the `group.method`
#' argument determines how the ability scale is divided into groups:
#' - `"equal.width"`: Examinees are grouped by intervals of equal width along
#'  the ability scale.
#' - `"equal.freq"`: Examinees are grouped such that each group contains
#'  (approximately) the same number of individuals.
#'
#' Note that `"equal.freq"` does not guarantee *exactly* equal frequencies across
#' all groups, since grouping is based on quantiles.
#'
#' When dividing the ability scale into intervals to compute the \eqn{\chi^2}
#' and \eqn{G^2} fit statistics, the intervals should be:
#' - **Wide enough** to ensure that each group contains a sufficient number of
#'  examinees (to avoid unstable estimates),
#' - **Narrow enough** to ensure that examinees within each group are relatively
#'  homogeneous in ability (Hambleton et al., 1991).
#'
#' Use `n.width` to change the number of groups (default 10). For reference:
#' - Yen (1981) used 10 groups of approximately equal size,
#' - Bock (1960) allowed for flexibility in the number of groups.
#'
#' With `loc.theta = "average"`, the expected probabilities of each group are
#' evaluated at the average ability estimate of the group. The probability at
#' this point approximates the average of the model probabilities of the
#' examinees in the group (e.g., Yen, 1981).
#'
#' Regarding degrees of freedom (*df*), let \eqn{G} be the number of ability
#' groups after collapsing, \eqn{K} the number of score categories, and
#' \eqn{m} the number of item parameters:
#' - The \eqn{\chi^2} statistic is approximately chi-square distributed with
#'   \eqn{G(K - 1) - m} degrees of freedom, which is \eqn{G - m} for a
#'   dichotomous item (Ames & Penfield, 2015).
#' - The \eqn{G^2} statistic is approximately chi-square distributed with
#'   \eqn{G(K - 1)} degrees of freedom, which is \eqn{G} for a dichotomous item
#'   (Ames & Penfield, 2015; Muraki & Bock, 2003).
#'
#' The number of item parameters \eqn{m} is determined by the IRT model of each
#' item: 1 for the 1PLM, 2 for the 2PLM, 3 for the 3PLM, and the number of score
#' categories \eqn{K} for the GRM and GPCM; an item listed in `pcm.loc` has
#' \eqn{m = K - 1}, and an item labeled `"DRM"` is counted as a 3PLM item. This
#' also applies to `est_irt` and `est_item` objects: parameters fixed during
#' estimation (for example, with `fix.g = TRUE`, or the items fixed in FIPC) are
#' still counted, and GPCM items estimated with `fix.a.gpcm = TRUE` are counted
#' as PCM items only when they are given in `pcm.loc`. When no degrees of
#' freedom remain for \eqn{\chi^2} after collapsing (`df.X2` of 0 or less),
#' `crit.val.X2` and `p.X2` are `NA`, and a warning names the item.
#'
#' Responses in `data` must be integer scores from 0 to the number of score
#' categories minus 1, with one column per item. Character and factor responses
#' are read as numbers, logical responses as 0 and 1, and other responses stop
#' the function with an error.
#'
#' For ability group \eqn{j} with \eqn{N_j} examinees, let \eqn{O_{jk}} and
#' \eqn{E_{jk}} be the observed proportion and the model-expected probability
#' of score category \eqn{k} at the theta point of the group (see
#' `loc.theta`). The statistics are
#' \eqn{\chi^2 = \sum_j \sum_k N_j (O_{jk} - E_{jk})^2 / E_{jk}} and
#' \eqn{G^2 = 2 \sum_j \sum_k N_j O_{jk} \log(O_{jk} / E_{jk})}, where cells
#' with \eqn{O_{jk} = 0} add nothing to \eqn{G^2}. For a dichotomous item,
#' \eqn{\chi^2} has the form of Yen's (1981) \eqn{Q_1} statistic
#' \eqn{\sum_j N_j (O_j - E_j)^2 / [E_j (1 - E_j)]}. The standardized
#' residual of a cell before collapsing is
#' \eqn{(O_{jk} - E_{jk}) / \sqrt{E_{jk} (1 - E_{jk}) / N_j}}.
#'
#' Infit and outfit use the category scores \eqn{0, \ldots, K - 1}. For
#' examinee \eqn{n}, the residual is \eqn{r_n = x_n - E(X \mid \theta_n)} and
#' the variance is \eqn{W_n = \sum_k (k - E(X \mid \theta_n))^2 P_k(\theta_n)},
#' where \eqn{\theta_n} is the ability estimate in `score`. Outfit is the
#' unweighted mean \eqn{\sum_n (r_n^2 / W_n) / N}, and infit is the
#' variance-weighted mean \eqn{\sum_n r_n^2 / \sum_n W_n}.
#'
#' Note that infit and outfit statistics should be interpreted with caution when
#' applied to non-Rasch models. The returned object, in particular its
#' contingency tables, can be passed to [irtQ::plot.irtfit()] to draw raw and
#' standardized residual plots (Hambleton et al., 1991).
#'
#' @return This function returns an object of class `irtfit`, which includes
#' the following components:
#'
#'   \item{fit_stat}{A data frame with one row per item and the columns `id`,
#'   `X2` (\eqn{\chi^{2}}), `G2` (\eqn{G^{2}}), `df.X2`, `df.G2`, `crit.val.X2`
#'   and `crit.val.G2` (critical values at `alpha`), `p.X2`, `p.G2`, `outfit`,
#'   `infit`, `N` (the number of examinees who responded to the item and have
#'   an ability estimate), and
#'   `overSR.prop` (the proportion of cells, before collapsing, whose absolute
#'   standardized residuals exceed `overSR`).}
#'
#'   \item{contingency.fitstat}{A list of contingency tables, one per item, used
#'   to compute \eqn{\chi^{2}} and \eqn{G^{2}}. Adjacent ability groups are
#'   merged until every expected frequency is at least `min.collapse`, the rows
#'   are numbered, and the interval labels are dropped. The columns are
#'   `total` and, for each score category, `obs.freq.*`, `exp.freq.*`,
#'   `obs.prop.*`, `exp.prob.*`, and `raw.rsd.*`.}
#'
#'   \item{contingency.plot}{A list of contingency tables, one per item, used by
#'   [irtQ::plot.irtfit()] to draw raw and standardized residual plots
#'   (Hambleton et al., 1991). The tables use the uncollapsed groups and contain
#'   `interval`, `point` (the theta point of the group), `total`, and, for each
#'   score category, `obs.freq.*`, `obs.prop.*`, `exp.prob.*`, `raw.rsd.*`,
#'   `se.*`, and `std.rsd.*`.}
#'
#'   \item{item_df}{The item metadata used in the analysis. Items with fewer
#'   than two responses are removed.}
#'
#'   \item{individual.info}{A list of data frames, one per item, with the
#'   residual (`resid`) and variance (`Var`) for each examinee who responded to
#'   the item. These values are used to compute infit and outfit.}
#'
#'   \item{ancillary}{A list with `range.score`, `alpha`, `overSR`, and
#'   `scale.D`, used by [irtQ::plot.irtfit()].}
#'
#'   \item{call}{The matched call.}
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [irtQ::plot.irtfit()], [irtQ::shape_df()], [irtQ::est_irt()],
#' [irtQ::est_item()]
#'
#' @references Ames, A. J., & Penfield, R. D. (2015). An NCME instructional
#'   module on item-fit statistics for item response theory models. *Educational
#'   Measurement: Issues and Practice, 34*(3), 39-48.
#'
#'   Bock, R. D. (1960). *Methods and applications of optimal scaling*. Chapel
#'   Hill, NC: L. L. Thurstone Psychometric Laboratory.
#'
#'   Hambleton, R. K., Swaminathan, H., & Rogers, H. J. (1991). *Fundamentals of
#'   item response theory*. Newbury Park, CA: Sage.
#'
#'   McKinley, R., & Mills, C. (1985). A comparison of several goodness-of-fit
#'   statistics. *Applied Psychological Measurement, 9*, 49-57.
#'
#'   Muraki, E., & Bock, R. D. (2003). PARSCALE 4: IRT item analysis and test
#'   scoring for rating scale data (Computer software). Chicago, IL: Scientific
#'   Software International. URL http://www.ssicentral.com
#'
#'   Yen, W. M. (1981). Using simulation results to choose a latent trait model.
#'   *Applied Psychological Measurement, 5*, 245-262.
#'
#' @examples
#' \donttest{
#' ## Example 1
#' ## Use the simulated CAT data
#' # Identify items whose summed item scores exceed 10,000
#' over10000 <- which(colSums(simCAT_MX$res.dat, na.rm = TRUE) > 10000)
#'
#' # Select these items
#' x <- simCAT_MX$item.prm[over10000, ]
#'
#' # Extract response data for the selected items
#' data <- simCAT_MX$res.dat[, over10000]
#'
#' # Extract examinees' ability estimates
#' score <- simCAT_MX$score
#'
#' # Compute item fit statistics
#' fit1 <- irtfit(
#'   x = x, score = score, data = data, group.method = "equal.width",
#'   n.width = 10, loc.theta = "average", range.score = NULL, D = 1, alpha = 0.05,
#'   missing = NA, overSR = 2
#' )
#'
#' # View the fit statistics
#' fit1$fit_stat
#'
#' # View the contingency tables used to compute fit statistics
#' fit1$contingency.fitstat
#'
#'
#' ## Example 2
#' ## Import the "-prm.txt" output file from flexMIRT
#' flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
#'
#' # Select the first two dichotomous items and the last polytomous item
#' x <- bring.flexmirt(file = flex_sam, "par")$Group1$full_df[c(1:2, 55), ]
#'
#' # Generate ability values from a standard normal distribution
#' set.seed(10)
#' score <- rnorm(1000, mean = 0, sd = 1)
#'
#' # Simulate response data
#' data <- simdat(x = x, theta = score, D = 1)
#'
#' # Compute item fit statistics
#' fit2 <- irtfit(
#'   x = x, score = score, data = data, group.method = "equal.freq",
#'   n.width = 11, loc.theta = "average", range.score = c(-4, 4), D = 1, alpha = 0.05
#' )
#'
#' # View the fit statistics
#' fit2$fit_stat
#'
#' # View the contingency tables used to compute fit statistics
#' fit2$contingency.fitstat
#'
#' # Plot raw and standardized residuals for the first item (dichotomous)
#' plot(x = fit2, item.loc = 1, type = "both", ci.method = "wald",
#'      show.table = TRUE, ylim.sr.adjust = TRUE)
#'
#' # Plot raw and standardized residuals for the third item (polytomous)
#' plot(x = fit2, item.loc = 3, type = "both", ci.method = "wald",
#'      show.table = FALSE, ylim.sr.adjust = TRUE)
#' }
#'
#' @export
irtfit <- function(x, ...) UseMethod("irtfit")

#' @describeIn irtfit Default method for computing traditional IRT item fit
#' statistics using a data frame `x` that contains item metadata.
#' @import dplyr
#' @export
#'
irtfit.default <- function(x,
                           score,
                           data,
                           group.method = c("equal.width", "equal.freq"),
                           n.width = 10,
                           loc.theta = "average",
                           range.score = NULL,
                           D = 1,
                           alpha = 0.05,
                           missing = NA,
                           overSR = 2,
                           min.collapse = 1,
                           pcm.loc = NULL,
                           ...) {

  # match.call
  cl <- match.call()

  # confirm and correct all item metadata information
  x <- confirm_df(x)

  # check the location of the theta point for each group
  loc.theta <- tolower(loc.theta)
  if (!(length(loc.theta) == 1L && loc.theta %in% c("average", "middle"))) {
    stop("'loc.theta' must be either \"average\" or \"middle\".", call. = FALSE)
  }

  # create a vector of PCM item indicators
  pcm.lg <- logical(nrow(x))
  pcm.lg[pcm.loc] <- TRUE

  # transform scores to a vector form
  if (is.matrix(score) | is.data.frame(score)) {
    score <- as.numeric(data.matrix(score))
  }

  # recode missing values
  if (!is.na(missing)) {
    data[data == missing] <- NA
  }

  # transform the response data to a numeric matrix and check the responses
  data <- resp_to_matrix(data, x$cats, x$id)

  # stop when the scores do not match the rows of the response data
  if (length(score) != nrow(data)) {
    stop("The length of 'score' must equal the number of rows in 'data'.", call. = FALSE)
  }

  # check if there are items which have zero or one response frequency
  n.score <- Rfast::colsums(!is.na(data))
  if (all(n.score %in% c(0L, 1L))) {
    stop("Every item has fewer than two responses. Each item must have at least two responses.", call. = FALSE)
  }

  if (any(n.score %in% c(0L, 1L))) {
    del_item <- which(n.score %in% c(0L, 1L))

    # remove items with fewer than two responses
    x_del <- x[del_item, ]
    x <- x[-del_item, ]
    data <- data[, -del_item, drop = FALSE]
    pcm.lg <- pcm.lg[-del_item]

    # warning message with the column numbers and ids of the removed items
    memo <- paste0(
      "The following items have fewer than two responses and are excluded from the analysis: ",
      paste0(x_del$id, " (column ", del_item, ")", collapse = ", "), "."
    )
    warning(memo, call. = FALSE)
  }

  # stop when the bounds of the score range are not valid
  if (!is.null(range.score) &&
    !isTRUE(is.numeric(range.score) && length(range.score) == 2L && range.score[1] < range.score[2])) {
    stop(
      "'range.score' must be a numeric vector of length two ",
      "with the first element smaller than the second element.",
      call. = FALSE
    )
  }

  # restrict the range of scores if required
  if (!is.null(range.score)) {
    score <- ifelse(score < range.score[1], range.score[1], score)
    score <- ifelse(score > range.score[2], range.score[2], score)
  } else {
    tmp.val <- max(ceiling(abs(range(score, na.rm = TRUE))))
    range.score <- c(-tmp.val, tmp.val)
  }

  # compute item fit statistics and obtain contingency tables across all items
  fits <-
    lapply(seq_len(nrow(x)), FUN = function(i) {
      itemfit(
        x_item = x[i, ], score = score, resp = data[, i], group.method = group.method,
        n.width = n.width, loc.theta = loc.theta, D = D, alpha = alpha, overSR = overSR,
        min.collapse = min.collapse, is.pcm = pcm.lg[i]
      )
    })

  # extract fit statistics
  fit_stat <-
    purrr::map(fits, .f = function(i) i$fit.stats) %>%
    do.call(what = "rbind")
  fit_stat <- data.frame(id = x$id, fit_stat)

  # extract the contingency tables used to compute the fit statistics
  contingency.fitstat <-
    purrr::map(fits, .f = function(i) i$contingency.fitstat)
  names(contingency.fitstat) <- x$id

  # extract the contingency tables to be used to draw residual plots
  contingency.plot <-
    purrr::map(fits, .f = function(i) i$contingency.plot)
  names(contingency.plot) <- x$id

  # extract the individual residuals and variances
  individual.info <-
    purrr::map(fits, .f = function(i) i$individual.info)
  names(individual.info) <- x$id

  # warn when an item has no degrees of freedom left for the chi-square test
  if (any(fit_stat$df.X2 <= 0)) {
    warning(
      "No degrees of freedom remain for the X2 statistic of item(s) ",
      paste(fit_stat$id[fit_stat$df.X2 <= 0], collapse = ", "),
      ". Their critical values and p-values are set to NA.",
      call. = FALSE
    )
  }

  # return results
  rst <- list(
    fit_stat = fit_stat, contingency.fitstat = contingency.fitstat,
    contingency.plot = contingency.plot,
    item_df = x, individual.info = individual.info,
    ancillary = list(range.score = range.score, alpha = alpha, overSR = overSR, scale.D = D)
  )
  class(rst) <- "irtfit"
  rst$call <- cl

  rst
}


#' @describeIn irtfit Method for an object of class `est_item` created by
#' [irtQ::est_item()]. The item parameter estimates, response data, ability
#' estimates, and scaling constant `D` are taken from the object.
#' @import dplyr
#' @export
#'
irtfit.est_item <- function(x,
                            group.method = c("equal.width", "equal.freq"),
                            n.width = 10,
                            loc.theta = "average",
                            range.score = NULL,
                            alpha = 0.05,
                            missing = NA,
                            overSR = 2,
                            min.collapse = 1,
                            pcm.loc = NULL,
                            ...) {
  # match.call
  cl <- match.call()

  # compute the fit statistics with the scores, data, and scaling constant stored in the object
  rst <- irtfit.default(
    x = x$par.est, score = x$score, data = x$data, group.method = group.method,
    n.width = n.width, loc.theta = loc.theta, range.score = range.score, D = x$scale.D,
    alpha = alpha, missing = missing, overSR = overSR, min.collapse = min.collapse,
    pcm.loc = pcm.loc
  )

  # keep the call of this method
  rst$call <- cl

  rst
}


#' @describeIn irtfit Method for an object of class `est_irt` created by
#' [irtQ::est_irt()]. The item parameter estimates, response data, and scaling
#' constant `D` are taken from the object, and `score` must be supplied.
#' @import dplyr
#' @export
#'
irtfit.est_irt <- function(x,
                           score,
                           group.method = c("equal.width", "equal.freq"),
                           n.width = 10,
                           loc.theta = "average",
                           range.score = NULL,
                           alpha = 0.05,
                           missing = NA,
                           overSR = 2,
                           min.collapse = 1,
                           pcm.loc = NULL,
                           ...) {
  # match.call
  cl <- match.call()

  # compute the fit statistics with the data and scaling constant stored in the object
  rst <- irtfit.default(
    x = x$par.est, score = score, data = x$data, group.method = group.method,
    n.width = n.width, loc.theta = loc.theta, range.score = range.score, D = x$scale.D,
    alpha = alpha, missing = missing, overSR = overSR, min.collapse = min.collapse,
    pcm.loc = pcm.loc
  )

  # keep the call of this method
  rst$call <- cl

  rst
}



#' @importFrom janitor adorn_totals
#' @importFrom tibble rownames_to_column remove_rownames column_to_rownames
#' @import dplyr
itemfit <- function(x_item, score, resp, group.method = c("equal.width", "equal.freq"),
                    n.width = 10, loc.theta = "average", D = 1, alpha = 0.05, overSR = 2, min.collapse = 1,
                    is.pcm = FALSE) {
  # break down the item metadata into several elements
  elm_item <- breakdown(x_item)

  # number of score categories
  cats <- elm_item$cats

  # delete missing responses and missing thetas
  na_lg <- is.na(resp) | is.na(score)
  resp <- resp[!na_lg]
  score <- score[!na_lg]

  # assign factor levels to the item responses
  resp <- factor(resp, levels = (seq_len(cats) - 1))

  # compute cut scores to divide score groups
  group.method <- match.arg(group.method)
  cutscore <- switch(group.method,
    equal.width = seq(from = min(score), to = max(score), length.out = n.width + 1),
    equal.freq = stats::quantile(score, probs = seq(0, 1, length.out = n.width + 1), type = 9)
  )

  # drop duplicate cut scores, which can occur with equal.freq
  cutscore <- unique(cutscore)

  # stop when the ability estimates do not span an interval
  if (length(cutscore) < 2L) {
    stop(
      "The ability estimates of the examinees who responded to item ", x_item$id,
      " have no variation.",
      call. = FALSE
    )
  }

  # assign score group variable to each score
  intv <- cut(score, breaks = cutscore, right = FALSE, include.lowest = TRUE, dig.lab = 7)

  # create a contingency table for the frequencies of score points
  obs.freq <-
    as.data.frame.matrix(table(intv, resp)) %>%
    tibble::rownames_to_column(var = "interval") %>%
    janitor::adorn_totals(where = "col", name = "total")
  delrow.lg <- obs.freq$total == 0
  obs.freq <- obs.freq[!delrow.lg, ]

  # create a contingency table for the category proportions
  obs.prop <- obs.freq

  # divide the category frequencies by the row totals
  obs.prop[, 2:(cats + 1)] <- obs.freq[, 2:(cats + 1)] / obs.freq$total

  # find a theta point for each score group
  loc.theta <- tolower(loc.theta)
  if (loc.theta == "middle") {
    theta <- purrr::map_dbl(.x = 2:length(cutscore), .f = function(i) mean(c(cutscore[i - 1], cutscore[i])))
    theta <- theta[!delrow.lg]
  } else {
    theta <-
      data.frame(intv, score) %>%
      dplyr::group_by(.data$intv) %>%
      dplyr::summarize(ave = mean(.data$score), .groups = "drop") %>%
      dplyr::pull(2)
  }

  # compute expected probabilities of endorsing an answer to each score category
  exp.prob <-
    trace(elm_item = elm_item, theta = theta, D = D, tcc = FALSE)$prob.cats[[1]] %>%
    data.frame()
  colnames(exp.prob) <- 0:(cats - 1)

  ## -------------------------------------------------------------------------
  # compute raw residuals (rr)
  rr <- obs.prop[2:(cats + 1)] - exp.prob

  # compute the standard errors (se)
  se <- sqrt((exp.prob * (1 - exp.prob)) / obs.freq$total)

  # standardize the residuals (sr)
  sr <- rr / se

  # compute the proportion of cells whose absolute standardized residuals
  # exceed overSR
  over_sr <- sum(abs(sr) > overSR)
  over_sr_prop <- round(over_sr / (cats * nrow(sr)), 3)

  # create a full contingency table to draw IRT residual plots
  ctg_tb <-
    data.frame(
      point = theta, obs.freq = obs.freq, obs.prop = obs.prop[, 2:(cats + 1)],
      exp.prob = exp.prob, raw.rsd = rr, se = se, std.rsd = sr
    ) %>%
    dplyr::relocate("obs.freq.interval", .before = "point") %>%
    dplyr::relocate("obs.freq.total", .after = "point") %>%
    dplyr::rename("interval" = "obs.freq.interval", "total" = "obs.freq.total")

  ## ------------------------------------------------------------------------------
  # collapsing the contingency tables to compute the chi-square fit statistics
  # check the number of expected frequency for all cells
  exp.freq <- exp.prob * obs.freq$total

  # collapse the expected and observed frequency tables
  ftable_info <-
    data.frame(exp.freq = exp.freq, obs.freq = obs.freq) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(var = "obs.freq.interval")
  for (i in 1:cats) {
    ftable_info <- collapse_ftable(x = ftable_info, col = i, min.collapse = min.collapse)
  }

  # new contingency tables after collapsing
  exp.freq.cp <- dplyr::select(ftable_info, dplyr::contains("exp.freq"))
  obs.freq.cp <- dplyr::select(ftable_info, dplyr::contains("obs.freq"))
  freq.tot.cp <- obs.freq.cp$obs.freq.total
  exp.prob.cp <- exp.freq.cp / freq.tot.cp
  obs.prop.cp <- (obs.freq.cp / freq.tot.cp)[, 1:cats]
  colnames(exp.prob.cp) <- paste0("exp.prob.", 0:(cats - 1))
  colnames(obs.prop.cp) <- paste0("obs.prop.", 0:(cats - 1))

  ## -------------------------------------------------------------------------
  # create a contingency table to compute the item fit statistics
  # first, compute raw residuals
  rr.cp <- obs.prop.cp - exp.prob.cp
  colnames(rr.cp) <- paste0("raw.rsd.", 0:(cats - 1))

  # create a full contingency table for chi-square fit statistic
  ctg_tb.cp <-
    data.frame(obs.freq.cp, exp.freq.cp, obs.prop.cp, exp.prob.cp, rr.cp) %>%
    dplyr::relocate("obs.freq.total", .before = "obs.freq.0") %>%
    dplyr::rename("total" = "obs.freq.total")
  rownames(ctg_tb.cp) <- 1:nrow(ctg_tb.cp)

  # compute the chi-square statistic (X2)
  x2 <- sum(freq.tot.cp * (rr.cp^2 / exp.prob.cp), na.rm = TRUE)

  # compute the likelihood ratio chi-square fit statistic (G2)
  g2 <- 2 * sum(obs.freq.cp[, 1:cats] * log(obs.prop.cp / exp.prob.cp), na.rm = TRUE)

  # find the number of parameters for each item
  model <- x_item$model
  if (is.pcm) model <- "PCM"
  count_prm <- NA
  count_prm[model %in% "1PLM"] <- 1
  count_prm[model %in% "2PLM"] <- 2
  count_prm[model %in% c("3PLM", "DRM")] <- 3
  count_prm[model %in% "PCM"] <- x_item[model %in% "PCM", 2] - 1
  count_prm[model %in% "GPCM"] <- x_item[model %in% "GPCM", 2]
  count_prm[model %in% "GRM"] <- x_item[model %in% "GRM", 2]

  # find a critical value and compute the p values
  df.x2 <- nrow(exp.freq.cp) * (ncol(exp.freq.cp) - 1) - count_prm
  df.g2 <- nrow(exp.freq.cp) * (ncol(exp.freq.cp) - 1)
  crtval.x2 <- if (df.x2 > 0) stats::qchisq(1 - alpha, df = df.x2, lower.tail = TRUE) else NA_real_
  crtval.g2 <- if (df.g2 > 0) stats::qchisq(1 - alpha, df = df.g2, lower.tail = TRUE) else NA_real_
  pval.x2 <- if (df.x2 > 0) 1 - stats::pchisq(x2, df = df.x2, lower.tail = TRUE) else NA_real_
  pval.g2 <- if (df.g2 > 0) 1 - stats::pchisq(g2, df = df.g2, lower.tail = TRUE) else NA_real_

  ## ------------------------------------------------------------------------------
  # infit & outfit
  # individual expected probabilities for each score category
  indiv_exp.prob <- trace(
    elm_item = elm_item, theta = score, D = D,
    tcc = FALSE
  )$prob.cats[[1]]

  # indicator matrix of the observed score category of each examinee
  n.resp <- length(resp)
  indiv_obs.prob <-
    table(1:n.resp, resp) %>%
    as.data.frame.matrix()

  # matrix of the category scores (0 to cats - 1)
  Emat <- matrix(0:(cats - 1), nrow(indiv_exp.prob), ncol(indiv_exp.prob), byrow = TRUE)

  # score residuals (observed minus expected item score)
  resid <- rowSums(indiv_obs.prob * Emat) - rowSums(Emat * indiv_exp.prob)

  # conditional variance of the item score
  Var <- rowSums((Emat - rowSums(Emat * indiv_exp.prob))^2 * indiv_exp.prob)

  # compute outfit & infit
  outfit <- sum(resid^2 / Var) / n.resp
  infit <- sum(resid^2) / sum(Var)

  ## ------------------------------------------------------------------------------
  # summary of fit statistics
  fitstats <- data.frame(
    X2 = round(x2, 3), G2 = round(g2, 3), df.X2 = df.x2, df.G2 = df.g2,
    crit.val.X2 = round(crtval.x2, 2), crit.val.G2 = round(crtval.g2, 2),
    p.X2 = round(pval.x2, 3), p.G2 = round(pval.g2, 3),
    outfit = round(outfit, 3), infit = round(infit, 3),
    N = n.resp, overSR.prop = over_sr_prop
  )

  ## ------------------------------------------------------------------------------
  # return results
  list(
    fit.stats = fitstats, contingency.fitstat = ctg_tb.cp, contingency.plot = ctg_tb,
    individual.info = data.frame(resid = resid, Var = Var)
  )
}
