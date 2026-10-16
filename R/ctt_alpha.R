#' Test-Level Reliability Summary (Cronbach's Alpha)
#'
#' Computes a test-level classical test theory (CTT) reliability summary
#' from scored item response data: Cronbach's alpha (raw and standardized
#' forms), the standard error of measurement (SEM), and the average item
#' difficulty and discrimination.
#'
#' @inheritParams ctt_item
#' @param correct Logical. Accepted for consistency with [ctt_item()] and
#'   [ctt()]. Both the raw and corrected mean item-total correlations are
#'   always reported and `ctt_alpha()` has no flagging step, so `correct`
#'   does not change the result.
#'
#' @details
#' Two forms of Cronbach's alpha are always computed and reported as separate
#' columns:
#'
#' Raw alpha (`alpha`) uses the standard variance-based formula,
#' `alpha = (k / (k - 1)) * (1 - sum(item variances) / total-score
#' variance)`, where k is the number of items. This formula is a general
#' reliability coefficient that applies unchanged to dichotomous and
#' polytomous item scores alike, and it reflects the reliability of the
#' unweighted total score obtained by summing the item scores, which is the
#' score most tests use for reporting and decisions.
#' In everyday terms, raw alpha asks: "if every item kept its own natural
#' scale and spread, how consistently do these items agree with each other?"
#'
#' Standardized alpha (`alpha_std`) is raw alpha computed after standardizing
#' every item with nonzero variance to unit variance. When no item is constant,
#' it equals `alpha_std = (k * r_bar) / (1 + (k - 1) * r_bar)`, where `r_bar`
#' is the average pairwise correlation among the items. A constant item has no
#' defined correlation and cannot be standardized, so it stays in the item
#' count k with zero variance, as in raw alpha, and `r_bar` is the average
#' correlation among the items that vary. When only one item varies,
#' `alpha_std` is 0, as is raw alpha; it is `NA` when no item varies or the
#' standardized total score is constant. In everyday terms, standardized alpha
#' asks: "if every item counted equally regardless of how much it happens to
#' vary in this particular sample, how consistently would these items agree
#' with each other?" Because it removes the influence of any single item's
#' variance, standardized alpha is most useful when items differ substantially
#' in scale or format (e.g., a mix of dichotomous and polytomous items with
#' very different score ranges); when all items share the same scale and format
#' (as with a dichotomous selected-response test scored 0/1), raw and
#' standardized alpha are typically close, and raw alpha remains the more
#' directly interpretable of the two since it matches the reliability of the
#' score actually used in practice. See Cronbach (1951) for the original
#' derivation of coefficient alpha, and Osburn (2000) for a discussion
#' contrasting the raw (covariance-based) and standardized (correlation-based)
#' forms.
#'
#' The standard error of measurement (SEM) is computed as
#' `SEM = SD(total score) * sqrt(1 - alpha)` (using raw alpha), following the
#' standard CTT relationship between test reliability and measurement
#' precision.
#'
#' Average difficulty and average discrimination are the simple means of the
#' per-item values computed by [ctt_item()], ignoring `NA` values (e.g., the
#' discrimination of a constant item). As in [ctt_item()], the raw and
#' corrected item-total correlations are averaged separately and both are
#' reported. `ctt_alpha()` has no flagging step, so `correct` does not
#' change the result.
#'
#' Average difficulty depends on `cats`. When `cats` is not supplied, it is
#' inferred for each item as the observed maximum score plus one, with a
#' minimum of two, as in [ctt_item()]. An item that every examinee scores 0
#' on therefore enters the average with a difficulty of 0. For a polytomous
#' item whose highest score category is not observed in the sample, `cats`
#' is inferred too small and the item's difficulty is biased upward, so
#' supply `cats` explicitly whenever the maximum possible score may not have
#' been observed.
#'
#' @return A one-row data frame with the following columns. Values other than
#'   the counts are rounded to three decimal places.
#' \item{n_examinee}{number of examinees included (after listwise deletion).}
#' \item{n_item}{number of items.}
#' \item{alpha}{Cronbach's alpha, raw (covariance-based) form; `NA` when the
#'   total score has zero variance.}
#' \item{alpha_std}{Cronbach's alpha, standardized (correlation-based) form:
#'   raw alpha computed after standardizing every item with nonzero variance
#'   to unit variance, with a constant item kept in the item count as in
#'   `alpha`; `NA` when it is undefined. It equals
#'   `k * r_bar / (1 + (k - 1) * r_bar)` when no item is constant. See
#'   **Details** for when it differs meaningfully from `alpha`.}
#' \item{sem}{the standard error of measurement (based on raw alpha); `NA`
#'   when `alpha` is `NA`.}
#' \item{mean_difficulty}{average item difficulty.}
#' \item{mean_discrimination_raw}{average raw (uncorrected) item-total
#'   correlation across items with a defined value; `NA` when there is none.}
#' \item{mean_discrimination_corrected}{average corrected (item-excluded)
#'   item-total correlation across items with a defined value; `NA` when
#'   there is none.}
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [ctt_item()], [irtQ::score_resp()]
#'
#' @references
#'   Cronbach, L. J. (1951). Coefficient alpha and the internal structure of
#'   tests. *Psychometrika, 16*(3), 297-334. \doi{10.1007/BF02310555}.
#'
#'   Osburn, H. G. (2000). Coefficient alpha and related internal consistency
#'   reliability coefficients. *Psychological Methods, 5*(3), 343-355.
#'   \doi{10.1037/1082-989X.5.3.343}.
#'
#' @keywords internal
ctt_alpha <- function(data, item.id = NULL, cats = NULL, correct = FALSE,
                       missing = NA) {

  # coerce to a plain data frame so column-wise access (data[[j]]) behaves
  # consistently for matrix/tibble/data.frame input alike
  data <- as.data.frame(data, stringsAsFactors = FALSE)

  # number of items (columns); at least two are required for alpha to be
  # mathematically defined
  n_item <- ncol(data)
  if (n_item < 2L) {
    stop("`data` must contain at least two items to compute alpha.",
         call. = FALSE)
  }

  # validate item.id length up front (also re-validated inside ctt_item(),
  # but checking here gives a clearer error before any other work is done)
  if (!is.null(item.id) && length(item.id) != n_item) {
    stop("length(item.id) must equal ncol(data): one ID per item.",
         call. = FALSE)
  }

  # recode a user-specified missing-value sentinel to NA before analysis,
  # mirroring the `missing` argument convention used elsewhere in irtQ
  if (!is.na(missing)) {
    data[data == missing] <- NA
  }

  # listwise-delete any examinee with a remaining missing response, so alpha
  # and the total-score variance are computed on a common sample
  complete_rows <- stats::complete.cases(data)
  n_dropped <- sum(!complete_rows)
  if (n_dropped > 0L) {
    warning(n_dropped, " examinee(s) with missing item responses were ",
            "excluded listwise from ctt_alpha().", call. = FALSE)
    data <- data[complete_rows, , drop = FALSE]
  }

  # stop when a score or a `cats` value is not valid, and round the scores
  data <- check_ctt_scores(data, cats, item.id)
  n_examinee <- nrow(data)

  # per-item variances and the total-score variance, used in the alpha
  # formula below
  item_var <- vapply(data, stats::var, numeric(1))  # variance of each item
  total <- rowSums(data)                             # total score per examinee
  total_var <- stats::var(total)                     # variance of total score

  # Cronbach's alpha; NA when the total score has zero variance (e.g., every
  # examinee has the same total), since the ratio is then undefined
  alpha <- if (total_var > 0) {
    (n_item / (n_item - 1)) * (1 - sum(item_var) / total_var)
  } else {
    NA_real_
  }

  # standard error of measurement: SD(total score) * sqrt(1 - alpha)
  sem <- if (!is.na(alpha)) stats::sd(total) * sqrt(1 - alpha) else NA_real_

  # standardized alpha: raw alpha computed after standardizing every item
  # with nonzero variance to unit variance; a constant item stays in the item
  # count with zero variance, as in raw alpha above
  # flag the items with nonzero variance, since only these can be standardized
  vary <- item_var > 0
  k_vary <- sum(vary)    # number of items that vary

  # mean off-diagonal correlation among the items that vary; NA when no item
  # varies, and 0 (unused) when one item varies, whose standardized total
  # score has variance 1 and gives a standardized alpha of 0
  r_bar <- NA_real_
  if (k_vary == 1L) {
    r_bar <- 0
  } else if (k_vary >= 2L) {
    item_cor <- stats::cor(data[vary])
    r_bar <- mean(item_cor[upper.tri(item_cor)])
  }

  # denominator of the standardized form; zero when the standardized total
  # score is constant
  denom_std <- 1 + (k_vary - 1) * r_bar

  # standardized alpha, or NA when it is undefined
  alpha_std <- if (!is.na(r_bar) && denom_std > sqrt(.Machine$double.eps)) {
    (n_item / (n_item - 1)) * ((k_vary - 1) * r_bar) / denom_std
  } else {
    NA_real_
  }

  # average item difficulty/discrimination, obtained by calling ctt_item()
  # on the same (already missing-recoded and listwise-deleted) data, so the
  # two functions always agree on the underlying per-item formulas; missing
  # is left at its default (NA) here since recoding already happened above
  item_stats <- ctt_item(data = data, item.id = item.id, cats = cats,
                          correct = correct, missing = NA, flag = FALSE)$item

  # average of the defined values, or NA when none is defined
  mean_or_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  mean_difficulty <- mean_or_na(item_stats$difficulty)
  mean_discrimination_raw <- mean_or_na(item_stats$discrimination_raw)
  mean_discrimination_corrected <- mean_or_na(item_stats$discrimination_corrected)

  # assemble the one-row test-level summary, rounded for readable reporting
  data.frame(
    n_examinee = n_examinee,
    n_item = n_item,
    alpha = round(alpha, 3),
    alpha_std = round(alpha_std, 3),
    sem = round(sem, 3),
    mean_difficulty = round(mean_difficulty, 3),
    mean_discrimination_raw = round(mean_discrimination_raw, 3),
    mean_discrimination_corrected = round(mean_discrimination_corrected, 3)
  )
}
