#' Frequency Distribution Table for Total Scores
#'
#' Computes a frequency distribution table for a vector of total (raw)
#' scores: the frequency, percentage, and cumulative percentage of each score
#' value, as commonly reported with classical test theory (CTT) item
#' analysis results. [ctt()] calls this function to build the total-score
#' frequency distribution in its output, and it can also be used on its own
#' for any vector of integer total scores.
#'
#' @param score A numeric vector of integer total (raw) scores, one value per
#'   examinee. A factor is converted through its labels (for example, the
#'   level `"12"` becomes the score 12), not through its level codes.
#' @param missing A value indicating missing scores in `score`, analogous to
#'   the `missing` argument in [irtQ::est_irt()] and [irtQ::score_resp()]. Any
#'   element equal to `missing` is recoded to `NA` before tabulation. Default
#'   is `NA`. Examinees with a (remaining) missing score are excluded from
#'   the table, with a warning reporting how many were dropped.
#'
#' @details
#' The table includes every integer score value spanning the observed
#' minimum to maximum score (not just the values that were actually
#' observed), so that a score with zero examinees still appears in the table
#' with a frequency of 0, matching how a raw-score frequency table is
#' conventionally reported. Percentages are computed relative to the number
#' of non-missing scores and rounded to two decimal places; cumulative
#' percentages are the running sum of the unrounded percentages, rounded to
#' two decimal places only in the final output, so rounding error does not
#' accumulate and the last cumulative percentage is 100.
#'
#' This function requires integer scores, as is standard for a raw total
#' score. It does not bin or group values, and it stops with an error when a
#' score is not a whole number. A score within 1e-8 of a whole number, such as
#' the result of a floating-point sum, is counted as that whole number.
#'
#' The output table always spans `min(score)` to `max(score)`, so its size
#' scales with the observed score *range*, not the sample size; a single
#' unusually large or small outlier score will produce a correspondingly
#' large table.
#'
#' @return A data frame with one row per score value from `min(score)` to
#'   `max(score)`, containing:
#' \item{score}{the score value.}
#' \item{freq}{the number of examinees with that score.}
#' \item{pct}{the percentage of examinees with that score.}
#' \item{cum_pct}{the cumulative percentage up to and including that score.}
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [ctt()]
#'
#' @examples
#' set.seed(1)
#' score <- rbinom(500, size = 20, prob = 0.6)
#' freq_score(score)
#'
#' @export
freq_score <- function(score, missing = NA) {

  # a factor is converted through its labels, not its level codes
  if (is.factor(score)) score <- as.character(score)

  # coerce to a plain numeric vector; a character score such as "12" becomes
  # a number
  score <- as.numeric(score)

  # recode a user-specified missing-value sentinel to NA before tabulation,
  # mirroring the `missing` argument convention used elsewhere in irtQ
  if (!is.na(missing)) {
    score[score == missing] <- NA
  }

  # drop any missing scores, warning how many were excluded
  n_dropped <- sum(is.na(score))
  if (n_dropped > 0L) {
    warning(n_dropped, " missing score(s) were excluded from freq_score().",
            call. = FALSE)
    score <- score[!is.na(score)]
  }

  # need at least one valid score to build a table
  if (length(score) == 0L) {
    stop("`score` has no non-missing values to tabulate.", call. = FALSE)
  }

  # this function only supports integer-valued scores (see @details); a
  # non-integer score would fall outside the integer bins below, so check
  # every score explicitly here
  if (any(!is.finite(score) | abs(score - round(score)) > 1e-8)) {
    stop("`score` must contain only integer-valued scores; freq_score() ",
         "does not bin or group non-integer values.", call. = FALSE)
  }

  # round each score so that it matches its integer bin exactly
  score <- round(score)

  # full range of integer score values (including any with zero observed
  # frequency), from the observed minimum to the observed maximum
  score_range <- seq.int(min(score), max(score))

  # count of examinees at each value in the full range; matching against
  # score_range (rather than using table() alone) ensures unobserved values
  # in the middle of the range still appear with frequency 0
  freq <- vapply(score_range, function(s) sum(score == s), integer(1))

  # percentage of (non-missing) examinees at each score value
  n_total <- length(score)
  pct <- 100 * freq / n_total

  # cumulative percentage: running sum of the unrounded percentages, with
  # rounding applied only once at the very end (see @details)
  cum_pct <- cumsum(pct)

  # assemble the frequency table, rounding percentages for readable reporting
  data.frame(
    score = score_range,
    freq = freq,
    pct = round(pct, 2),
    cum_pct = round(cum_pct, 2)
  )
}
