# This function computes a chi-square fit statistic with its degrees of freedom
# and the critical value
#' @importFrom Rfast rowsums
chisq_stat <- function(exp.freq, obs.freq, count.prm, alpha) {
  # transform the two frequency tables to the matrix forms
  exp.freq2 <- data.matrix(exp.freq)
  obs.freq2 <- data.matrix(obs.freq)

  # replace NA with 0 for the two frequency tables
  exp.freq2[is.na(exp.freq2)] <- 0
  obs.freq2[is.na(obs.freq2)] <- 0

  # create the two proportion tables
  exp.prop <- prop.table(exp.freq2, margin = 1)
  obs.prop <- prop.table(obs.freq2, margin = 1)

  # compute the chi-square statistic
  chisq_fit <- sum(Rfast::rowsums(obs.freq2) * ((obs.prop - exp.prop)^2 / exp.prop), na.rm = TRUE)

  # copy the proportion tables and set the collapsed cells to NA
  exp.prop2 <- exp.prop
  obs.prop2 <- obs.prop
  exp.prop2[is.na(exp.freq)] <- NA
  obs.prop2[is.na(obs.freq)] <- NA

  # compute degrees of freedom
  # the number of collapsed cells
  counted_NA <- sum(is.na(exp.freq))

  # degrees of freedom
  df <- nrow(exp.freq) * (ncol(exp.freq) - 1) - count.prm - counted_NA

  # compute the critical value only when degrees of freedom remain
  crtval_cen <- if (df > 0) stats::qchisq(1 - alpha, df = df, lower.tail = TRUE) else NA_real_

  # return results
  list(
    chisq_fit = chisq_fit, df = df, crtval_cen = crtval_cen,
    exp.prop = exp.prop2, obs.prop = obs.prop2
  )
}
