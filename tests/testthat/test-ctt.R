# ctt_item() and ctt_alpha() are internal (not exported); access via :::,
# following the same convention as test-confirm_df.R. ctt(), print.ctt(), and
# summary.ctt() are exported and tested directly.
ctt_item <- irtQ:::ctt_item
ctt_alpha <- irtQ:::ctt_alpha

# ---- ctt_item(): dichotomous -------------------------------------------------

test_that("ctt_item() computes difficulty and item-total correlations matching a from-scratch calculation", {
  set.seed(1)
  dat <- data.frame(matrix(rbinom(50 * 6, 1, 0.6), nrow = 50))
  out <- ctt_item(data = dat, flag = FALSE)$item

  total <- rowSums(dat)
  for (j in seq_len(ncol(dat))) {
    expect_equal(out$difficulty[j], round(mean(dat[[j]]), 3))
    expect_equal(out$discrimination_raw[j], round(stats::cor(dat[[j]], total), 3))
    expect_equal(out$discrimination_corrected[j],
                 round(stats::cor(dat[[j]], total - dat[[j]]), 3))
  }
})

test_that("ctt_item() generalizes correctly to polytomous items", {
  set.seed(2)
  dat <- data.frame(matrix(sample(0:3, 40 * 4, replace = TRUE), nrow = 40))
  out <- ctt_item(data = dat, cats = rep(4, 4), flag = FALSE)$item

  expect_equal(out$difficulty, unname(round(vapply(dat, mean, numeric(1)) / 3, 3)))
  expect_true(all(out$cats == 4))
})

test_that("ctt_item() flags extreme difficulty as documented, at the correct boundary", {
  dat <- data.frame(
    easy = c(rep(1, 19), 0),   # difficulty = 19/20 = 0.95 -> boundary, NOT > 0.95
    hard = c(rep(0, 19), 1),   # difficulty = 1/20  = 0.05 -> < 0.10, should flag
    ok   = rep(c(1, 0), 10)    # difficulty = 0.50 -> within bounds
  )
  out <- ctt_item(data = dat, item.id = c("easy", "hard", "ok"),
                   crit.p = c(0.10, 0.95), crit.dis = 0.20)$item
  expect_match(out$flag[out$item == "hard"], "difficulty too low")
  expect_false(grepl("difficulty", out$flag[out$item == "easy"]))
  expect_false(grepl("difficulty", out$flag[out$item == "ok"]))
})

test_that("ctt_item() uses custom item.id and does not fall back to colnames", {
  dat <- data.frame(foo = c(1, 0, 1, 0), bar = c(1, 1, 0, 0))
  out <- ctt_item(data = dat, item.id = c("Q1", "Q2"), flag = FALSE)$item
  expect_equal(out$item, c("Q1", "Q2"))

  out_default <- ctt_item(data = dat, flag = FALSE)$item
  expect_equal(out_default$item, c("V1", "V2"))  # not "foo"/"bar"
})

test_that("ctt_item() errors when item.id length does not match ncol(data)", {
  dat <- data.frame(a = c(1, 0), b = c(0, 1))
  expect_error(ctt_item(data = dat, item.id = c("Q1", "Q2", "Q3")),
               "length\\(item.id\\) must equal ncol\\(data\\)")
})

test_that("ctt_item() excludes examinees with missing responses listwise, with a warning", {
  dat <- data.frame(I1 = c(1, 0, NA, 1), I2 = c(1, 1, 0, 0), I3 = c(0, 1, 1, 1))
  expect_warning(out <- ctt_item(data = dat, flag = FALSE), "excluded listwise")
  expect_equal(nrow(out$item), 3)  # still 3 items
})

test_that("ctt_item() recodes a custom missing sentinel before listwise deletion", {
  dat <- data.frame(I1 = c(1, 0, -9, 1), I2 = c(1, 1, 0, 0), I3 = c(0, 1, 1, 1))
  expect_warning(out <- ctt_item(data = dat, missing = -9, flag = FALSE),
                 "excluded listwise")
  expect_equal(out$item$difficulty, ctt_item(
    data = data.frame(I1 = c(1, 0, 1), I2 = c(1, 1, 0), I3 = c(0, 1, 1)),
    flag = FALSE
  )$item$difficulty)
})

test_that("ctt_item() errors with fewer than two items", {
  dat <- data.frame(I1 = c(1, 0, 1))
  expect_error(ctt_item(data = dat), "at least two items")
})

# ---- ctt_alpha() --------------------------------------------------------------

test_that("ctt_alpha() raw alpha matches the standard variance-based formula", {
  set.seed(4)
  dat <- data.frame(matrix(rbinom(100 * 5, 1, 0.55), nrow = 100))
  out <- ctt_alpha(data = dat)

  k <- ncol(dat)
  item_var <- vapply(dat, stats::var, numeric(1))
  total_var <- stats::var(rowSums(dat))
  alpha_ref <- (k / (k - 1)) * (1 - sum(item_var) / total_var)

  expect_equal(out$alpha, round(alpha_ref, 3))
  expect_equal(out$sem, round(stats::sd(rowSums(dat)) * sqrt(1 - alpha_ref), 3))
})

test_that("ctt_alpha() standardized alpha matches the mean-inter-item-correlation formula", {
  set.seed(5)
  dat <- data.frame(matrix(rbinom(100 * 5, 1, 0.5), nrow = 100))
  out <- ctt_alpha(data = dat)

  k <- ncol(dat)
  item_cor <- stats::cor(dat)
  r_bar <- mean(item_cor[upper.tri(item_cor)])
  alpha_std_ref <- (k * r_bar) / (1 + (k - 1) * r_bar)

  expect_equal(out$alpha_std, round(alpha_std_ref, 3))
})

test_that("ctt_alpha() mean discrimination columns match ctt_item()'s own output", {
  set.seed(6)
  dat <- data.frame(matrix(rbinom(80 * 6, 1, 0.6), nrow = 80))
  item_out <- ctt_item(data = dat, flag = FALSE)$item
  alpha_out <- ctt_alpha(data = dat)

  expect_equal(alpha_out$mean_discrimination_raw,
               round(mean(item_out$discrimination_raw, na.rm = TRUE), 3))
  expect_equal(alpha_out$mean_discrimination_corrected,
               round(mean(item_out$discrimination_corrected, na.rm = TRUE), 3))
})

test_that("ctt_alpha() returns NA alpha when total score has zero variance", {
  dat <- data.frame(I1 = c(1, 1, 1), I2 = c(0, 0, 0), I3 = c(1, 1, 1))
  # total score is constant (2,2,2) for all examinees
  out <- ctt_alpha(data = dat)
  expect_true(is.na(out$alpha))
  expect_true(is.na(out$sem))
})

# ---- ctt(): combined wrapper + print/summary S3 methods ----------------------

test_that("ctt() bundles the same item/alpha content as calling the pieces separately", {
  set.seed(7)
  dat <- data.frame(matrix(rbinom(150 * 8, 1, 0.6), nrow = 150))
  out <- ctt(data = dat)

  expect_s3_class(out, "ctt")
  expect_identical(out$item, ctt_item(data = dat)$item)
  expect_identical(out$alpha, ctt_alpha(data = dat))
})

test_that("ctt()'s internally derived total score matches freq_score()'s tabulated scores", {
  set.seed(8)
  dat <- data.frame(matrix(rbinom(60 * 5, 1, 0.5), nrow = 60))
  out <- ctt(data = dat)
  ref_freq <- freq_score(rowSums(dat))
  expect_equal(out$freq, ref_freq)
})

test_that("print.ctt() and summary.ctt()/print.summary.ctt() run without error and return invisibly", {
  set.seed(9)
  dat <- data.frame(matrix(rbinom(40 * 5, 1, 0.6), nrow = 40))
  out <- ctt(data = dat)

  expect_output(print(out), "Classical Test Theory")
  expect_true(inherits(print(out), "ctt"))

  smry <- summary(out)
  expect_s3_class(smry, "summary.ctt")
  expect_output(print(smry), "Item-Level Statistics")
  expect_output(print(smry), "Test-Level Reliability Summary")
  expect_output(print(smry), "Total-Score Frequency Distribution")
})

test_that("ctt() works on the bundled LSAT6 dataset without error", {
  skip_if_not(exists("LSAT6"), "LSAT6 dataset not available")
  out <- ctt(data = as.data.frame(LSAT6))
  expect_s3_class(out, "ctt")
  expect_true(out$alpha$alpha >= 0 && out$alpha$alpha <= 1)
})

# ---- score and cats validation -----------------------------------------------

test_that("the CTT functions stop when an item score is out of range", {
  dat <- data.frame(I1 = c(0, 1, 1, 0, 1), I2 = c(1, 0, 3, 1, 0),
                    I3 = c(1, 1, 0, 0, 1))

  # a score above cats - 1
  expect_error(ctt_item(data = dat, cats = c(2, 2, 2)), "between 0 and cats - 1")
  expect_error(ctt_alpha(data = dat, cats = c(2, 2, 2)), "V2")
  expect_error(ctt(data = dat, cats = c(2, 2, 2)), "between 0 and cats - 1")

  # a negative score
  dat_neg <- dat
  dat_neg$I2[3] <- -1
  expect_error(ctt_item(data = dat_neg), "V2")
  expect_error(ctt(data = dat_neg), "between 0 and cats - 1")

  # a non-integer score
  dat_frac <- dat
  dat_frac$I2[3] <- 0.5
  expect_error(ctt_item(data = dat_frac), "whole numbers")
  expect_error(ctt_alpha(data = dat_frac), "whole numbers")

  # a character column
  dat_chr <- dat
  dat_chr$I2 <- as.character(dat_chr$I2)
  expect_error(ctt_item(data = dat_chr), "numeric")
  expect_error(ctt(data = dat_chr), "V2")
})

test_that("the CTT functions stop when cats is not a whole number of at least 2", {
  dat <- data.frame(I1 = c(0, 1, 1, 0, 1), I2 = c(1, 0, 1, 1, 0),
                    I3 = c(1, 1, 0, 0, 1))
  expect_error(ctt_item(data = dat, cats = c(1, 1, 1)), "at least 2")
  expect_error(ctt_alpha(data = dat, cats = c(2, 2, 1)), "at least 2")
  expect_error(ctt(data = dat, cats = c(2, 2.5, 2)), "whole numbers")
  expect_error(ctt_item(data = dat, cats = c(2, 2)), "length(cats)",
               fixed = TRUE)
})

test_that("an item that every examinee scores 0 on has two inferred categories", {
  dat <- data.frame(I1 = c(0, 1, 1, 0, 1, 1), I2 = c(0, 0, 0, 0, 0, 0),
                    I3 = c(1, 1, 0, 0, 1, 0))
  out <- ctt_item(data = dat)$item
  expect_equal(out$cats[2], 2)
  expect_equal(out$difficulty[2], 0)
  expect_match(out$flag[2], "difficulty too low")

  # the all-zero item enters the mean difficulty as 0
  alpha_out <- ctt_alpha(data = dat)
  expect_equal(alpha_out$mean_difficulty, round(mean(out$difficulty), 3))
})

# ---- listwise deletion, row names, and constant items ------------------------

test_that("ctt() reports listwise deletion once and matches the result on complete data", {
  set.seed(21)
  dat <- data.frame(matrix(rbinom(120 * 6, 1, 0.6), nrow = 120))
  dat[c(3, 10), 2] <- NA
  dat[7, 5] <- 9

  # one warning that names ctt()
  w <- testthat::capture_warnings(out <- ctt(data = dat, missing = 9))
  expect_length(w, 1)
  expect_match(w, "3 examinee(s) with missing item responses were excluded listwise from ctt().",
               fixed = TRUE)

  # same item, alpha, and frequency results as on the complete rows
  dat_complete <- dat[-c(3, 7, 10), ]
  out_complete <- ctt(data = dat_complete)
  expect_identical(out$item, out_complete$item)
  expect_identical(out$alpha, out_complete$alpha)
  expect_identical(out$freq, out_complete$freq)
})

test_that("ctt() stops when no examinee has complete responses", {
  dat <- data.frame(I1 = c(0, NA, 1), I2 = c(NA, 1, 1), I3 = c(1, 0, NA))
  expect_error(suppressWarnings(ctt(data = dat)), "No examinee has complete")
})

test_that("ctt_item() returns default row names when cats is inferred", {
  set.seed(22)
  dat <- data.frame(matrix(rbinom(80 * 4, 1, 0.5), nrow = 80))
  out <- ctt_item(data = dat)$item
  expect_identical(rownames(out), as.character(1:4))
})

test_that("ctt_alpha() gives no warning for a constant item", {
  dat <- data.frame(I1 = c(0, 1, 1, 0, 1, 1), I2 = c(1, 1, 1, 1, 1, 1),
                    I3 = c(1, 1, 0, 0, 1, 0))
  expect_no_warning(out <- ctt_alpha(data = dat))
  expect_false(is.na(out$alpha_std))
})

# ---- input checks ------------------------------------------------------------

test_that("the CTT functions stop when fewer than two examinees remain", {
  dat <- data.frame(a = 1, b = 0, c = 1)
  expect_error(ctt(data = dat), "two examinees")
  expect_error(ctt_item(data = dat), "two examinees")
  expect_error(ctt_alpha(data = dat), "two examinees")
})

test_that("ctt_item() checks the flagging thresholds only when flagging", {
  set.seed(31)
  dat <- data.frame(matrix(rbinom(30 * 4, 1, 0.5), nrow = 30))
  expect_error(ctt_item(data = dat, crit.p = 0.1), "crit.p")
  expect_error(ctt_item(data = dat, crit.p = c(0.95, 0.1)), "crit.p")
  expect_error(ctt_item(data = dat, crit.p = c(0.1, NA)), "crit.p")
  expect_error(ctt_item(data = dat, crit.dis = NA_real_), "crit.dis")
  expect_error(ctt_item(data = dat, crit.dis = c(0.2, 0.3)), "crit.dis")
  expect_error(ctt(data = dat, crit.p = 0.1), "crit.p")
  expect_no_error(ctt_item(data = dat, flag = FALSE, crit.p = 0.1))
})

# ---- standardized alpha with a constant item ---------------------------------

test_that("ctt_alpha() keeps a constant item in the item count of the standardized alpha", {
  set.seed(3)
  th <- rnorm(200)
  dat <- data.frame(sapply(1:6, function(j) as.integer(th + rnorm(200) > (j - 3) / 2)))
  dat$K <- 1L
  out <- ctt_alpha(data = dat)

  # raw alpha of the standardized items, with the constant item at zero variance
  z <- as.data.frame(scale(dat[, 1:6]))
  z$K <- 0
  alpha_ref <- 7 / 6 * (1 - 6 / stats::var(rowSums(z)))
  expect_equal(out$alpha_std, round(alpha_ref, 3))
  expect_equal(out$alpha_std, 0.674)
  expect_equal(ctt(dat)$alpha$alpha_std, 0.674)
})

test_that("ctt_alpha() standardized alpha is at most 1 and NA when undefined", {
  a <- c(0, 0, 1, 1, 1, 0, 1, 0, 1, 0)
  b <- 1 - a
  b[1] <- 0
  out <- ctt_alpha(data = data.frame(a = a, b = b, c = 1))
  expect_equal(out$alpha, -6.667)
  expect_equal(out$alpha_std, -6.674)

  # a zero denominator of the standardized form
  dat_zero <- data.frame(a = c(1, NA, 0), b = c(0, 1, 1))
  out_zero <- suppressWarnings(ctt(data = dat_zero))
  expect_true(is.na(out_zero$alpha$alpha_std))

  # every item constant: all summaries are NA, not NaN
  out_const <- ctt_alpha(data = data.frame(a = c(1, 1, 1), b = c(0, 0, 0)))
  expect_true(is.na(out_const$alpha_std))
  expect_false(is.nan(out_const$mean_discrimination_raw))
  expect_true(is.na(out_const$mean_discrimination_raw))
  expect_true(is.na(out_const$mean_discrimination_corrected))
})

# ---- scores close to whole numbers -------------------------------------------

test_that("ctt() treats scores within tolerance of whole numbers as whole numbers", {
  dat <- data.frame(a = c(1, 0, 1, 1, 0, 1), b = c(1, 1, 0, 1, 0, 1),
                    c = c(1, 0, 1, 1, 0, 1))
  dat_near <- dat
  dat_near$c[1] <- 1 + 3e-9

  out <- ctt(data = dat_near)
  out_exact <- ctt(data = dat)
  expect_equal(sum(out$freq$freq), nrow(dat))
  expect_equal(out$freq$cum_pct[nrow(out$freq)], 100)
  expect_identical(out$freq, out_exact$freq)
  expect_identical(out$item, out_exact$item)
  expect_identical(out$alpha, out_exact$alpha)
})

# ---- flag for an undefined discrimination ------------------------------------

test_that("ctt_item() flags an item whose discrimination is undefined", {
  set.seed(1)
  dat <- data.frame(matrix(rbinom(40 * 4, 1, 0.5), nrow = 40))
  dat$P <- 1
  out <- ctt_item(data = dat, cats = c(2, 2, 2, 2, 3))$item
  expect_equal(out$difficulty[5], 0.5)
  expect_true(is.na(out$discrimination_raw[5]))
  expect_equal(out$flag[5], "discrimination undefined")
  expect_false(any(grepl("undefined", out$flag[1:4])))

  # the flag does not depend on which correlation is used for flagging
  out_corr <- ctt_item(data = dat, cats = c(2, 2, 2, 2, 3), correct = TRUE)$item
  expect_equal(out_corr$flag[5], "discrimination undefined")

  # no flag column when flagging is off
  expect_false("flag" %in% names(ctt_item(data = dat, flag = FALSE)$item))

  # the undefined item counts in the number of flagged items
  expect_output(print(ctt(data = dat, cats = c(2, 2, 2, 2, 3))),
                "Flagged items: 1 of 5")

  # an extreme item can carry both a difficulty and an undefined flag
  dat_zero <- data.frame(I1 = c(0, 1, 1, 0, 1, 1), I2 = c(0, 0, 0, 0, 0, 0),
                         I3 = c(1, 1, 0, 0, 1, 0))
  expect_equal(ctt_item(data = dat_zero)$item$flag[2],
               "difficulty too low; discrimination undefined")
})

# ---- correlation used for flagging -------------------------------------------

test_that("ctt() records which item-total correlation is used for flagging", {
  set.seed(32)
  dat <- data.frame(matrix(rbinom(50 * 5, 1, 0.6), nrow = 50))

  expect_identical(names(ctt(data = dat)$crit), c("crit.p", "crit.dis", "correct"))
  expect_identical(ctt(data = dat)$crit$correct, FALSE)
  expect_identical(ctt(data = dat, correct = TRUE)$crit$correct, TRUE)
  expect_identical(summary(ctt(data = dat, correct = TRUE))$crit$correct, TRUE)
  expect_identical(ctt_item(data = dat, flag = FALSE)$crit$correct, FALSE)
  expect_identical(ctt_item(data = dat)$crit$crit.dis, 0.2)

  # the full report names the correlation
  expect_output(print(summary(ctt(data = dat))),
                "discrimination >= 0.2 (raw item-total correlation)",
                fixed = TRUE)
  expect_output(print(summary(ctt(data = dat, correct = TRUE))),
                "discrimination >= 0.2 (corrected item-total correlation)",
                fixed = TRUE)

  # an object without crit$correct omits the correlation name
  smry <- summary(ctt(data = dat))
  smry$crit$correct <- NULL
  lines <- utils::capture.output(print(smry))
  expect_true(any(grepl("discrimination >= 0.2$", lines)))
  expect_false(any(grepl("item-total correlation)", lines, fixed = TRUE)))
})

# ---- hand-checked values for a mixed-format data set -------------------------

dat_mixed <- data.frame(
  I1 = c(1, 1, 1, 0, 1, 0, 1, 0, 0, 1),
  I2 = c(1, 1, 0, 1, 1, 0, 0, 0, 0, 1),
  I3 = c(2, 1, 2, 0, 1, 0, 1, 1, 0, 2),
  I4 = c(1, 0, 1, 0, 1, 1, 0, 0, 0, 1)
)

test_that("ctt_item() gives alpha with the item removed and both discrimination columns", {
  out <- ctt_item(data = dat_mixed, flag = FALSE)$item

  # values agree with psych::alpha() (alpha.drop$raw_alpha, raw.r, r.drop)
  expect_equal(out$alpha_removed, c(0.583, 0.786, 0.605, 0.718))
  expect_equal(out$discrimination_raw, c(0.866, 0.574, 0.889, 0.688))
  expect_equal(out$discrimination_corrected, c(0.761, 0.330, 0.696, 0.484))
})

test_that("ctt_alpha() summarizes a mixed-format data set", {
  out <- ctt_alpha(data = dat_mixed)
  expect_equal(out$alpha, 0.746)
  expect_equal(out$alpha_std, 0.751)
  expect_equal(out$sem, 0.927)
  expect_equal(out$mean_difficulty, 0.525)
  expect_equal(out$mean_discrimination_raw, 0.754)
  expect_equal(out$mean_discrimination_corrected, 0.568)
})

test_that("correct selects the correlation used for the discrimination flag", {
  raw_flag <- ctt_item(data = dat_mixed, crit.dis = 0.5)$item$flag
  corrected_flag <- ctt_item(data = dat_mixed, crit.dis = 0.5,
                             correct = TRUE)$item$flag
  expect_equal(raw_flag, rep("", 4))
  expect_equal(corrected_flag,
               c("", "discrimination too low", "", "discrimination too low"))
})

test_that("print.ctt() omits the flagged items line when flag is FALSE", {
  out <- ctt(data = dat_mixed, flag = FALSE)
  lines <- utils::capture.output(print(out))
  expect_false(any(grepl("Flagged items", lines, fixed = TRUE)))
  expect_true(any(grepl("Flagged items: 0 of 4", utils::capture.output(print(ctt(data = dat_mixed))),
                        fixed = TRUE)))
})

test_that("print.summary.ctt() applies digits to the item table", {
  lines <- utils::capture.output(print(summary(ctt(data = dat_mixed)), digits = 2))
  expect_true("    V1     2         0.6                0.87                      0.76" %in% lines)
})

# ---- boundary cases ----------------------------------------------------------

test_that("a score within rounding error of 0 is counted as 0", {
  dat <- data.frame(a = c(1, 0, 1, 0.3 - (0.1 + 0.2)), b = c(1, 1, 0, 0),
                    c = c(0, 1, 1, 1))
  out <- ctt(data = dat)
  expect_equal(out$freq$freq[out$freq$score == 1], 1L)
  expect_equal(sum(out$freq$freq), 4)
  expect_equal(ctt_distr(data = dat)$distr$freq[1], 2L)
  expect_error(ctt(data = data.frame(a = c(1, 0, -1e-6), b = c(1, 1, 0))),
               "whole numbers")
})

test_that("ctt_alpha() gives a standardized alpha of 0 when one item varies", {
  dat <- data.frame(a = c(0, 1, 1, 0, 1), b = 1, c = 1)
  out <- ctt_alpha(data = dat)
  expect_equal(out$alpha, 0)
  expect_equal(out$alpha_std, 0)
  expect_equal(ctt(dat)$alpha$alpha_std, 0)
})
