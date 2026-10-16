# ctt_distr() computes option/category response distributions and
# point-biserial correlations, in both selected-response (with a key) and
# scored-category (already-scored, no key) modes.

test_that("ctt_distr() selected-response mode reproduces a hand-worked distribution", {
  raw <- data.frame(
    V1 = c("1", "2", "1", "1", "3", "2", "1", "4"),
    V2 = c("4", "4", "2", "2,4", "4", "1", "4", "4")
  )
  key <- c(1, 4)
  out <- ctt_distr(data = raw, key = key)

  # V1: option 1 chosen by rows 1,3,4,7 -> freq 4 of 8 -> 50%
  v1_opt1 <- out$distr[out$distr$item == "V1" & out$distr$option == 1, ]
  expect_equal(v1_opt1$freq, 4L)
  expect_equal(v1_opt1$pct, 50)
  expect_true(v1_opt1$is_key)

  # V2: option 4 is the key; row 4 ("2,4") is double-marked and must not count
  # as choosing either option 2 or option 4
  v2_opt4 <- out$distr[out$distr$item == "V2" & out$distr$option == 4, ]
  expect_equal(v2_opt4$freq, 5L)  # rows 1,2,5,7,8
  expect_false(v2_opt4$freq == 6L)

  # omission/double-marking: V2 has exactly one double-marked response
  expect_equal(out$omit$pct_double[out$omit$item == "V2"], round(100 * 1 / 8, 2))
  expect_equal(out$omit$pct_blank[out$omit$item == "V1"], 0)
})

test_that("ctt_distr() flags a distractor whose point-biserial correlation is positive", {
  # 10 examinees ranked by an externally supplied total score; the wrong
  # option "2" is chosen only by the 5 highest-scoring examinees, and the
  # keyed option "1" only by the 5 lowest-scoring examinees, so option 2's
  # point-biserial correlation with total must be strongly positive.
  raw <- data.frame(target = c("2", "2", "2", "2", "2", "1", "1", "1", "1", "1"))
  total <- c(10, 9, 8, 7, 6, 5, 4, 3, 2, 1)
  out <- ctt_distr(data = raw, key = c(1), total = total)
  distractor_row <- out$distr[out$distr$option == 2, ]
  expect_gt(distractor_row$pb_raw, 0)
  expect_match(distractor_row$flag, "more attractive")
})

test_that("ctt_distr() scored-category mode reports category frequencies with no key concept", {
  scored <- data.frame(
    I1 = c(0, 1, 2, 1, 0, 2, 1, 1),
    I2 = c(1, 1, 0, 2, 2, 1, 0, 1)
  )
  out <- ctt_distr(data = scored)
  expect_null(out$omit)
  expect_true(all(is.na(out$distr$is_key)))

  i1_cat0 <- out$distr[out$distr$item == "V1" & out$distr$option == 0, ]
  expect_equal(i1_cat0$freq, sum(scored$I1 == 0))
})

test_that("ctt_distr() errors when a supplied total score has the wrong length", {
  raw <- data.frame(V1 = c("1", "2", "1"))
  expect_error(ctt_distr(data = raw, key = c(1), total = c(1, 2)),
               "length\\(total\\) must equal nrow\\(data\\)")
})

test_that("ctt_distr() pb_raw/pb_corrected are NA for a constant (all-key) option", {
  raw <- data.frame(V1 = c("1", "1", "1", "1"))
  out <- ctt_distr(data = raw, key = c(1))
  row <- out$distr[out$distr$item == "V1" & out$distr$option == 1, ]
  expect_true(is.na(row$pb_raw))
  expect_true(is.na(row$pb_corrected))
})

test_that("ctt_distr() scores a letter-coded (A-E) test form identically to its numeric analog", {
  # same responses/keys as the hand-worked numeric example above, but
  # expressed with letter option labels (1->A, 2->B, 3->C, 4->D); every
  # numeric output column should match row-for-row once options are mapped
  raw_num <- data.frame(
    V1 = c("1", "2", "1", "1", "3", "2", "1", "4"),
    V2 = c("4", "4", "2", "2,4", "4", "1", "4", "4")
  )
  out_num <- ctt_distr(data = raw_num, key = c(1, 4))

  raw_let <- data.frame(
    V1 = c("A", "B", "A", "A", "C", "B", "A", "D"),
    V2 = c("D", "D", "B", "B,D", "D", "A", "D", "D")
  )
  out_let <- ctt_distr(data = raw_let, key = c("A", "D"))

  option_map <- c("A" = 1, "B" = 2, "C" = 3, "D" = 4)
  out_let$distr$option_num <- option_map[out_let$distr$option]
  merged <- merge(out_num$distr, out_let$distr,
                   by.x = c("item", "option"), by.y = c("item", "option_num"),
                   suffixes = c("_num", "_let"))

  expect_equal(nrow(merged), nrow(out_num$distr))
  expect_equal(merged$freq_num, merged$freq_let)
  expect_equal(merged$pct_num, merged$pct_let)
  expect_equal(merged$pb_raw_num, merged$pb_raw_let)
  expect_equal(merged$pb_corrected_num, merged$pb_corrected_let)
  expect_equal(merged$is_key_num, merged$is_key_let)
  expect_equal(merged$flag_num, merged$flag_let)
  expect_equal(out_num$omit$pct_blank, out_let$omit$pct_blank)
  expect_equal(out_num$omit$pct_double, out_let$omit$pct_double)
})

test_that("ctt_distr() matches letter-coded options case-insensitively", {
  raw <- data.frame(V1 = c("a", "A", "b", "B", "a", "B", "A", "a"))
  out <- ctt_distr(data = raw, key = c("A"))
  row_a <- out$distr[out$distr$option == "A", ]
  expect_equal(row_a$freq, 5L)
  expect_true(row_a$is_key)
})

test_that("ctt_distr() scores a general-label (non-Latin) test form correctly", {
  # option labels are Korean syllables (romanized ga/na/da), written with
  # \uXXXX escapes for ASCII source portability, as in score_resp()
  raw <- data.frame(V1 = c("\uAC00", "\uB098", "\uB2E4", "\uAC00",
                            "\uB098", "\uB2E4", "\uAC00", "\uAC00"))
  out <- ctt_distr(data = raw, key = c("\uAC00"))
  row_key <- out$distr[out$distr$option == "\uAC00", ]
  expect_equal(row_key$freq, 4L)
  expect_true(row_key$is_key)
})

test_that("ctt_distr() respects an explicit character opt= for a letter-coded item", {
  raw <- data.frame(
    V1 = c("A", "B", "A", "A", "C", "B", "A", "D"),
    V2 = c("D", "D", "B", "B,D", "D", "A", "D", "D")
  )
  out <- ctt_distr(data = raw, key = c("A", "D"), opt = c("A", "B", "C", "D"))
  expect_equal(sort(unique(out$distr$option[out$distr$item == "V1"])),
               c("A", "B", "C", "D"))
})

test_that("ctt_distr() matches a key data frame by item number for a character item column", {
  # twelve items so that string order differs from numeric order
  set.seed(12)
  n_item <- 12
  raw <- as.data.frame(
    matrix(sample(as.character(1:4), 30 * n_item, replace = TRUE), 30, n_item)
  )
  key_vec <- rep(c(1, 4, 3, 2), length.out = n_item)
  out_vec <- ctt_distr(data = raw, key = key_vec)

  # the same key rows in a shuffled order, with character and factor item
  # columns
  shuffled <- c(5, 12, 1, 9, 3, 11, 7, 2, 10, 6, 8, 4)
  key_chr <- data.frame(item = as.character(shuffled), key = key_vec[shuffled])
  key_fct <- data.frame(item = factor(shuffled), key = key_vec[shuffled])
  expect_identical(ctt_distr(data = raw, key = key_chr)$distr$is_key,
                   out_vec$distr$is_key)
  expect_identical(ctt_distr(data = raw, key = key_fct)$distr$is_key,
                   out_vec$distr$is_key)
})

test_that("ctt_distr() scored-category mode stops when a score is out of range", {
  dat <- data.frame(I1 = c(0, 1, 1, 0, 1), I2 = c(1, 0, 3, 1, 0))
  expect_error(ctt_distr(data = dat, cats = c(2, 2)), "between 0 and cats - 1")
  dat$I2[3] <- -1
  expect_error(ctt_distr(data = dat), "V2")
  dat$I2[3] <- 0.5
  expect_error(ctt_distr(data = dat), "whole numbers")
  expect_error(ctt_distr(data = data.frame(I1 = c(0, 1, 1), I2 = c(1, 0, 1)),
                         cats = c(2, 1)), "at least 2")
})

test_that("ctt_distr() scored-category mode gives an all-zero item two categories", {
  dat <- data.frame(I1 = c(0, 1, 1, 0, 1, 1), I2 = c(0, 0, 0, 0, 0, 0))
  out <- ctt_distr(data = dat)$distr
  out_i2 <- out[out$item == "V2", ]
  expect_equal(out_i2$option, c(0, 1))
  expect_equal(out_i2$freq, c(6L, 0L))
})

# ---- input checks ------------------------------------------------------------

test_that("ctt_distr() stops when numeric and non-numeric keys are mixed", {
  raw <- data.frame(V1 = c("1", "2", "1", "3"), V2 = c("A", "B", "A", "C"))
  expect_error(ctt_distr(data = raw, key = c(1, "A")), "numeric for every item")
})

test_that("ctt_distr() stops when total is not numeric or has missing values", {
  raw <- data.frame(V1 = c("1", "2", "1", "3"))
  expect_error(ctt_distr(data = raw, key = 1, total = c(3, NA, 2, 1)),
               "no missing values")
  expect_error(ctt_distr(data = raw, key = 1, total = c("3", "2", "2", "1")),
               "numeric")

  dat <- data.frame(a = c(0, 1, 1, 0), b = c(1, 1, 0, 0))
  expect_error(ctt_distr(data = dat, total = c(1, NA, 1, 0)),
               "no missing values")

  # a missing total for an examinee deleted listwise is allowed
  dat_na <- data.frame(a = c(0, 1, 1, 0, NA), b = c(1, 1, 0, 0, 1))
  expect_warning(
    out <- ctt_distr(data = dat_na, total = c(1, 2, 1, 0, NA)),
    "excluded listwise"
  )
  expect_s3_class(out$distr, "data.frame")
})

test_that("ctt_distr() stops when fewer than two examinees are analyzed", {
  expect_error(ctt_distr(data = data.frame(V1 = "1", V2 = "2"), key = c(1, 2)),
               "two examinees")
  expect_error(ctt_distr(data = data.frame(a = 1, b = 0, c = 1)),
               "two examinees")
})

test_that("ctt_distr() checks crit.distractor in selected-response mode", {
  raw <- data.frame(V1 = c("1", "2", "1", "3"))
  expect_error(ctt_distr(data = raw, key = 1, crit.distractor = NA_real_),
               "crit.distractor")
  expect_error(ctt_distr(data = raw, key = 1, crit.distractor = c(0, 0.1)),
               "crit.distractor")
})

test_that("ctt_distr() scored-category mode counts near-integer scores in their category", {
  dat <- data.frame(I1 = c(0, 1, 2, (0.1 + 0.2) * 10 - 1, 1, 1),
                    I2 = c(1, 0, 1, 1, 0, 1))
  out <- ctt_distr(data = dat)$distr
  out_i1 <- out[out$item == "V1", ]
  expect_equal(sum(out_i1$freq), 6)
  expect_equal(out_i1$freq, c(1L, 3L, 2L))
})

test_that("ctt_distr() lists an option written as 1 and 01 once", {
  raw <- data.frame(V1 = c("1", "01", "2", "1", "3"),
                    V2 = c("2", "2", "1", "2", "2"))
  out <- ctt_distr(data = raw, key = c(1, 2))$distr
  expect_equal(nrow(out), 6)
  expect_equal(out$option[out$item == "V1"], c(1, 2, 3))
  expect_equal(out$freq[out$item == "V1"], c(3L, 1L, 1L))
  expect_equal(out$pct[out$item == "V1"], c(60, 20, 20))
  expect_equal(sum(out$pct[out$item == "V2"]), 100)
})

test_that("ctt_distr() treats a missing code with surrounding spaces as an omission", {
  raw <- data.frame(V1 = c("1", " 9", "9 ", "2"), V2 = c("1", "2", "1", "2"))
  out <- ctt_distr(data = raw, key = c(1, 2), missing = "9")
  expect_equal(out$omit$pct_blank, c(50, 0))
  expect_equal(out$distr$option[out$distr$item == "V1"], c(1, 2))
  expect_equal(out$distr$freq[out$distr$item == "V1"], c(1L, 1L))
})

# ---- hand-checked values for the selected-response mode ----------------------

raw_hand <- data.frame(
  V1 = c("1", "2", "1", "1", "3", "2", "1", "4"),
  V2 = c("4", "4", "2", "2,4", "4", "1", "4", "4")
)

test_that("ctt_distr() selected-response mode reproduces the point-biserial correlations", {
  out <- ctt_distr(data = raw_hand, key = c(1, 4))$distr
  pick <- function(item, option) out[out$item == item & out$option == option, ]

  expect_equal(pick("V1", 1)$pb_raw, 0.626)
  expect_equal(pick("V1", 1)$pb_corrected, -0.258)
  expect_equal(pick("V1", 2)$pb_raw, -0.602)
  expect_equal(pick("V1", 2)$pb_corrected, -0.149)
  expect_equal(pick("V2", 4)$pb_raw, 0.592)
  expect_equal(pick("V2", 4)$pb_corrected, -0.258)

  # an option chosen by no one has an undefined correlation
  expect_true(is.na(pick("V2", 3)$pb_raw))
  expect_true(is.na(pick("V2", 3)$pb_corrected))
})

test_that("correct selects the correlation used for the distractor flag", {
  flag_raw <- ctt_distr(data = raw_hand, key = c(1, 4))$distr$flag
  flag_corrected <- ctt_distr(data = raw_hand, key = c(1, 4),
                              correct = TRUE)$distr$flag
  expect_equal(flag_raw, rep("", 8))
  expect_equal(which(flag_corrected != ""), c(3L, 4L, 6L))
  expect_equal(unique(flag_corrected[flag_corrected != ""]),
               "distractor more attractive to high scorers")
})

test_that("ctt_distr() counts a missing code as an omission in selected-response mode", {
  raw9 <- raw_hand
  raw9$V1[c(2, 5)] <- "9"
  out <- ctt_distr(data = raw9, key = c(1, 4), missing = "9")
  expect_equal(out$omit$pct_blank, c(25, 0))
  v1 <- out$distr[out$distr$item == "V1", ]
  expect_equal(v1$option, c(1, 2, 4))
  expect_equal(v1$freq, c(4L, 1L, 1L))
})

test_that("ctt_distr() scored-category mode aligns total with the retained rows", {
  dat <- data.frame(I1 = c(0, 1, 2, 1, -9, 2), I2 = c(1, 0, 1, 1, 0, 1))
  expect_warning(
    out <- ctt_distr(data = dat, total = c(5, 3, 6, 4, 1, 7), missing = -9),
    "1 examinee(s) with missing item responses were excluded listwise",
    fixed = TRUE
  )
  row <- out$distr[out$distr$item == "V1" & out$distr$option == 2, ]
  expect_equal(row$freq, 2L)
  expect_equal(row$pct, 40)
  expect_equal(row$pb_raw, 0.866)
  expect_equal(row$pb_corrected, 0.490)
})

test_that("ctt_distr() lists a supplied option written as 1 and 01 once, in the supplied order", {
  raw <- data.frame(V1 = c("1", "01", "2", "3", "1"),
                    V2 = c("2", "2", "1", "3", "2"))
  out <- ctt_distr(data = raw, key = c(1, 2), opt = c("3", "1", "01", "2"))$distr
  v1 <- out[out$item == "V1", ]
  expect_equal(v1$option, c(3, 1, 2))
  expect_equal(v1$freq, c(1L, 3L, 1L))
  expect_equal(sum(v1$pct), 100)
})
