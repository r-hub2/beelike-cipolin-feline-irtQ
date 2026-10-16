# Simulate a small dichotomous data set with missing responses, an examinee who responded only to
# item 3, and an examinee without any response
dc_sim <- function(seed = 8) {
  set.seed(seed)
  x <- shape_df(
    par.drm = list(a = runif(8, 0.8, 2), b = rnorm(8), g = c(rep(0, 4), runif(4, 0.05, 0.25))),
    cats = 2, model = c(rep("2PLM", 4), rep("3PLM", 4))
  )
  resp <- simdat(x, theta = rnorm(300), D = 1)
  resp[sample(length(resp), 0.3 * length(resp))] <- NA
  resp[1, ] <- NA
  resp[1, 3] <- 1
  resp[2, ] <- NA
  list(x = x, resp = resp)
}

test_that("dif_rescore() rescoring only the examinees who responded to the removed item matches a full rescoring", {
  skip_on_cran()
  sim <- dc_sim()
  for (m in c("ML", "MLF", "WL", "MAP", "EAP")) {
    score <- suppressWarnings(est_score(sim$x, sim$resp, D = 1, method = m)$est.theta)

    # remove item 3 and rescore
    x_puri <- sim$x[-3, ]
    data_puri <- sim$resp[, -3]
    loc_resp <- which(!is.na(sim$resp[, 3]))
    part <- dif_rescore(
      x = x_puri, data = data_puri, score = score, loc_resp = loc_resp, first = FALSE, D = 1,
      method = m, range = c(-5, 5), norm.prior = c(0, 1), nquad = 41, weights = NULL, ncore = 1,
      quiet_na = TRUE
    )
    full <- suppressWarnings(est_score(x_puri, data_puri, D = 1, method = m)$est.theta)
    expect_identical(part, full)

    # the examinee who responded only to the removed item has no ability estimate left
    expect_true(is.na(part[1]))
    expect_true(is.na(part[2]))
  }
})

test_that("dif_rescore() rescores all examinees in the first iteration and mutes only the requested warnings", {
  sim <- dc_sim()
  score <- rep(0, nrow(sim$resp))
  args <- list(
    x = sim$x, data = sim$resp, score = score, loc_resp = 1L, first = TRUE, D = 1, method = "ML",
    range = c(-5, 5), norm.prior = c(0, 1), nquad = 41, weights = NULL, ncore = 1
  )

  # all examinees are rescored when first = TRUE, whatever loc_resp is
  full <- suppressWarnings(est_score(sim$x, sim$resp, D = 1)$est.theta)
  expect_identical(suppressWarnings(do.call(dif_rescore, args)), full)

  # the warning about the examinees without a response is muted only when requested
  expect_warning(do.call(dif_rescore, args), "NA values are returned")
  expect_no_warning(do.call(dif_rescore, c(args, list(quiet_na = TRUE))))
})

test_that("check_crit_val() stops for invalid values", {
  expect_null(check_crit_val(NULL))
  expect_no_error(check_crit_val(0))
  expect_no_error(check_crit_val(0.002))
  for (bad in list("a", -1, c(1, 2), NA, Inf)) {
    expect_error(check_crit_val(bad), "crit.val")
  }
})

test_that("check_min_resp(), check_score(), and check_pattern_method() stop for invalid values", {
  expect_null(check_min_resp(NULL))
  expect_no_error(check_min_resp(0))
  expect_no_error(check_min_resp(7.5))
  for (bad in list("a", -1, c(1, 2), NA, Inf)) {
    expect_error(check_min_resp(bad), "min.resp")
  }

  expect_identical(check_score(c(0.5, NA, 1), 3), c(0.5, NA, 1))
  expect_identical(check_score(matrix(c(0.5, 1), ncol = 1), 2), c(0.5, 1))
  expect_error(check_score(c("a", "b"), 2), "numeric")
  expect_error(check_score(c(0, 1), 3), "length")

  expect_no_error(check_pattern_method("ML", fn = "ripd()"))
  for (bad in list("EAP.SUM", "INV.TCC", c("ML", "EAP"), NULL)) {
    expect_error(check_pattern_method(bad, fn = "ripd()"), "not supported by ripd()", fixed = TRUE)
  }
  expect_error(check_pattern_method("INV.TCC", fn = "rdif()"), "not supported by rdif()",
               fixed = TRUE)
})
