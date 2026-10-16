# irtfit() computes observed proportions, residuals, and fit statistics for each item

prm_file <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
x_mix <- bring.flexmirt(file = prm_file, "par")$Group1$full_df[c(1:15, 51:55), ]
set.seed(2026)
theta_fit <- rnorm(1000)
resp_fit <- simdat(x = x_mix, theta = theta_fit, D = 1)

test_that("irtfit() observed proportions are the category frequencies divided by the group size", {
  skip_on_cran()
  for (gm in c("equal.width", "equal.freq")) {
    fit <- irtfit(
      x = x_mix, score = theta_fit, data = resp_fit, group.method = gm,
      n.width = 8, loc.theta = "average", range.score = c(-4, 4), D = 1
    )
    for (i in seq_along(fit$contingency.plot)) {
      tb <- fit$contingency.plot[[i]]
      freq <- tb[, grep("^obs\\.freq\\.[0-9]+$", names(tb)), drop = FALSE]
      prop <- tb[, grep("^obs\\.prop\\.[0-9]+$", names(tb)), drop = FALSE]
      # proportions in each interval sum to one
      expect_equal(unname(rowSums(prop)), rep(1, nrow(tb)))
      # proportions equal the frequencies over the interval totals
      expect_equal(unname(as.matrix(prop)), unname(as.matrix(freq)) / tb$total)
    }
  }
})

test_that("irtfit() raw residuals equal observed proportions minus expected probabilities", {
  fit <- irtfit(
    x = x_mix, score = theta_fit, data = resp_fit, group.method = "equal.freq",
    n.width = 8, loc.theta = "average", range.score = c(-4, 4), D = 1
  )
  tb <- fit$contingency.plot[[1]]
  prop <- tb[, grep("^obs\\.prop\\.[0-9]+$", names(tb))]
  expp <- tb[, grep("^exp\\.prob\\.[0-9]+$", names(tb))]
  rsd <- tb[, grep("^raw\\.rsd\\.[0-9]+$", names(tb))]
  expect_equal(unname(as.matrix(rsd)), unname(as.matrix(prop - expp)))
})

test_that("irtfit() drops empty score groups and keeps proportions that sum to one", {
  # two extreme scores leave empty intervals in the equal-width grouping
  sc <- theta_fit
  sc[1] <- -6
  sc[2] <- 6
  fit <- irtfit(
    x = x_mix, score = sc, data = resp_fit, group.method = "equal.width",
    n.width = 10, loc.theta = "average", range.score = c(-7, 7), D = 1
  )
  n_group <- vapply(fit$contingency.plot, nrow, integer(1))
  expect_true(any(n_group < 10L))
  for (tb in fit$contingency.plot) {
    prop <- tb[, grep("^obs\\.prop\\.[0-9]+$", names(tb))]
    expect_equal(unname(rowSums(prop)), rep(1, nrow(tb)))
  }
})

# reads the most recent ggplot with the function that the installed ggplot2 provides
last_ggplot <- function() {
  getter <- if ("get_last_plot" %in% getNamespaceExports("ggplot2")) "get_last_plot" else "last_plot"
  getExportedValue("ggplot2", getter)()
}

test_that("plot.irtfit() Wald intervals use the two-sided critical value", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  fit <- irtfit(
    x = x_mix, score = theta_fit, data = resp_fit, group.method = "equal.freq",
    n.width = 8, loc.theta = "average", range.score = c(-4, 4), D = 1, alpha = 0.05
  )
  plot(x = fit, item.loc = 1, type = "icc", ci.method = "wald", show.table = FALSE)
  built <- ggplot2::ggplot_build(last_ggplot())
  # the segment layer holds the interval (y = upper limit, yend = lower limit)
  seg <- Filter(function(d) all(c("y", "yend", "xend") %in% names(d)), built$data)[[1]]
  unclipped <- seg$yend > 0 & seg$y < 1
  half_width <- (seg$y - seg$yend)[unclipped] / 2
  tb <- fit$contingency.plot[[1]]
  se_all <- c(tb$se.0, tb$se.1)
  # each half width equals the 97.5th percentile of the normal times a standard error
  expect_true(length(half_width) > 0)
  expect_true(all(vapply(half_width, function(h) any(abs(h / stats::qnorm(0.975) - se_all) < 1e-6), logical(1))))
})

# ---- Input checks ---------------------------------------------------------------

x_ft <- x_mix[c(1:3, 20), ]
resp_ft <- resp_fit[, c(1:3, 20)]

test_that("irtfit() stops for data that do not match the items or the score categories", {
  # a column is missing
  expect_error(irtfit(x_ft, score = theta_fit, data = resp_ft[, 1:3], D = 1), "number of columns")

  # responses above the highest category and non-integer responses
  r <- resp_ft
  r[1, 1] <- 2
  expect_error(irtfit(x_ft, score = theta_fit, data = r, D = 1), "outside the score categories")
  r <- resp_ft
  r[2, 4] <- 0.5
  expect_error(irtfit(x_ft, score = theta_fit, data = r, D = 1), "outside the score categories")

  # a code declared with 'missing' is accepted and gives the result for NA
  r <- resp_ft
  r[1, 1] <- -9
  expect_error(irtfit(x_ft, score = theta_fit, data = r, D = 1), "outside the score categories")
  expect_identical(
    irtfit(x_ft, score = theta_fit, data = r, D = 1, missing = -9)$fit_stat,
    irtfit(x_ft, score = theta_fit, data = replace(r, r == -9, NA), D = 1)$fit_stat
  )
})

test_that("irtfit() reads character and factor responses like numeric responses", {
  r_na <- resp_ft
  r_na[1, 2] <- NA
  r_chr <- matrix(as.character(r_na), nrow(r_na))
  r_chr[is.na(r_na)] <- "."
  df_chr <- as.data.frame(r_chr, stringsAsFactors = FALSE)
  df_fac <- as.data.frame(lapply(df_chr, factor))
  ref <- irtfit(x_ft, score = theta_fit, data = r_na, D = 1)
  for (dat in list(r_chr, df_chr, df_fac)) {
    fit <- irtfit(x_ft, score = theta_fit, data = dat, D = 1, missing = ".")
    expect_identical(fit[names(fit) != "call"], ref[names(ref) != "call"])
  }

  # a response that cannot be read as a number stops
  r_chr[1, 1] <- "a"
  expect_error(irtfit(x_ft, score = theta_fit, data = r_chr, D = 1, missing = "."), "numeric scores")
})

test_that("irtfit() drops examinees with missing ability estimates", {
  sc <- theta_fit
  sc[c(5, 10)] <- NA
  fit <- irtfit(x_ft, score = sc, data = resp_ft, D = 1)
  expect_equal(fit$fit_stat$N, rep(998, 4))
  ref <- irtfit(x_ft, score = theta_fit[-c(5, 10)], data = resp_ft[-c(5, 10), ], D = 1)
  expect_equal(fit$fit_stat, ref$fit_stat)
})

test_that("irtfit() keeps working when only one item remains after excluding items", {
  resp2 <- resp_ft[, 1:2]
  resp2[, 2] <- NA
  expect_warning(
    fit <- irtfit(x_ft[1:2, ], score = theta_fit, data = resp2, D = 1),
    "fewer than two responses.*V2|fewer than two responses.*column 2"
  )
  expect_equal(fit$fit_stat$id, x_ft$id[1])
  expect_identical(fit$fit_stat, irtfit(x_ft[1, ], score = theta_fit, data = resp_ft[, 1], D = 1)$fit_stat)

  # the warning names the id and the column of the excluded item
  expect_warning(irtfit(x_ft[1:2, ], score = theta_fit, data = resp2, D = 1), paste0(x_ft$id[2], " [(]column 2[)]"))

  # every item without responses stops
  resp3 <- resp_ft[, 1:2]
  resp3[] <- NA
  expect_error(irtfit(x_ft[1:2, ], score = theta_fit, data = resp3, D = 1), "fewer than two responses")
})

test_that("irtfit() stops for invalid grouping and score arguments", {
  # a misspelled location of the theta point
  expect_error(irtfit(x_ft, score = theta_fit, data = resp_ft, D = 1, loc.theta = "midle"), "must be either")
  expect_identical(
    irtfit(x_ft, score = theta_fit, data = resp_ft, D = 1, loc.theta = "MIDDLE")$fit_stat,
    irtfit(x_ft, score = theta_fit, data = resp_ft, D = 1, loc.theta = "middle")$fit_stat
  )

  # scores that do not match the rows of the data
  expect_error(irtfit(x_ft, score = theta_fit[1:400], data = resp_ft, D = 1), "length of 'score'")

  # scores without any variation
  expect_error(irtfit(x_ft, score = rep(0.5, 1000), data = resp_ft, D = 1), "no variation")

  # a reversed or incomplete score range
  expect_error(irtfit(x_ft, score = theta_fit, data = resp_ft, D = 1, range.score = c(3, -3)), "range.score")
  expect_error(irtfit(x_ft, score = theta_fit, data = resp_ft, D = 1, range.score = 3), "range.score")
})

test_that("plot.irtfit() stops for an unknown type or an item location outside the evaluated items", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  fit <- irtfit(x_ft, score = theta_fit, data = resp_ft, D = 1)
  expect_error(plot(fit, item.loc = 1, type = "bothh"), "must be one of")
  expect_error(plot(fit, item.loc = 9), "item.loc")
  expect_error(plot(fit, item.loc = 0), "item.loc")
  expect_error(plot(fit, item.loc = 1.5), "item.loc")
  expect_no_error(plot(fit, item.loc = 4, type = "SR", show.table = FALSE))
})

test_that("irtfit() reports NA critical values and p-values when no degrees of freedom remain", {
  expect_warning(
    fit <- irtfit(x = x_ft, score = theta_fit[1:30], data = resp_ft[1:30, ], n.width = 10, D = 1),
    "No degrees of freedom"
  )
  no_df <- fit$fit_stat$df.X2 <= 0
  expect_true(any(no_df))
  expect_true(all(is.na(fit$fit_stat$p.X2[no_df])))
  expect_true(all(is.na(fit$fit_stat$crit.val.X2[no_df])))

  # the G2 statistic and the items with degrees of freedom are not affected
  expect_false(anyNA(fit$fit_stat$p.G2))
  expect_false(anyNA(fit$fit_stat$p.X2[!no_df]))

  # no warning when every item keeps degrees of freedom
  expect_no_warning(irtfit(x = x_ft, score = theta_fit, data = resp_ft, D = 1))
})

test_that("plot.irtfit() uses xlab.text for type = 'both'", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  fit <- irtfit(x = x_ft, score = theta_fit, data = resp_ft, D = 1)
  for (tp in c("both", "icc", "sr")) {
    plot(fit, item.loc = 1, type = tp, show.table = FALSE, xlab.text = "Ability")
    expect_identical(last_ggplot()$labels$x, "Ability")
  }

  # the default label is theta
  plot(fit, item.loc = 1, type = "both", show.table = FALSE)
  expect_false(identical(last_ggplot()$labels$x, "Ability"))
})

# ---- Statistics against direct computations ---------------------------------------

x_chk <- x_mix[c(1:3, 20), ]
set.seed(2027)
theta_chk <- rnorm(1500)
resp_chk <- simdat(x = x_chk, theta = theta_chk, D = 1)

test_that("irtfit() X2, G2, and df agree with the uncollapsed contingency tables", {
  fit <- irtfit(x = x_chk, score = theta_chk, data = resp_chk, n.width = 8, min.collapse = 0, D = 1)
  npar <- c(3, 3, 3, 5)
  for (i in 1:4) {
    tb <- fit$contingency.plot[[i]]
    k <- x_chk$cats[i]
    o <- data.matrix(tb[, paste0("obs.freq.", 0:(k - 1))])
    p <- data.matrix(tb[, paste0("exp.prob.", 0:(k - 1))])
    e <- p * tb$total
    expect_equal(fit$fit_stat$X2[i], round(sum((o - e)^2 / e), 3))
    expect_equal(fit$fit_stat$G2[i], round(2 * sum(ifelse(o > 0, o * log(o / e), 0)), 3))
    expect_equal(fit$fit_stat$df.X2[i], nrow(tb) * (k - 1) - npar[i])
    expect_equal(fit$fit_stat$df.G2[i], nrow(tb) * (k - 1))
  }
})

test_that("irtfit() infit and outfit are the weighted and unweighted mean squares", {
  fit <- irtfit(x = x_chk, score = theta_chk, data = resp_chk, D = 1)
  for (i in 1:4) {
    p <- traceline(x = x_chk[i, ], theta = theta_chk, D = 1)$prob.cats[[1]]
    k <- 0:(x_chk$cats[i] - 1)
    ex <- as.vector(p %*% k)
    v <- as.vector(p %*% k^2) - ex^2
    z <- resp_chk[, i] - ex
    expect_equal(fit$fit_stat$outfit[i], round(mean(z^2 / v), 3))
    expect_equal(fit$fit_stat$infit[i], round(sum(z^2) / sum(v), 3))
  }
})

test_that("irtfit() evaluates the middle of each interval when loc.theta = 'middle'", {
  fit <- irtfit(x = x_chk, score = theta_chk, data = resp_chk, n.width = 10, loc.theta = "middle", D = 1)
  cuts <- seq(min(theta_chk), max(theta_chk), length.out = 11)
  expect_equal(fit$contingency.plot[[1]]$point, (cuts[-1] + cuts[-11]) / 2)
})

test_that("irtfit() for est_item objects uses the scores, data, and scaling constant of the object", {
  mod <- suppressWarnings(
    est_item(x = x_chk, data = resp_chk, score = theta_chk, D = 1.702, use.gprior = TRUE, verbose = FALSE)
  )
  f1 <- irtfit(mod)
  f2 <- irtfit(mod$par.est, score = theta_chk, data = resp_chk, D = 1.702)
  expect_identical(f1$fit_stat, f2$fit_stat)
  expect_equal(f1$ancillary$scale.D, 1.702)
})

test_that("irtfit() for est_irt objects uses the data and scaling constant of the object", {
  skip_on_cran()
  mod <- est_irt(x = x_chk, data = resp_chk, D = 1.702, use.gprior = TRUE, verbose = FALSE)
  f1 <- irtfit(mod, score = theta_chk)
  f2 <- irtfit(mod$par.est, score = theta_chk, data = resp_chk, D = 1.702)
  expect_identical(f1$fit_stat, f2$fit_stat)
  expect_equal(f1$ancillary$scale.D, 1.702)
})

test_that("irtfit() accepts only the defined values of loc.theta, ignoring case", {
  for (bad in c("mid", "m", "avg", "midle", "")) {
    expect_error(
      irtfit(x_ft, score = theta_fit, data = resp_ft, D = 1, loc.theta = bad),
      "must be either"
    )
  }
  expect_error(
    irtfit(x_ft, score = theta_fit, data = resp_ft, D = 1, loc.theta = c("average", "middle")),
    "must be either"
  )
  expect_identical(
    irtfit(x_ft, score = theta_fit, data = resp_ft, D = 1, loc.theta = "Average")$fit_stat,
    irtfit(x_ft, score = theta_fit, data = resp_ft, D = 1)$fit_stat
  )
})

test_that("plot.irtfit() accepts only the defined values of type, ignoring case", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  fit <- irtfit(x_ft, score = theta_fit, data = resp_ft, D = 1)
  for (bad in c("bo", "i", "s", "")) {
    expect_error(plot(fit, item.loc = 1, type = bad), "must be one of")
  }
  expect_no_error(plot(fit, item.loc = 1, type = "ICC", show.table = FALSE))
})

test_that("irtfit() names the items and columns of responses outside the score categories", {
  r <- resp_ft
  r[1, 2] <- 2
  expect_error(irtfit(x_ft, score = theta_fit, data = r, D = 1), paste0(x_ft$id[2], " [(]column 2[)]"))
})

test_that("irtfit() reads logical responses and a single remaining group", {
  # a logical response matrix gives the result for the 0/1 matrix
  x_dich <- x_mix[1:3, ]
  r01 <- resp_fit[, 1:3]
  rlg <- r01 == 1
  expect_identical(
    irtfit(x_dich, score = theta_fit, data = rlg, D = 1)[c("fit_stat", "contingency.fitstat", "individual.info")],
    irtfit(x_dich, score = theta_fit, data = r01, D = 1)[c("fit_stat", "contingency.fitstat", "individual.info")]
  )

  # every ability group is merged into one row
  expect_warning(
    fit <- irtfit(x_ft, score = theta_fit[1:60], data = resp_ft[1:60, ], D = 1, min.collapse = 50),
    "No degrees of freedom"
  )
  expect_true(all(vapply(fit$contingency.fitstat, nrow, integer(1)) == 1L))
  expect_true(all(fit$fit_stat$df.X2 <= 0))
  expect_true(all(is.na(fit$fit_stat$p.X2)))
  expect_false(anyNA(fit$fit_stat$p.G2))
})

test_that("plot.irtfit() keeps the colors of the standardized residuals when all or none exceed overSR", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  layer_look <- function(over) {
    fit <- irtfit(x_ft, score = theta_fit, data = resp_ft, D = 1, overSR = over)
    suppressWarnings(plot(fit, item.loc = 1, type = "sr", show.table = FALSE))
    d <- suppressWarnings(ggplot2::ggplot_build(last_ggplot()))$data[[1]]
    list(colour = unique(d$colour), shape = unique(d$shape))
  }

  # every residual exceeds the threshold: red circles
  all_over <- layer_look(0)
  expect_identical(all_over$colour, "red")
  expect_identical(as.numeric(all_over$shape), 1)

  # no residual exceeds the threshold: blue crosses
  none_over <- layer_look(100)
  expect_identical(none_over$colour, "blue")
  expect_identical(as.numeric(none_over$shape), 4)
})

test_that("plot.irtfit() with type = 'both' keeps the colors of the standardized residuals", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  # type = 'both' draws the two plots with grid.arrange(), so the plots are
  # taken from the arguments of that call
  layer_look_both <- function(over) {
    store <- new.env()
    ns <- asNamespace("gridExtra")
    suppressMessages(base::trace(
      "grid.arrange", where = ns, print = FALSE,
      tracer = bquote(assign("plots", list(...), envir = .(store)))
    ))
    on.exit(suppressMessages(base::untrace("grid.arrange", where = ns)), add = TRUE)
    fit <- irtfit(x_ft, score = theta_fit, data = resp_ft, D = 1, overSR = over)
    suppressWarnings(plot(fit, item.loc = 1, type = "both", show.table = FALSE))
    d <- suppressWarnings(ggplot2::ggplot_build(store$plots[[2]]))$data[[1]]
    list(colour = unique(d$colour), shape = unique(d$shape))
  }

  all_over <- layer_look_both(0)
  expect_identical(all_over$colour, "red")
  expect_identical(as.numeric(all_over$shape), 1)

  none_over <- layer_look_both(100)
  expect_identical(none_over$colour, "blue")
  expect_identical(as.numeric(none_over$shape), 4)
})
