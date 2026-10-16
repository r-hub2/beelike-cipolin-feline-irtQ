# ---- helpers -------------------------------------------------------------------

# item metadata with a 3PLM item, a GRM item, and a GPCM item
make_mixed_meta <- function() {
  data.frame(
    id = c("D3", "G4", "P5"), cats = c(2, 4, 5), model = c("3PLM", "GRM", "GPCM"),
    par.1 = c(1.2, 1.4, 0.9), par.2 = c(0.4, -1.2, -1), par.3 = c(0.18, 0.1, 0),
    par.4 = c(NA, 1.3, 0.7), par.5 = c(NA, NA, 1.6)
  )
}

# ---- estimation objects --------------------------------------------------------

test_that("traceline() for est_irt and est_item objects uses the stored metadata and D", {
  meta <- make_mixed_meta()
  theta <- seq(-3, 3, 1)
  ref <- traceline(meta, theta, D = 1.702)
  fit_irt <- structure(list(par.est = meta, scale.D = 1.702), class = c("est_irt", "list"))
  fit_item <- structure(list(par.est = meta, scale.D = 1.702), class = c("est_item", "list"))
  expect_identical(traceline(fit_irt, theta), ref)
  expect_identical(traceline(fit_item, theta), ref)

  # a D supplied through ... does not replace the stored scaling constant
  expect_identical(traceline(fit_irt, theta, D = 1), ref)
})

# ---- structure of the results ---------------------------------------------------

test_that("traceline() keeps one-row matrices for a single theta", {
  meta <- make_mixed_meta()
  tr <- traceline(meta, 0)
  expect_true(all(vapply(tr$prob.cats, is.matrix, logical(1))))
  expect_equal(vapply(tr$prob.cats, ncol, numeric(1)), c(D3 = 2, G4 = 4, P5 = 5))
  expect_equal(vapply(tr$prob.cats, nrow, numeric(1)), c(D3 = 1, G4 = 1, P5 = 1))

  # the test characteristic curve can be plotted for a single theta
  grDevices::pdf(NULL)
  p <- plot(tr)
  grDevices::dev.off()
  expect_s3_class(p, "ggplot")
})

test_that("traceline() ICC and TCC are expected scores of the category probabilities", {
  meta <- make_mixed_meta()
  tr <- traceline(meta, seq(-3, 3, 1))
  es <- sapply(tr$prob.cats, function(p) p %*% (0:(ncol(p) - 1)))
  expect_equal(unname(tr$icc), unname(es), tolerance = 1e-8)
  expect_equal(tr$tcc, rowSums(tr$icc))
  expect_true(all(abs(sapply(tr$prob.cats, rowSums) - 1) < 1e-8))
})

test_that("traceline() ICC ignores the categories an item does not have", {
  # a GRM item with 3 categories next to a GRM item with 5 categories
  meta <- data.frame(
    id = c("A", "B"), cats = c(3, 5), model = "GRM", par.1 = c(1, 1.2), par.2 = c(-1, -1),
    par.3 = c(1, 0), par.4 = c(NA, 1), par.5 = c(NA, 2)
  )
  tr <- traceline(meta, c(-2, 0, 1.5))
  es <- sapply(tr$prob.cats, function(p) p %*% (0:(ncol(p) - 1)))
  expect_equal(unname(tr$icc), unname(es), tolerance = 1e-12)
  expect_equal(tr$tcc, rowSums(es), tolerance = 1e-12)

  # the same holds for GPCM items
  meta$model <- "GPCM"
  tr <- traceline(meta, c(-2, 0, 1.5))
  es <- sapply(tr$prob.cats, function(p) p %*% (0:(ncol(p) - 1)))
  expect_equal(unname(tr$icc), unname(es), tolerance = 1e-12)
})

test_that("traceline() computes a GRM or GPCM item with two categories as a 2PLM item", {
  theta <- c(-2, -1, 0, 0.5, 1.5)
  two_cat <- data.frame(
    id = c("G", "P"), cats = 2, model = c("GRM", "GPCM"), par.1 = c(1.2, 0.8),
    par.2 = c(0.3, -0.4), par.3 = NA
  )
  tr <- traceline(two_cat, theta, D = 1.702)

  # the probabilities equal those of prm() with a single threshold
  expect_equal(tr$prob.cats$G[, 2], prm(theta, 1.2, 0.3, 1.702, "GRM")[, 2], tolerance = 1e-12)
  expect_equal(tr$prob.cats$P[, 2], prm(theta, 0.8, -0.4, 1.702, "GPCM")[, 2], tolerance = 1e-12)

  # and the 2PLM items with the same parameters
  twopl <- data.frame(
    id = c("G", "P"), cats = 2, model = "2PLM", par.1 = c(1.2, 0.8),
    par.2 = c(0.3, -0.4), par.3 = NA
  )
  expect_equal(tr$icc, traceline(twopl, theta, D = 1.702)$icc, tolerance = 1e-12)
})

test_that("traceline() treats NA guessing parameters of 3PLM and DRM items as zeros", {
  na_g <- data.frame(
    id = c("A", "B"), cats = 2, model = c("3PLM", "DRM"), par.1 = 1, par.2 = 0, par.3 = NA
  )
  tr <- suppressWarnings(traceline(na_g, c(-1, 0, 1)))
  expect_equal(unname(tr$icc[, 1]), plogis(c(-1, 0, 1)), tolerance = 1e-12)
  expect_equal(unname(tr$icc[, 2]), plogis(c(-1, 0, 1)), tolerance = 1e-12)
})

# ---- plot labels ---------------------------------------------------------------

test_that("plot.traceline() labels items and score panels as given", {
  meta <- data.frame(
    id = c("1", "item-2", "CR 3"), cats = c(2, 2, 3), model = c("3PLM", "2PLM", "GRM"),
    par.1 = c(1, 1.2, 1), par.2 = c(0, 0.5, -1), par.3 = c(0.2, NA, 1)
  )
  tr <- traceline(meta, seq(-2, 2, 0.5))
  grDevices::pdf(NULL)
  p_curve <- plot(tr, item.loc = 1:3, score.curve = TRUE, overlap = TRUE)
  p_panel <- plot(tr, item.loc = 3)
  p_over <- plot(tr, item.loc = 3, overlap = TRUE)
  grDevices::dev.off()

  # item IDs are not converted to syntactic names
  expect_setequal(as.character(unique(p_curve$data$Item)), c("1", "item-2", "CR 3"))

  # the panel titles show the score value only
  expect_setequal(unique(p_panel$data$score), c("Score: 0", "Score: 1", "Score: 2"))
  expect_setequal(as.character(unique(p_over$data$score)), c("0", "1", "2"))
})

# ---- a flexMIRT test form ------------------------------------------------------

test_that("traceline() of the flexMIRT sample form gives expected scores and keeps matrices", {
  meta <- bring.flexmirt(
    file = system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ"), "par"
  )$Group1$full_df

  # ICC and TCC are expected scores, and the category probabilities sum to 1
  tr <- traceline(meta, seq(-3, 3, 1))
  es <- sapply(tr$prob.cats, function(p) p %*% (0:(ncol(p) - 1)))
  expect_lt(max(abs(tr$icc - es)), 1e-12)
  expect_equal(tr$tcc, rowSums(tr$icc))
  expect_true(all(abs(sapply(tr$prob.cats, rowSums) - 1) < 1e-8))

  # a single theta keeps one-row matrices for every item
  tr1 <- traceline(meta, 0)
  expect_true(all(vapply(tr1$prob.cats, is.matrix, logical(1))))
  expect_equal(dim(tr1$prob.cats[[39]]), c(1L, 5L))
})

test_that("traceline() computes a two-category GRM or GPCM item next to polytomous items of the same model", {
  theta <- c(-1, 0, 1)
  for (m in c("GRM", "GPCM")) {
    meta <- data.frame(
      id = c("g2", "g3", "g4"), cats = c(2, 3, 4), model = m, par.1 = c(1.2, 1.1, 0.8),
      par.2 = c(-0.5, -1, -1), par.3 = c(NA, 1, 0), par.4 = c(NA, NA, 1)
    )
    tr <- traceline(meta, theta, D = 1.702)
    expect_equal(vapply(tr$prob.cats, ncol, numeric(1)), c(g2 = 2, g3 = 3, g4 = 4))

    # the two-category item has the probabilities of a 2PLM item
    twopl <- data.frame(id = "g2", cats = 2, model = "2PLM", par.1 = 1.2, par.2 = -0.5, par.3 = NA)
    expect_equal(tr$prob.cats$g2, traceline(twopl, theta, D = 1.702)$prob.cats$g2, tolerance = 1e-12)

    # the polytomous items do not depend on the two-category item
    expect_equal(tr$prob.cats[2:3], traceline(meta[2:3, ], theta, D = 1.702)$prob.cats, tolerance = 1e-12)
  }
})
