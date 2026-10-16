# confirm_df() is not exported; access via :::
confirm_df <- irtQ:::confirm_df

# ---- helpers -------------------------------------------------------------------

make_drm_df <- function(model = "3PLM") {
  data.frame(
    id    = c("I1", "I2"),
    cats  = c(2L, 2L),
    model = model,
    par.1 = c(1.0, 1.2),
    par.2 = c(-0.5, 0.5),
    par.3 = c(0.2, 0.15),
    stringsAsFactors = FALSE
  )
}

make_grm_df <- function() {
  data.frame(
    id    = c("I1", "I2"),
    cats  = c(4L, 4L),
    model = "GRM",
    par.1 = c(1.0, 1.2),
    par.2 = c(-1.0, -0.5),
    par.3 = c(0.0, 0.5),
    par.4 = c(1.0, 1.5),
    stringsAsFactors = FALSE
  )
}

# ---- valid inputs --------------------------------------------------------------

test_that("confirm_df() accepts a valid 3PLM data frame", {
  x <- make_drm_df("3PLM")
  result <- confirm_df(x)
  expect_s3_class(result, "data.frame")
  expect_equal(colnames(result)[1:3], c("id", "cats", "model"))
})

test_that("confirm_df() accepts a valid 2PLM data frame and sets par.3 = 0", {
  x <- data.frame(
    id = "I1", cats = 2L, model = "2PLM",
    par.1 = 1.0, par.2 = 0.0, par.3 = NA_real_,
    stringsAsFactors = FALSE
  )
  result <- confirm_df(x)
  expect_equal(result$par.3, 0)
})

test_that("confirm_df() accepts a valid 1PLM data frame and sets par.3 = 0", {
  x <- data.frame(
    id = "I1", cats = 2L, model = "1PLM",
    par.1 = 1.0, par.2 = 0.5, par.3 = NA_real_,
    stringsAsFactors = FALSE
  )
  result <- confirm_df(x)
  expect_equal(result$par.3, 0)
})

test_that("confirm_df() accepts a valid GRM data frame", {
  x <- make_grm_df()
  result <- confirm_df(x)
  expect_equal(result$model, c("GRM", "GRM"))
})

test_that("confirm_df() accepts a valid GPCM data frame", {
  x <- data.frame(
    id = "I1", cats = 4L, model = "GPCM",
    par.1 = 1.2, par.2 = -0.5, par.3 = 0.0, par.4 = 0.8,
    stringsAsFactors = FALSE
  )
  expect_no_error(confirm_df(x))
})

test_that("confirm_df() converts DRM to 3PLM with a warning", {
  x <- make_drm_df("DRM")
  expect_warning(result <- confirm_df(x), "treated as '3PLM'")
  expect_true(all(result$model == "3PLM"))
})

test_that("confirm_df() normalises model names to upper case", {
  x <- make_drm_df("3plm")
  result <- confirm_df(x)
  expect_equal(result$model, c("3PLM", "3PLM"))
})

test_that("confirm_df() renames columns to id/cats/model/par.1/par.2/...", {
  x <- make_drm_df("2PLM")
  x$par.3 <- NA_real_
  colnames(x) <- c("item", "ncat", "irt", "a", "b", "g")
  result <- confirm_df(x)
  expect_equal(colnames(result)[1:3], c("id", "cats", "model"))
  expect_true("par.1" %in% colnames(result))
})

test_that("confirm_df() g2na=TRUE sets par.3 to NA for 1PLM/2PLM items", {
  x <- make_drm_df("2PLM")
  x$par.3 <- NA_real_
  result <- confirm_df(x, g2na = TRUE)
  expect_true(all(is.na(result$par.3)))
})

test_that("confirm_df() drops all-NA trailing parameter columns", {
  x <- make_drm_df("3PLM")
  x$par.4 <- NA_real_
  x$par.5 <- NA_real_
  result <- confirm_df(x)
  expect_false("par.5" %in% colnames(result))
})

test_that("confirm_df() accepts a factor model column and converts it", {
  x <- make_drm_df("3PLM")
  x$model <- factor(x$model)
  result <- confirm_df(x)
  expect_true(is.character(result$model))
})

# ---- invalid inputs ------------------------------------------------------------

test_that("confirm_df() errors on unknown model name", {
  x <- make_drm_df("4PLM")
  expect_error(confirm_df(x), "mis-specified")
})

test_that("confirm_df() errors when cats < 1 (cats = 0)", {
  x <- make_drm_df("3PLM")
  x$cats[1] <- 0L
  expect_error(confirm_df(x), "score category")
})

test_that("confirm_df() errors when par.3 column missing and model is not all-2PLM", {
  x <- data.frame(
    id = "I1", cats = 2L, model = "3PLM",
    par.1 = 1.0, par.2 = 0.0,
    stringsAsFactors = FALSE
  )
  expect_error(confirm_df(x), "par.3")
})

test_that("confirm_df() errors when a parameter column holds text", {
  x <- make_drm_df("2PLM")
  x$par.1 <- c("1.2", "0.9")
  expect_error(confirm_df(x), "numeric")

  # an all-NA logical column is not text
  y <- make_drm_df("2PLM")
  y$par.3 <- NA
  expect_no_error(confirm_df(y))
})

test_that("confirm_df() errors when cats is missing or less than 2", {
  x <- make_drm_df("3PLM")
  x$cats[2] <- 1L
  expect_error(confirm_df(x), "score category")

  x$cats[2] <- NA
  expect_error(confirm_df(x), "score category")
})

test_that("confirm_df() errors when a dichotomous model has cats other than 2", {
  x <- make_drm_df("3PLM")
  x$cats[1] <- 3L
  expect_error(confirm_df(x), "cats = 2")
})

test_that("confirm_df() errors when a polytomous item has more thresholds than cats - 1", {
  # two GRM items with 3 and 4 categories, where the first has a third threshold
  x <- data.frame(
    id = c("A", "B"), cats = c(3L, 4L), model = "GRM", par.1 = c(1, 1.2),
    par.2 = c(-1, -1), par.3 = c(0, 0), par.4 = c(1, 1)
  )
  expect_error(confirm_df(x), "more threshold parameters")

  # NA in the unused threshold cell is accepted
  x$par.4[1] <- NA
  expect_no_error(confirm_df(x))
})

test_that("confirm_df() accepts a GRM item with two categories", {
  x <- data.frame(
    id = c("A", "B"), cats = c(2L, 3L), model = "GRM", par.1 = c(1, 1.2),
    par.2 = c(0.3, -1), par.3 = c(NA, 0.5)
  )
  r <- confirm_df(x)
  expect_equal(r$cats, c(2, 3))
})

test_that("confirm_df() adds par.3 when every item is 1PLM or 2PLM", {
  # all 1PLM items without a par.3 column
  x1 <- data.frame(id = c("A", "B"), cats = 2L, model = "1PLM", par.1 = 1, par.2 = c(0, 1))
  r1 <- confirm_df(x1)
  expect_equal(r1$par.3, c(0, 0))

  # a mixture of 1PLM and 2PLM items without a par.3 column
  x2 <- data.frame(
    id = c("A", "B"), cats = 2L, model = c("1PLM", "2PLM"),
    par.1 = c(1, 1.3), par.2 = c(0, 1)
  )
  r2 <- confirm_df(x2)
  expect_equal(r2$par.3, c(0, 0))
  expect_equal(r2$par.2, c(0, 1))
})

test_that("confirm_df() errors when a two-category GRM or GPCM item has a value in par.3", {
  for (m in c("GRM", "GPCM")) {
    x <- data.frame(
      id = c("A", "B"), cats = c(2L, 3L), model = m, par.1 = c(1, 1.2),
      par.2 = c(0.3, -1), par.3 = c(0.2, 0.5)
    )
    expect_error(confirm_df(x), "more threshold parameters")

    # the cell must be NA
    x$par.3[1] <- NA
    expect_no_error(confirm_df(x))
  }
})

test_that("confirm_df() reads a cats column stored as text or factor labels", {
  x <- make_grm_df()
  x$cats <- as.character(x$cats)
  r1 <- confirm_df(x)
  expect_true(is.numeric(r1$cats))
  expect_equal(r1$cats, c(4, 4))

  # factor labels are read as the numbers they show
  x$cats <- factor(c(4, 4))
  expect_equal(confirm_df(x)$cats, c(4, 4))

  # the result equals that of numeric cats
  expect_equal(confirm_df(x)$par.2, confirm_df(make_grm_df())$par.2)

  # text that is not a number stops with a message
  x$cats <- c("4", "many")
  expect_error(confirm_df(x), "cats column")
})

test_that("info() and simdat() accept item metadata with a text cats column", {
  meta <- data.frame(
    id = c("A", "B"), cats = c("2", "4"), model = c("2PLM", "GRM"), par.1 = c(1, 1.2),
    par.2 = c(0, -1), par.3 = c(NA, 0), par.4 = c(NA, 1)
  )
  num <- meta
  num$cats <- c(2, 4)
  expect_equal(info(meta, c(-1, 1))$iif, info(num, c(-1, 1))$iif)
  expect_equal(dim(simdat(meta, theta = c(-1, 0, 1))), c(3L, 2L))
})
