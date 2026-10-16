# shape_df() builds the item metadata, including the default parameter values

test_that("shape_df() with default.par repeats a single cats value for every item", {
  models <- c("1PLM", "2PLM", "2PLM", "3PLM", "3PLM")
  meta1 <- shape_df(cats = 2, model = models, default.par = TRUE)
  # the same input written out in full
  meta2 <- shape_df(cats = rep(2, 5), model = models, default.par = TRUE)
  expect_identical(meta1, meta2)
  # item ids follow the number of items
  expect_identical(meta1$id, paste0("V", 1:5))
  # guessing parameters are NA for 1PLM and 2PLM items and 0.2 for 3PLM items
  expect_identical(meta1$par.3, c(NA, NA, NA, 0.2, 0.2))
})

test_that("shape_df() with default.par repeats a single cats value for polytomous items", {
  meta1 <- shape_df(cats = 3, model = c("GRM", "GPCM"), default.par = TRUE)
  meta2 <- shape_df(cats = c(3, 3), model = c("GRM", "GPCM"), default.par = TRUE)
  expect_identical(meta1, meta2)
  expect_identical(meta1$id, c("V1", "V2"))
  expect_identical(meta1$cats, c(3, 3))
})

test_that("shape_df() with default.par takes the number of items from item.id", {
  meta1 <- shape_df(cats = 2, model = "3PLM", item.id = paste0("it", 1:3), default.par = TRUE)
  expect_identical(meta1$id, paste0("it", 1:3))
  expect_identical(meta1$par.3, rep(0.2, 3))
})

test_that("shape_df() with default.par repeats a single model for every item", {
  meta1 <- shape_df(cats = rep(2, 4), model = "3PLM", default.par = TRUE)
  meta2 <- shape_df(cats = rep(2, 4), model = rep("3PLM", 4), default.par = TRUE)
  expect_identical(meta1, meta2)
  expect_identical(meta1$par.3, rep(0.2, 4))
})

test_that("shape_df() with default.par gives g = 0.2 to every 3PLM item", {
  meta <- shape_df(cats = 2, model = c("3PLM", "2PLM", "3PLM", "1PLM"), default.par = TRUE)
  expect_equal(meta$par.3, c(0.2, NA, 0.2, NA))
})

test_that("shape_df() errors when cats, model, or item.id do not match the number of items", {
  # model length 2 does not match 4 items
  expect_error(
    shape_df(cats = rep(2, 4), model = c("3PLM", "2PLM"), default.par = TRUE),
    "length 1 or one value per item"
  )
  expect_error(
    shape_df(cats = rep(2, 4), model = "3PLM", item.id = c("a", "b"), default.par = TRUE),
    "item.id"
  )

  # model length 2 does not match the 3 items given by the parameters
  expect_error(
    shape_df(
      par.drm = list(a = c(1, 1.2, 0.9), b = c(0, 1, -1), g = rep(0.2, 3)),
      cats = rep(2, 3), model = c("3PLM", "2PLM")
    ),
    "length 1 or one value per item"
  )
})

test_that("shape_df() errors when the parameter vectors do not match the items", {
  # slope vector shorter than the polytomous items
  expect_error(
    shape_df(
      par.prm = list(a = 1, d = list(c(0, 1), c(-1, 0, 1))),
      cats = c(3, 4), model = "GRM"
    ),
    "par.prm"
  )

  # threshold vector length does not equal cats - 1
  expect_error(
    shape_df(par.prm = list(a = 1, d = list(c(-1, 0, 1))), cats = 3, model = "GRM"),
    "par.prm"
  )

  # difficulty vector shorter than the slope vector
  expect_error(
    shape_df(par.drm = list(a = c(1, 1.2), b = 0), cats = c(2, 2), model = "2PLM"),
    "par.drm"
  )
})

test_that("shape_df() accepts par.drm without the g element", {
  meta <- shape_df(par.drm = list(a = c(1, 1.2), b = c(0, 1)), cats = 2, model = "2PLM")
  expect_equal(meta$par.1, c(1, 1.2))
  expect_equal(meta$par.2, c(0, 1))
  expect_true(all(is.na(meta$par.3)))

  # a NULL g element gives zero guessing values for 3PLM items
  meta3 <- shape_df(par.drm = list(a = 1, b = 0, g = NULL), cats = 2, model = "3PLM")
  expect_equal(meta3$par.3, 0)
})

test_that("startval_df() repeats a single cats value like a full vector", {
  models <- c("1PLM", "3PLM", "DRM")
  expect_identical(
    irtQ:::startval_df(cats = 2, model = models),
    irtQ:::startval_df(cats = rep(2, 3), model = models)
  )
})

test_that("shape_df() uses a single g value for every dichotomous item", {
  meta <- shape_df(
    par.drm = list(a = c(0.8, 1.5, 0.8, 1.5), b = c(0, 0, -0.5, -0.5), g = 0.15),
    cats = 2, model = "3PLM"
  )
  expect_equal(meta$par.3, rep(0.15, 4))

  # a single value also works next to polytomous items
  mixed <- shape_df(
    par.drm = list(a = c(1, 1.2), b = c(0, 0.5), g = 0.2),
    par.prm = list(a = 1, d = list(c(-1, 1))), cats = c(2, 3, 2), model = c("3PLM", "GRM", "3PLM")
  )
  expect_equal(mixed$par.3[c(1, 3)], c(0.2, 0.2))

  # any other wrong length still stops
  expect_error(
    shape_df(
      par.drm = list(a = c(0.8, 1.5, 0.8), b = c(0, 0, -0.5), g = c(0.1, 0.2)),
      cats = 2, model = "3PLM"
    ),
    "par.drm"
  )
})
