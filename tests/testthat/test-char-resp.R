# est_irt() and est_mg() read character response columns by their values

drop_time <- function(fit) {
  fit[c("EMtime", "SEtime", "TotalTime", "call", "data")] <- NULL
  fit
}

test_that("est_irt() fits character response columns as the numeric ones", {
  skip_on_cran()
  set.seed(7)
  x <- shape_df(par.drm = list(a = c(1, 1.3, 0.8, 1.1), b = c(-1, 0, 0.5, 1), g = 0), cats = 2, model = "2PLM")
  dat <- as.data.frame(simdat(x = x, theta = rnorm(400), D = 1))
  chr <- dat
  chr[, 1:2] <- lapply(chr[, 1:2], as.character)
  f_num <- est_irt(data = dat, D = 1, model = "2PLM", cats = 2, verbose = FALSE)
  f_chr <- est_irt(data = chr, D = 1, model = "2PLM", cats = 2, verbose = FALSE)
  expect_identical(drop_time(f_chr), drop_time(f_num))
  # the fixed item parameter calibration as well
  g_num <- est_irt(x = x, data = dat, D = 1, fipc = TRUE, fix.loc = 1:2, verbose = FALSE)
  g_chr <- est_irt(x = x, data = chr, D = 1, fipc = TRUE, fix.loc = 1:2, verbose = FALSE)
  expect_identical(drop_time(g_chr), drop_time(g_num))
})

test_that("est_mg() fits a character response column in only one group", {
  set.seed(8)
  x <- shape_df(par.drm = list(a = c(1, 1.3, 0.8, 1.1), b = c(-1, 0, 0.5, 1), g = 0), cats = 2, model = "2PLM")
  d1 <- as.data.frame(simdat(x = x, theta = rnorm(300), D = 1))
  d2 <- as.data.frame(simdat(x = x, theta = rnorm(300, 0.3), D = 1))
  d2c <- d2
  d2c[, 1] <- as.character(d2c[, 1])
  f_num <- est_mg(
    x = list(x, x), data = list(d1, d2), group.name = c("A", "B"), D = 1, free.group = 2, verbose = FALSE
  )
  f_chr <- est_mg(
    x = list(x, x), data = list(d1, d2c), group.name = c("A", "B"), D = 1, free.group = 2, verbose = FALSE
  )
  expect_identical(drop_time(f_chr), drop_time(f_num))
})

test_that("est_irt() fits a character column with a non-numeric missing code as the numeric one", {
  skip_on_cran()
  set.seed(9)
  x <- shape_df(par.drm = list(a = c(1, 1.3, 0.8, 1.1), b = c(-1, 0, 0.5, 1), g = 0), cats = 2, model = "2PLM")
  dat <- as.data.frame(simdat(x = x, theta = rnorm(400), D = 1))
  dat[c(5, 20, 77), 1] <- NA
  chr <- dat
  chr[, 1] <- as.character(chr[, 1])
  chr[is.na(dat[, 1]), 1] <- "."
  f_num <- est_irt(data = dat, D = 1, model = "2PLM", cats = 2, verbose = FALSE)
  f_chr <- est_irt(data = chr, D = 1, model = "2PLM", cats = 2, missing = ".", verbose = FALSE)
  expect_identical(drop_time(f_chr), drop_time(f_num))
  # the fixed item parameter calibration as well
  g_num <- est_irt(x = x, data = dat, D = 1, fipc = TRUE, fix.loc = 1:2, verbose = FALSE)
  g_chr <- est_irt(x = x, data = chr, D = 1, fipc = TRUE, fix.loc = 1:2, missing = ".", verbose = FALSE)
  expect_identical(drop_time(g_chr), drop_time(g_num))
})

test_that("est_irt() and est_mg() recode a missing code in a data frame that also has NA", {
  skip_on_cran()
  set.seed(13)
  x <- shape_df(
    par.drm = list(a = c(1, 1.3, 0.8, 1.1), b = c(-1, 0, 0.5, 1), g = 0),
    cats = 2, model = "2PLM"
  )
  mat <- simdat(x = x, theta = rnorm(400), D = 1)
  mat[c(5, 20), 1] <- NA
  mat[c(7, 60), 2] <- 99
  dat <- as.data.frame(mat)
  f_mat <- est_irt(data = mat, D = 1, model = "2PLM", cats = 2, missing = 99,
                   verbose = FALSE)
  f_dat <- est_irt(data = dat, D = 1, model = "2PLM", cats = 2, missing = 99,
                   verbose = FALSE)
  expect_identical(drop_time(f_dat), drop_time(f_mat))
  # the same data split into two groups
  grp <- rep(1:2, length.out = nrow(mat))
  fit_mg <- function(d) {
    est_mg(
      x = list(x, x), data = list(d[grp == 1, ], d[grp == 2, ]), group.name = c("A", "B"),
      D = 1, fipc = TRUE, fix.loc = list(1:2, 1:2), free.group = 2, missing = 99,
      verbose = FALSE
    )
  }
  expect_identical(drop_time(fit_mg(dat)), drop_time(fit_mg(mat)))
})
