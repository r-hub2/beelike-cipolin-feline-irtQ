# write.flexmirt() writes item lines that bring.flexmirt() and flexMIRT can read

# step parameters of GPCM items with 3, 4, 5, and 6 score categories
gpcm_items <- function() {
  shape_df(
    par.prm = list(
      a = c(1.1, 0.9, 1.3, 0.7),
      d = list(c(-0.5, 0.6), c(-1, 0.2, 0.9), c(-1.2, -0.3, 0.4, 1.1), c(-1.5, -0.6, 0, 0.5, 1.4))
    ),
    cats = c(3, 4, 5, 6), model = "GPCM", item.id = paste0("G", 3:6)
  )
}

test_that("write.flexmirt() writes GPCM items with any number of score categories", {
  gp <- gpcm_items()
  prm_file <- tempfile(fileext = "-prm.txt")
  write.flexmirt(gp, file = prm_file, rePar = FALSE)
  lines <- strsplit(readLines(prm_file), "\t")
  # one line for each item and one line for the group
  expect_length(lines, 5L)
  theta <- c(-2, 0, 1.5)
  for (i in 1:4) {
    K <- gp$cats[i]
    row <- suppressWarnings(as.numeric(lines[[i]]))
    # six leading fields, a zero, K - 1 scoring values, the slope, a zero, and K - 1 intercepts
    expect_length(row, 7 + 2 * K)
    # the scoring function matrix of the nominal model
    Tmat <- matrix(0, K, K - 1)
    Tmat[, 1] <- 0:(K - 1)
    for (k in 2:(K - 1)) for (j in 2:(K - 1)) Tmat[k, j] <- sin(pi * (j - 1) * (k - 1) / (K - 1))
    alpha <- row[8:(8 + K - 2)]
    slope <- row[8 + K - 1]
    gamma <- row[(8 + K + 1):(8 + 2 * K - 1)]
    # category probabilities of the nominal model computed from the values in the file
    logit <- outer(theta, slope * drop(Tmat %*% alpha)) + matrix(drop(Tmat %*% gamma), length(theta), K, byrow = TRUE)
    prob_file <- exp(logit) / rowSums(exp(logit))
    prob_irtq <- traceline(gp[i, ], theta = theta, D = 1)$prob.cats[[1]]
    expect_equal(unname(prob_file), unname(prob_irtq), tolerance = 1e-6)
  }
  # the step parameters are recovered by bring.flexmirt()
  back <- bring.flexmirt(prm_file, "par")$Group1$full_df
  expect_equal(as.matrix(back[, 4:9]), as.matrix(gp[, 4:9]), tolerance = 1e-6, ignore_attr = TRUE)
})

test_that("write.flexmirt() writes the files that bring.flexmirt() reads for GRM and 3PLM items", {
  grm <- shape_df(
    par.prm = list(a = c(1.1, 0.8), d = list(c(-1.2, -0.3, 0.4, 1.1), c(-1, 0, 0.6, 1.3))),
    cats = 5, model = "GRM"
  )
  drm <- shape_df(par.drm = list(a = c(1.2, 0.8), b = c(0.5, -0.4), g = c(0.2, 0.1)), cats = 2, model = "3PLM")
  for (meta in list(grm, drm)) {
    prm_file <- tempfile(fileext = "-prm.txt")
    write.flexmirt(meta, file = prm_file, rePar = FALSE)
    back <- bring.flexmirt(prm_file, "par")$Group1$full_df
    npar <- ncol(meta) - 3
    expect_equal(as.matrix(back[, 4:(3 + npar)]), as.matrix(meta[, 4:(3 + npar)]), tolerance = 1e-6, ignore_attr = TRUE)
  }
})

test_that("write.flexmirt() writes the same file whatever the console width", {
  gp <- gpcm_items()
  files <- character(2)
  for (k in 1:2) {
    testthat::local_reproducible_output(width = c(40, 10000)[k])
    files[k] <- tempfile(fileext = "-prm.txt")
    write.flexmirt(gp, file = files[k], rePar = FALSE)
  }
  expect_identical(readLines(files[1]), readLines(files[2]))
  # one line for each of the four items and one line for the group
  expect_length(readLines(files[1]), 5L)
})

test_that("write.flexmirt() converts the nominal parameters of GPCM items back with rePar = TRUE", {
  gp <- shape_df(par.prm = list(a = 0.9, d = list(c(-1, 0.2, 0.9))), cats = 4, model = "GPCM")
  first_file <- tempfile(fileext = "-prm.txt")
  write.flexmirt(gp, file = first_file, rePar = FALSE)
  # the nominal parameters of the first file, read without conversion
  nominal <- bring.flexmirt(first_file, "par", rePar.gpc = FALSE)$Group1$full_df
  second_file <- tempfile(fileext = "-prm.txt")
  write.flexmirt(nominal, file = second_file, rePar = TRUE)
  expect_identical(readLines(second_file), readLines(first_file))
})

test_that("write.flexmirt() writes the groups of a multiple-group data set", {
  prm_file <- tempfile(fileext = "-prm.txt")
  write.flexmirt(simMG$item.prm,
    file = prm_file, rePar = FALSE, mgroup = TRUE,
    norm.pop = list(c(0, 1), c(0.5, 0.8), c(-0.3, 1.3)), group.name = c("A", "B", "C")
  )
  lines <- strsplit(readLines(prm_file), "\t")
  group_lines <- lines[vapply(lines, function(z) z[1] == "0", logical(1))]
  expect_length(group_lines, 3L)
  # the second value of norm.pop is written to the variance field
  expect_equal(as.numeric(vapply(group_lines, `[`, character(1), 7)), c(1, 0.8, 1.3))
  back <- bring.flexmirt(prm_file, "par")
  expect_equal(unname(vapply(back, function(z) z$cov[1, 1], numeric(1))), c(1, 0.8, 1.3))
})

test_that("write.flexmirt() writes a two-category GPCM item as a 2PL item", {
  gp <- data.frame(id = "G2", cats = 2, model = "GPCM", par.1 = 1.2, par.2 = 0.4)
  prm_file <- tempfile(fileext = "-prm.txt")
  write.flexmirt(gp, file = prm_file, rePar = FALSE)
  row <- as.numeric(strsplit(readLines(prm_file)[1], "\t")[[1]][-2])
  # one scoring value of 1, the slope, and the intercept -a * b
  expect_equal(row[7:10], c(1, 1.2, 0, -1.2 * 0.4))
  back <- bring.flexmirt(prm_file, "par")$Group1$full_df
  expect_equal(unname(unlist(back[1, 4:5])), c(1.2, 0.4), tolerance = 1e-6)
})

test_that("write.flexmirt() closes the file when it stops with an error", {
  gp <- data.frame(id = "G2", cats = 3, model = "GPCM", par.1 = 1.2, par.2 = 0.4)
  prm_file <- tempfile(fileext = "-prm.txt")
  # getAllConnections() does not run gc(), which would close a left-open connection
  n_open <- length(getAllConnections())
  expect_error(write.flexmirt(gp, file = prm_file, rePar = FALSE))
  expect_equal(length(getAllConnections()), n_open)
})

test_that("write.flexmirt() stops when the file is not given", {
  gp <- shape_df(par.drm = list(a = 1, b = 0, g = 0), cats = 2, model = "2PLM")
  expect_error(write.flexmirt(gp, file = NULL), "destination of the file")
})
