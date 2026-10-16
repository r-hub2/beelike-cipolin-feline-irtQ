# ---- Fixtures ------------------------------------------------------------------

prm_file <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
x_full   <- bring.flexmirt(file = prm_file, "par")$Group1$full_df  # 55 items

## 1. 3PLM-only: 10 dichotomous items (from flexMIRT sample)
x_drm <- x_full[1:10, ]
set.seed(99)
theta_drm  <- rnorm(200)
resp_drm   <- simdat(x = x_drm, theta = theta_drm, D = 1)

## 2. GRM-only: 5 items, 4 categories (scored 0-3)
x_grm <- shape_df(
  par.prm = list(
    a = c(0.9, 1.1, 1.2, 1.0, 0.8),
    d = list(c(-1.2, -0.2, 0.8), c(-1.0, 0.0, 1.0),
             c(-0.8,  0.2, 1.2), c(-1.5,-0.3, 0.7),
             c(-0.9,  0.1, 0.9))
  ),
  cats = rep(4L, 5), model = "GRM"
)
set.seed(11)
theta_grm  <- rnorm(300)
resp_grm   <- simdat(x = x_grm, theta = theta_grm, D = 1)

## 3. GPCM-only: 5 items, 4 categories
x_gpcm <- shape_df(
  par.prm = list(
    a = c(1.0, 1.2, 0.9, 1.1, 1.0),
    d = list(c(-1.0, 0.0, 0.8), c(-0.8, 0.2, 1.0),
             c(-1.2,-0.2, 0.9), c(-0.5, 0.5, 1.2),
             c(-1.0,-0.1, 0.7))
  ),
  cats = rep(4L, 5), model = "GPCM"
)
set.seed(22)
theta_gpcm <- rnorm(300)
resp_gpcm  <- simdat(x = x_gpcm, theta = theta_gpcm, D = 1)

## 4. Mixed: 5 x 3PLM + 2 x GRM (4 cats)  -> max sum score = 5 + 6 = 11
x_mix_grm <- shape_df(
  par.drm = list(a = c(1.0,1.2,0.9,1.1,1.0),
                 b = c(-1.0,-0.5,0.0,0.5,1.0),
                 g = c(0.2, 0.15,0.1,0.2,0.15)),
  par.prm = list(a = c(1.0,1.2),
                 d = list(c(-1.0,0.0,1.0), c(-0.5,0.5,1.5))),
  cats  = c(2L,2L,2L,2L,2L,4L,4L),
  model = c(rep("3PLM",5), "GRM","GRM")
)
set.seed(33)
theta_mix_grm  <- rnorm(300)
resp_mix_grm   <- simdat(x = x_mix_grm, theta = theta_mix_grm, D = 1)

## 5. Mixed: 5 x 3PLM + 2 x GPCM (4 cats) -> max sum score = 5 + 6 = 11
x_mix_gpcm <- shape_df(
  par.drm = list(a = c(1.0,1.2,0.9,1.1,1.0),
                 b = c(-1.0,-0.5,0.0,0.5,1.0),
                 g = c(0.2, 0.15,0.1,0.2,0.15)),
  par.prm = list(a = c(1.0,1.2),
                 d = list(c(-1.0,0.0,1.0), c(-0.5,0.5,1.5))),
  cats  = c(2L,2L,2L,2L,2L,4L,4L),
  model = c(rep("3PLM",5), "GPCM","GPCM")
)
set.seed(44)
theta_mix_gpcm  <- rnorm(300)
resp_mix_gpcm   <- simdat(x = x_mix_gpcm, theta = theta_mix_gpcm, D = 1)


# ---- Helpers -------------------------------------------------------------------

check_pointwise <- function(result, n = 200) {
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), n)
  expect_true("est.theta" %in% colnames(result))
  expect_true("se.theta"  %in% colnames(result))
  expect_true(all(!is.nan(result$se.theta) | is.na(result$se.theta)))
}

check_sumtable <- function(result, max_ss) {
  expect_type(result, "list")
  expect_true("est.par"     %in% names(result))
  expect_true("score.table" %in% names(result))
  expect_true(all(c("sum.score","est.theta","se.theta") %in% colnames(result$est.par)))
  expect_equal(nrow(result$score.table), max_ss + 1L)
  expect_true(all(diff(result$score.table$est.theta) >= 0))
}


# ==============================================================================
# 1. 3PLM ONLY (10 dichotomous items, max sum = 10)
# ==============================================================================

base_drm <- list(x = x_drm, data = resp_drm, D = 1)

test_that("3PLM | ML: returns est.theta/se.theta data frame, correlates with truth", {
  res <- do.call(est_score, c(base_drm, list(method = "ML", range = c(-6, 6))))
  check_pointwise(res)
  expect_gt(cor(res$est.theta, theta_drm), 0.7)
})

test_that("3PLM | MLF: finite estimates within fence range", {
  res <- do.call(est_score, c(base_drm, list(method = "MLF", range = c(-5, 5))))
  check_pointwise(res)
  expect_true(all(is.finite(res$est.theta)))
  expect_true(all(res$est.theta >= -5 & res$est.theta <= 5))
})

test_that("3PLM | WL: returns est.theta/se.theta data frame", {
  res <- do.call(est_score, c(base_drm, list(method = "WL", range = c(-6, 6))))
  check_pointwise(res)
})

test_that("3PLM | MAP: returns est.theta/se.theta data frame", {
  res <- do.call(est_score, c(base_drm, list(method = "MAP",
                                              norm.prior = c(0,1), range = c(-6,6))))
  check_pointwise(res)
})

test_that("3PLM | EAP: estimates within quadrature range", {
  res <- do.call(est_score, c(base_drm, list(method = "EAP",
                                              norm.prior = c(0,1), nquad = 41L)))
  check_pointwise(res)
  expect_true(all(res$est.theta >= -4 & res$est.theta <= 4))
})

test_that("3PLM | EAP.SUM: list structure, score.table has 11 rows (0-10)", {
  res <- do.call(est_score, c(base_drm, list(method = "EAP.SUM",
                                              norm.prior = c(0,1), nquad = 41L)))
  check_sumtable(res, max_ss = 10L)
  expect_equal(nrow(res$est.par), 200L)
})

test_that("3PLM | INV.TCC: list structure, score.table monotone", {
  res <- do.call(est_score, c(base_drm, list(method = "INV.TCC", range.tcc = c(-7,7))))
  check_sumtable(res, max_ss = 10L)
})

test_that("3PLM | ML handles partial NA responses", {
  resp_na <- resp_drm; resp_na[1,1] <- NA; resp_na[2,2:3] <- NA
  res <- est_score(x=x_drm, data=resp_na, D=1, method="ML", missing=NA, range=c(-6,6))
  check_pointwise(res)
})

test_that("3PLM | EAP handles partial NA responses", {
  resp_na <- resp_drm; resp_na[5, 1:3] <- NA
  res <- est_score(x=x_drm, data=resp_na, D=1, method="EAP",
                   missing=NA, norm.prior=c(0,1))
  check_pointwise(res)
})

test_that("3PLM | ML/WL/MAP/EAP rank-correlate > 0.95", {
  ml  <- est_score(x=x_drm, data=resp_drm, D=1, method="ML",  range=c(-6,6))$est.theta
  wl  <- est_score(x=x_drm, data=resp_drm, D=1, method="WL",  range=c(-6,6))$est.theta
  map <- est_score(x=x_drm, data=resp_drm, D=1, method="MAP", range=c(-6,6),
                   norm.prior=c(0,1))$est.theta
  eap <- est_score(x=x_drm, data=resp_drm, D=1, method="EAP",
                   norm.prior=c(0,1))$est.theta
  expect_gt(cor(ml, wl,  method="spearman"), 0.95)
  expect_gt(cor(ml, map, method="spearman"), 0.95)
  expect_gt(cor(ml, eap, method="spearman"), 0.95)
})


# ==============================================================================
# 2. GRM ONLY (5 items, 4 cats -> max sum = 15)
# ==============================================================================

base_grm <- list(x = x_grm, data = resp_grm, D = 1)

test_that("GRM | ML: returns correct data frame structure (n=300)", {
  res <- do.call(est_score, c(base_grm, list(method = "ML", range = c(-6,6))))
  check_pointwise(res, n = 300L)
})

test_that("GRM | WL: returns correct data frame structure", {
  res <- do.call(est_score, c(base_grm, list(method = "WL", range = c(-6,6))))
  check_pointwise(res, n = 300L)
})

test_that("GRM | MAP: returns correct data frame structure", {
  res <- do.call(est_score, c(base_grm, list(method = "MAP",
                                              norm.prior = c(0,1), range = c(-6,6))))
  check_pointwise(res, n = 300L)
})

test_that("GRM | EAP: returns correct data frame structure", {
  res <- do.call(est_score, c(base_grm, list(method = "EAP",
                                              norm.prior = c(0,1), nquad = 41L)))
  check_pointwise(res, n = 300L)
})

test_that("GRM | EAP.SUM: list structure, score.table has 16 rows (0-15)", {
  res <- do.call(est_score, c(base_grm, list(method = "EAP.SUM",
                                              norm.prior = c(0,1), nquad = 41L)))
  check_sumtable(res, max_ss = 15L)
  expect_equal(nrow(res$est.par), 300L)
})

test_that("GRM | INV.TCC: list structure, score.table monotone (16 rows)", {
  res <- do.call(est_score, c(base_grm, list(method = "INV.TCC", range.tcc = c(-7,7))))
  check_sumtable(res, max_ss = 15L)
})

test_that("GRM | ML estimates correlate with true theta", {
  res <- do.call(est_score, c(base_grm, list(method = "ML", range = c(-6,6))))
  # 5 polytomous items: lower information than 10 dichotomous -> threshold 0.65
  expect_gt(cor(res$est.theta, theta_grm), 0.65)
})


# ==============================================================================
# 3. GPCM ONLY (5 items, 4 cats -> max sum = 15)
# ==============================================================================

base_gpcm <- list(x = x_gpcm, data = resp_gpcm, D = 1)

test_that("GPCM | ML: returns correct data frame structure (n=300)", {
  res <- do.call(est_score, c(base_gpcm, list(method = "ML", range = c(-6,6))))
  check_pointwise(res, n = 300L)
})

test_that("GPCM | WL: returns correct data frame structure", {
  res <- do.call(est_score, c(base_gpcm, list(method = "WL", range = c(-6,6))))
  check_pointwise(res, n = 300L)
})

test_that("GPCM | MAP: returns correct data frame structure", {
  res <- do.call(est_score, c(base_gpcm, list(method = "MAP",
                                               norm.prior = c(0,1), range = c(-6,6))))
  check_pointwise(res, n = 300L)
})

test_that("GPCM | EAP: returns correct data frame structure", {
  res <- do.call(est_score, c(base_gpcm, list(method = "EAP",
                                               norm.prior = c(0,1), nquad = 41L)))
  check_pointwise(res, n = 300L)
})

test_that("GPCM | EAP.SUM: list structure, score.table has 16 rows (0-15)", {
  res <- do.call(est_score, c(base_gpcm, list(method = "EAP.SUM",
                                               norm.prior = c(0,1), nquad = 41L)))
  check_sumtable(res, max_ss = 15L)
  expect_equal(nrow(res$est.par), 300L)
})

test_that("GPCM | INV.TCC: list structure, score.table monotone (16 rows)", {
  res <- do.call(est_score, c(base_gpcm, list(method = "INV.TCC", range.tcc = c(-7,7))))
  check_sumtable(res, max_ss = 15L)
})

test_that("GPCM | ML estimates correlate with true theta", {
  res <- do.call(est_score, c(base_gpcm, list(method = "ML", range = c(-6,6))))
  expect_gt(cor(res$est.theta, theta_gpcm), 0.7)
})


# ==============================================================================
# 4. MIXED: 5 x 3PLM + 2 x GRM (max sum = 5*1 + 2*3 = 11)
# ==============================================================================

base_mix_grm <- list(x = x_mix_grm, data = resp_mix_grm, D = 1)

test_that("Mixed(3PLM+GRM) | ML: returns correct data frame structure (n=300)", {
  res <- do.call(est_score, c(base_mix_grm, list(method = "ML", range = c(-6,6))))
  check_pointwise(res, n = 300L)
})

test_that("Mixed(3PLM+GRM) | WL: returns correct data frame structure", {
  res <- do.call(est_score, c(base_mix_grm, list(method = "WL", range = c(-6,6))))
  check_pointwise(res, n = 300L)
})

test_that("Mixed(3PLM+GRM) | MAP: returns correct data frame structure", {
  res <- do.call(est_score, c(base_mix_grm, list(method = "MAP",
                                                  norm.prior = c(0,1), range = c(-6,6))))
  check_pointwise(res, n = 300L)
})

test_that("Mixed(3PLM+GRM) | EAP: returns correct data frame structure", {
  res <- do.call(est_score, c(base_mix_grm, list(method = "EAP",
                                                  norm.prior = c(0,1), nquad = 41L)))
  check_pointwise(res, n = 300L)
})

test_that("Mixed(3PLM+GRM) | EAP.SUM: list structure, score.table has 12 rows (0-11)", {
  res <- do.call(est_score, c(base_mix_grm, list(method = "EAP.SUM",
                                                  norm.prior = c(0,1), nquad = 41L)))
  check_sumtable(res, max_ss = 11L)
  expect_equal(nrow(res$est.par), 300L)
})

test_that("Mixed(3PLM+GRM) | INV.TCC: list structure, score.table monotone", {
  res <- do.call(est_score, c(base_mix_grm, list(method = "INV.TCC", range.tcc = c(-7,7))))
  check_sumtable(res, max_ss = 11L)
})

test_that("Mixed(3PLM+GRM) | ML estimates correlate with true theta", {
  res <- do.call(est_score, c(base_mix_grm, list(method = "ML", range = c(-6,6))))
  expect_gt(cor(res$est.theta, theta_mix_grm), 0.7)
})


# ==============================================================================
# 5. MIXED: 5 x 3PLM + 2 x GPCM (max sum = 5*1 + 2*3 = 11)
# ==============================================================================

base_mix_gpcm <- list(x = x_mix_gpcm, data = resp_mix_gpcm, D = 1)

test_that("Mixed(3PLM+GPCM) | ML: returns correct data frame structure (n=300)", {
  res <- do.call(est_score, c(base_mix_gpcm, list(method = "ML", range = c(-6,6))))
  check_pointwise(res, n = 300L)
})

test_that("Mixed(3PLM+GPCM) | WL: returns correct data frame structure", {
  res <- do.call(est_score, c(base_mix_gpcm, list(method = "WL", range = c(-6,6))))
  check_pointwise(res, n = 300L)
})

test_that("Mixed(3PLM+GPCM) | MAP: returns correct data frame structure", {
  res <- do.call(est_score, c(base_mix_gpcm, list(method = "MAP",
                                                   norm.prior = c(0,1), range = c(-6,6))))
  check_pointwise(res, n = 300L)
})

test_that("Mixed(3PLM+GPCM) | EAP: returns correct data frame structure", {
  res <- do.call(est_score, c(base_mix_gpcm, list(method = "EAP",
                                                   norm.prior = c(0,1), nquad = 41L)))
  check_pointwise(res, n = 300L)
})

test_that("Mixed(3PLM+GPCM) | EAP.SUM: list structure, score.table has 12 rows (0-11)", {
  res <- do.call(est_score, c(base_mix_gpcm, list(method = "EAP.SUM",
                                                   norm.prior = c(0,1), nquad = 41L)))
  check_sumtable(res, max_ss = 11L)
  expect_equal(nrow(res$est.par), 300L)
})

test_that("Mixed(3PLM+GPCM) | INV.TCC: list structure, score.table monotone", {
  res <- do.call(est_score, c(base_mix_gpcm, list(method = "INV.TCC", range.tcc = c(-7,7))))
  check_sumtable(res, max_ss = 11L)
})

test_that("Mixed(3PLM+GPCM) | ML estimates correlate with true theta", {
  res <- do.call(est_score, c(base_mix_gpcm, list(method = "ML", range = c(-6,6))))
  expect_gt(cor(res$est.theta, theta_mix_gpcm), 0.7)
})


# ==============================================================================
# 6. INPUT CHECKS AND PARALLEL SCORING
# ==============================================================================

test_that("est_score() stops when data and item metadata have different numbers of items", {
  for (m in c("ML", "EAP", "EAP.SUM", "INV.TCC")) {
    expect_error(est_score(x_drm, resp_drm[1:3, 1:9], D = 1, method = m), "number of columns")
    expect_error(est_score(x_drm, cbind(resp_drm[1:3, ], 1), D = 1, method = m), "number of columns")
  }

  # a single response vector is checked in the same way
  expect_error(est_score(x_drm, resp_drm[1, 1:5], D = 1), "number of columns")
})

test_that("est_score() stops for responses that are not integer scores within the categories", {
  r <- resp_drm[1:3, ]

  # a score above the highest category of a dichotomous item
  r_high <- r
  r_high[1, 2] <- 2
  expect_error(est_score(x_drm, r_high, D = 1), "Responses must be integers")
  expect_error(est_score(x_drm, r_high, D = 1, method = "EAP"), "Responses must be integers")
  expect_error(est_score(x_drm, r_high, D = 1, method = "EAP.SUM"), "Responses must be integers")

  # a non-integer score
  r_frac <- r
  r_frac[2, 3] <- 0.7
  expect_error(est_score(x_drm, r_frac, D = 1), "Responses must be integers")

  # a negative code that is not declared as missing
  r_neg <- r
  r_neg[3, 1] <- -9
  expect_error(est_score(x_drm, r_neg, D = 1), "Responses must be integers")

  # the same code is accepted when it is declared as missing
  expect_identical(
    est_score(x_drm, r_neg, D = 1, missing = -9),
    est_score(x_drm, replace(r_neg, r_neg == -9, NA), D = 1)
  )

  # a polytomous item accepts its own categories only
  r_grm <- resp_grm[1:3, ]
  expect_no_error(est_score(x_grm, r_grm, D = 1))
  r_grm[1, 1] <- 4
  expect_error(est_score(x_grm, r_grm, D = 1), "Responses must be integers")
})

test_that("est_score() reads character and factor responses like numeric responses", {
  r <- resp_drm[1:6, ]
  r[2, 3] <- NA

  # a character matrix with a character code for missing responses
  r_chr <- matrix(as.character(r), nrow(r))
  r_chr[is.na(r)] <- "."

  # a data frame with character columns and a data frame with factor columns
  df_chr <- as.data.frame(r_chr, stringsAsFactors = FALSE)
  df_fac <- as.data.frame(lapply(df_chr, factor))

  for (m in c("ML", "MLF", "WL", "MAP", "EAP")) {
    ref <- est_score(x_drm, r, D = 1, method = m)
    expect_identical(est_score(x_drm, r_chr, D = 1, method = m, missing = "."), ref)
    expect_identical(est_score(x_drm, df_chr, D = 1, method = m, missing = "."), ref)
    expect_identical(est_score(x_drm, df_fac, D = 1, method = m, missing = "."), ref)
  }

  # the summed-score methods read the labels in the same way
  r0 <- resp_drm[1:6, ]
  r0_chr <- matrix(as.character(r0), nrow(r0))
  df0_fac <- as.data.frame(lapply(as.data.frame(r0_chr, stringsAsFactors = FALSE), factor))
  for (m in c("EAP.SUM", "INV.TCC")) {
    ref0 <- est_score(x_drm, r0, D = 1, method = m)
    expect_identical(est_score(x_drm, r0_chr, D = 1, method = m), ref0)
    expect_identical(est_score(x_drm, df0_fac, D = 1, method = m), ref0)
  }

  # a response that cannot be read as a number stops with a pointer to 'missing'
  expect_error(est_score(x_drm, r_chr, D = 1), "check the 'missing' argument")
  expect_error(est_score(x_drm, df_fac, D = 1, method = "EAP"), "check the 'missing' argument")

  # numeric labels outside the categories are still rejected
  r_chr_high <- r_chr
  r_chr_high[1, 1] <- "2"
  expect_error(est_score(x_drm, r_chr_high, D = 1, missing = "."), "Responses must be integers")
})

test_that("est_score() stops when the iteration limit or the tolerance is not valid", {
  r <- resp_drm[1:3, ]
  for (m in c("ML", "MLF", "WL", "MAP")) {
    for (bad in list(0, -1, 2.5, NA, Inf, c(10, 20), "10")) {
      expect_error(est_score(x_drm, r, D = 1, method = m, max.iter = bad), "'max.iter' must be")
    }
  }
  for (m in c("ML", "MLF", "WL", "MAP", "INV.TCC")) {
    for (bad in list(0, -1e-4, NA, Inf, c(1e-4, 1e-3), "1e-4")) {
      expect_error(est_score(x_drm, r, D = 1, method = m, tol = bad), "'tol' must be")
    }
  }

  # the limits of the valid range are accepted
  expect_no_error(est_score(x_drm, r, D = 1, method = "ML", max.iter = 1))
  expect_no_error(est_score(x_drm, r, D = 1, method = "ML", tol = 1e-12))
})

test_that("est_score() stops for an unknown method or starting value option", {
  r <- resp_drm[1:3, ]
  expect_error(est_score(x_drm, r, D = 1, method = "ml"), "'method' must be one of")
  expect_error(est_score(x_drm, r, D = 1, method = c("ML", "WL")), "'method' must be one of")
  expect_error(est_score(x_drm, r, D = 1, stval.opt = 4), "'stval.opt' must be 1, 2, or 3")
  expect_error(est_score(x_drm, r, D = 1, stval.opt = c(1, 2)), "'stval.opt' must be 1, 2, or 3")
})

test_that("est_score() for an est_irt object uses the stored data and D and rejects data or D", {
  fit <- structure(list(par.est = x_drm, data = resp_drm[1:5, ], scale.D = 1.702), class = "est_irt")
  expect_identical(
    est_score(fit, method = "ML"),
    est_score(x_drm, resp_drm[1:5, ], D = 1.702, method = "ML")
  )
  expect_error(est_score(fit, data = resp_drm[1:2, ], method = "ML"), "cannot be supplied")
  expect_error(est_score(fit, D = 1, method = "ML"), "cannot be supplied")
})

test_that("parallel scoring gives the serial result for any number of examinees per core", {
  skip_on_cran()
  n_conn <- nrow(showConnections())
  # use at most 2 cores when R CMD check --as-cran limits the cores of parallel
  limit_cores <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
  nc_set <- if (nzchar(limit_cores) && limit_cores != "false") 2 else c(2, 4)
  for (n in c(1, 2, 3, 5)) {
    r <- resp_drm[seq_len(n), , drop = FALSE]
    serial <- est_score(x_drm, r, D = 1, method = "ML")
    for (nc in nc_set) {
      expect_warning(par <- est_score(x_drm, r, D = 1, method = "ML", ncore = nc), "ncore > 1")
      expect_equal(par, serial)
    }
  }

  # the clusters are closed after scoring
  expect_equal(nrow(showConnections()), n_conn)
})

test_that("lwrc() accepts a single theta for a test with polytomous items", {
  x_pm <- x_full[c(1, 53), ]
  one <- lwrc(x_pm, theta = 0.4)
  two <- lwrc(x_pm, theta = c(0.4, 1))
  expect_equal(dim(one), c(6L, 1L))
  expect_equal(unname(one[, 1]), unname(two[, 1]))
  expect_equal(sum(one), 1)
})

test_that("llike_score() uses default fences for MLF and keeps examinees with duplicated row names apart", {
  x_ll <- x_full[c(1:6, 51:53), ]
  set.seed(8)
  d_ll <- simdat(x_ll, rnorm(3), D = 1)
  th <- seq(-2, 2, 1)

  # the default fences are c(-5, 5)
  expect_identical(
    llike_score(x_ll, d_ll, th, method = "MLF"),
    llike_score(x_ll, d_ll, th, method = "MLF", fence.b = c(-5, 5))
  )

  # two identical response rows with the same row name give two identical columns
  d_dup <- rbind(d_ll[1, ], d_ll[1, ], d_ll[2, ])
  rownames(d_dup) <- c("a", "a", "b")
  ll_dup <- llike_score(x_ll, d_dup, th)
  ref <- llike_score(x_ll, d_ll[1:2, ], th)
  expect_equal(unname(as.matrix(ll_dup)), unname(as.matrix(ref[, c(1, 1, 2)])))
})

test_that("llike_score() returns NA for an examinee without any observed response", {
  x_ll <- x_full[c(1:6, 51:53), ]
  set.seed(8)
  d_ll <- simdat(x_ll, rnorm(3), D = 1)
  th <- seq(-2, 2, 1)
  d_na <- rbind(d_ll[1, ], rep(NA, ncol(d_ll)), d_ll[2, ])
  for (m in c("ML", "MLF", "MAP")) {
    res <- llike_score(x_ll, d_na, th, method = m)
    expect_true(all(is.na(res[[2]])))
    expect_equal(res[, c(1, 3)], llike_score(x_ll, d_ll[1:2, ], th, method = m), ignore_attr = TRUE)

    # the other columns equal the values computed without that examinee
    expect_identical(
      unname(as.matrix(res[, c(1, 3)])),
      unname(as.matrix(llike_score(x_ll, d_ll[1:2, ], th, method = m)))
    )
  }

  # a data frame gives the same NA column for MLF
  res_df <- llike_score(x_ll, as.data.frame(d_na), th, method = "MLF")
  expect_true(all(is.na(res_df[[2]])))
  expect_false(anyNA(res_df[[1]]))
})

test_that("bisection() warns when the bounds do not bracket a root and stops at max.it", {
  f <- function(t, p) p - drm(theta = t, a = 1, b = 0.2, D = 1)

  # the root of f lies below lb, so both bounds have the same sign
  expect_warning(r <- bisection(f, p = 0.2, lb = 1, ub = 10), "opposite signs")
  expect_true(is.finite(r$root))

  # no warning when the bounds bracket the root
  expect_no_warning(bisection(f, p = 0.2, lb = -10, ub = 10))

  # exactly max.it iterations are done when the interval is still wide
  expect_warning(r <- bisection(f, p = 0.2, lb = -10, ub = 10, tol = 1e-12, max.it = 10), "maximum number")
  expect_equal(r$iter, 10)
  expect_equal(r$delta, 20 / 2^10)

  # reaching the tolerance in exactly max.it iterations gives no warning
  expect_no_warning(r <- bisection(f, p = 0.2, lb = -10, ub = 10, tol = 20 / 2^10, max.it = 10))
  expect_equal(r$iter, 10)
})

test_that("INV.TCC stops with a clear message when the TCC does not reach a sum score", {
  # the guessing parameters leave no sum score between their sum and the maximum
  x_g <- shape_df(
    par.drm = list(a = c(1, 1), b = c(0, 0.5), g = c(0.6, 0.6)),
    cats = 2, model = "3PLM"
  )
  expect_error(
    est_score(x_g, rbind(c(1, 0)), D = 1, method = "INV.TCC"),
    "does not reach the sum score"
  )

  # very flat items keep the TCC from reaching the lowest nonzero sum score
  x_flat <- shape_df(
    par.drm = list(a = rep(0.1, 10), b = rep(0, 10), g = rep(0, 10)),
    cats = 2, model = "3PLM"
  )
  expect_error(
    est_score(x_flat, rbind(rep(1, 10)), D = 1, method = "INV.TCC"),
    "does not reach the sum score"
  )

  # the message is not wrapped in an index error of the loop over the sum scores
  msg <- tryCatch(
    est_score(x_flat, rbind(rep(1, 10)), D = 1, method = "INV.TCC"),
    error = function(e) conditionMessage(e)
  )
  expect_false(grepl("index", msg))
})


# ==============================================================================
# 7. SCORING VALUES
# ==============================================================================

test_that("ML, WL, and MAP reach the stationary point when Fisher scoring cycles", {
  # Fisher scoring alternates between -2 and -3 for an all-incorrect pattern under WL
  wl <- est_score(x_drm, rbind(rep(0, 10)), D = 1, method = "WL")
  expect_equal(wl$est.theta, -2.591107, tolerance = 1e-3)

  # the same happens for ML with this pattern
  ml <- est_score(x_drm, rbind(c(1, 0, 0, 0, 1, 0, 0, 0, 0, 0)), D = 1, method = "ML")
  expect_equal(ml$est.theta, -2.761422, tolerance = 1e-3)

  # the estimate does not depend on the parity of max.iter
  for (m in c("WL", "ML", "MLF")) {
    a1 <- est_score(x_drm, rbind(rep(0, 10)), D = 1, method = m, max.iter = 100)
    a2 <- est_score(x_drm, rbind(rep(0, 10)), D = 1, method = m, max.iter = 101)
    expect_equal(a1$est.theta, a2$est.theta, tolerance = 1e-3)

    # a large iteration limit gives the same estimate
    a3 <- est_score(x_drm, rbind(rep(0, 10)), D = 1, method = m, max.iter = 3000)
    expect_equal(a1$est.theta, a3$est.theta, tolerance = 1e-3)
  }

  # MAP for a 15-item 3PLM test with only the first item answered correctly
  x15 <- shape_df(
    par.drm = list(
      a = c(1.0700, 2.7576, 1.3941, 1.7162, 1.3461, 2.7103, 0.8456, 0.9710, 0.6149, 1.3543,
            1.0006, 0.9705, 1.1978, 1.8814, 1.5150),
      b = c(1.1640, -0.9606, 0.5895, -0.1602, 0.1801, 0.7633, 1.0025, 1.1542, 0.5137, 0.0823,
            -0.5269, 1.1872, 1.1288, 0.0991, -0.1955),
      g = c(0.2637, 0.1888, 0.1845, 0.3310, 0.2532, 0.1545, 0.1319, 0.2345, 0.1529, 0.1468,
            0.1708, 0.1099, 0.1766, 0.2699, 0.1671)
    ),
    cats = 2, model = "3PLM"
  )
  r15 <- rbind(c(1, rep(0, 14)))
  map100 <- est_score(x15, r15, D = 1.702, method = "MAP", max.iter = 100)
  map101 <- est_score(x15, r15, D = 1.702, method = "MAP", max.iter = 101)
  expect_equal(map100$est.theta, -1.398218, tolerance = 1e-3)
  expect_equal(map101$est.theta, -1.398218, tolerance = 1e-3)
  expect_equal(map100$se.theta, 0.6907323, tolerance = 1e-3)
})

test_that("ML estimates are local maxima of the log-likelihood for every response pattern", {
  # every response pattern of 10 items
  pat <- as.matrix(expand.grid(rep(list(0:1), 10)))
  ml <- est_score(x_drm, pat, D = 1, method = "ML")

  # log-likelihood of each pattern at its own theta value, written out independently
  ll_vec <- function(th) {
    p <- sapply(seq_len(10), function(j) {
      x_drm$par.3[j] + (1 - x_drm$par.3[j]) / (1 + exp(-x_drm$par.1[j] * (th - x_drm$par.2[j])))
    })
    rowSums(pat * log(p) + (1 - pat) * log(1 - p))
  }

  # patterns whose estimate lies inside the range cannot be improved by a small move
  inner <- ml$se.theta != 99.9999
  ll0 <- ll_vec(ml$est.theta)
  expect_true(all(ll0[inner] >= ll_vec(ml$est.theta - 1e-3)[inner] - 1e-9))
  expect_true(all(ll0[inner] >= ll_vec(ml$est.theta + 1e-3)[inner] - 1e-9))
})

test_that("INV.TCC SEs are the conditional SDs over the scores with an estimate when intpol = FALSE", {
  x20 <- x_full[1:20, ]
  res <- suppressWarnings(est_score(x20, rbind(rep(1, 20)), D = 1, method = "INV.TCC", intpol = FALSE))
  st <- res$score.table
  ok <- !is.na(st$est.theta)
  expect_true(any(!ok))
  expect_true(all(is.finite(st$se.theta[ok])))
  expect_true(all(is.na(st$se.theta[!ok])))

  # conditional SD over the scores that have an estimate
  lk <- lwrc(x20, theta = st$est.theta[ok])[ok, , drop = FALSE]
  lk <- t(t(lk) / colSums(lk))
  mu <- colSums(lk * st$est.theta[ok])
  expect_equal(unname(st$se.theta[ok]), unname(sqrt(colSums(lk * st$est.theta[ok]^2) - mu^2)), tolerance = 1e-12)

  # the SEs do not zigzag between neighboring scores
  expect_true(all(abs(diff(st$se.theta[ok])) < 0.5))
})

test_that("INV.TCC SEs are rescaled when the interpolation range is rejected", {
  x20 <- x_full[1:20, ]
  d20 <- rbind(rep(1, 20))
  expect_warning(
    res <- est_score(x20, d20, D = 1, method = "INV.TCC", range.tcc = c(-1, 7)),
    "lower bound"
  )
  st <- res$score.table
  ok <- !is.na(st$est.theta)
  expect_true(any(!ok))
  lk <- lwrc(x20, theta = st$est.theta[ok])[ok, , drop = FALSE]
  lk <- t(t(lk) / colSums(lk))
  mu <- colSums(lk * st$est.theta[ok])
  expect_equal(unname(st$se.theta[ok]), unname(sqrt(colSums(lk * st$est.theta[ok]^2) - mu^2)), tolerance = 1e-12)
})

test_that("examinees with all missing responses get NA under every pointwise method", {
  d <- rbind(rep(NA, 10), c(1, 0, 1, 0, 1, 0, 1, 0, 1, 0))
  for (m in c("ML", "MLF", "WL", "MAP", "EAP")) {
    expect_warning(res <- est_score(x_drm, d, D = 1, method = m), "all missing")
    expect_true(all(is.na(res[1, ])))
    expect_false(anyNA(res[2, ]))

    # the other examinee gets the result of scoring alone
    alone <- est_score(x_drm, d[2, , drop = FALSE], D = 1, method = m)
    expect_equal(res[2, ], alone[1, ], ignore_attr = TRUE)
  }

  # a data frame is handled in the same way for MLF
  expect_warning(res_df <- est_score(x_drm, as.data.frame(d), D = 1, method = "MLF"), "all missing")
  expect_true(all(is.na(res_df[1, ])))
})

test_that("EAP gives finite estimates for a very long test", {
  # 1500 items: the product of the item probabilities underflows to 0 at every node
  x_long <- x_full[rep(1:50, 30), ]
  x_long$id <- paste0("item", seq_len(nrow(x_long)))
  set.seed(31)
  d_long <- simdat(x_long, c(-1, 0.5), D = 1)

  eap <- est_score(x_long, d_long, D = 1, method = "EAP")
  ml <- est_score(x_long, d_long, D = 1, method = "ML")
  expect_false(anyNA(eap))
  expect_true(all(is.finite(eap$est.theta)))
  expect_true(all(eap$se.theta > 0))

  # with this many items the EAP estimate is close to the ML estimate
  expect_lt(max(abs(eap$est.theta - ml$est.theta)), 0.1)
})

test_that("est_score() has a method for est_item objects that matches the default method", {
  # calibrate the pretest items of a small test with the fixed abilities
  utils::capture.output(
    ei <- est_item(x = x_drm[1:6, ], data = resp_drm[1:120, 1:6], score = theta_drm[1:120], D = 1.702)
  )
  expect_s3_class(ei, "est_item")

  for (m in c("ML", "WL", "EAP", "EAP.SUM", "INV.TCC")) {
    expect_identical(
      est_score(ei, method = m),
      est_score(ei$par.est, ei$data, D = ei$scale.D, method = m)
    )
  }

  # data and D come from the object
  expect_error(est_score(ei, data = resp_drm[1:3, 1:6], method = "ML"), "cannot be supplied")
  expect_error(est_score(ei, D = 1, method = "ML"), "cannot be supplied")
})

test_that("est_score() for an est_item object skips the items excluded for having no responses", {
  # the third item has no responses and is excluded from the calibration
  d_del <- resp_drm[1:120, 1:6]
  d_del[, 3] <- NA
  expect_warning(
    ei_del <- est_item(
      x = x_drm[1:6, ], data = d_del, score = theta_drm[1:120], D = 1.702, verbose = FALSE
    ),
    "excluded"
  )
  expect_equal(ei_del$deleted.item, 3)

  for (m in c("ML", "MAP", "EAP", "EAP.SUM")) {
    expect_identical(
      est_score(ei_del, method = m),
      est_score(ei_del$par.est[-3, ], ei_del$data, D = ei_del$scale.D, method = m)
    )
  }
})


# ==============================================================================
# 8. VALUES AGAINST INDEPENDENT COMPUTATIONS
# ==============================================================================

# mixed-format test: 3 x 3PLM, 2 x GPCM, 1 x GRM
x_mix <- shape_df(
  par.drm = list(a = c(1.0, 1.2, 0.9), b = c(-1.0, 0, 1.0), g = c(0.2, 0.15, 0.1)),
  par.prm = list(a = c(1.0, 1.2, 0.8), d = list(c(-1.0, 0.0, 1.0), c(-0.5, 0.5, 1.5), c(-1, 0.3))),
  cats = c(2, 2, 2, 4, 4, 3), model = c(rep("3PLM", 3), "GPCM", "GRM", "GPCM")
)

# category probabilities of one item written out independently of the package
cat_prob <- function(th, it, D = 1) {
  a <- it$par.1
  if (it$model == "3PLM") {
    p <- it$par.3 + (1 - it$par.3) / (1 + exp(-D * a * (th - it$par.2)))
    return(cbind(1 - p, p))
  }
  d <- as.numeric(unlist(it[paste0("par.", 2:it$cats)]))
  if (it$model == "GRM") {
    ps <- matrix(sapply(d, function(b) plogis(D * a * (th - b))), nrow = length(th))
    return(cbind(1, ps) - cbind(ps, 0))
  }
  z <- matrix(sapply(d, function(b) D * a * (th - b)), nrow = length(th))
  e <- exp(t(apply(cbind(0, z), 1, cumsum)))
  e / rowSums(e)
}

# log-likelihood of one response pattern
loglik_ref <- function(th, x, resp, D = 1) {
  out <- 0
  for (j in seq_len(nrow(x))) {
    if (!is.na(resp[j])) out <- out + log(cat_prob(th, x[j, ], D)[, resp[j] + 1])
  }
  out
}

# test information by numerical differentiation of the category probabilities
info_ref <- function(th, x, h = 1e-4) {
  sum(sapply(seq_len(nrow(x)), function(j) {
    dp <- (cat_prob(th + h, x[j, ]) - cat_prob(th - h, x[j, ])) / (2 * h)
    sum(dp^2 / cat_prob(th, x[j, ]))
  }))
}

set.seed(5)
resp_mix <- simdat(x_mix, rnorm(8), D = 1)
resp_mix[3, c(2, 5)] <- NA

test_that("ML and MAP agree with optimize() on a mixed-format test", {
  ml <- est_score(x_mix, resp_mix, D = 1, method = "ML", tol = 1e-9)
  map <- est_score(x_mix, resp_mix, D = 1, method = "MAP", norm.prior = c(0.3, 1.2), tol = 1e-9)
  for (i in seq_len(nrow(resp_mix))) {
    r <- resp_mix[i, ]
    ok <- !is.na(r)
    ref_ml <- optimize(function(t) loglik_ref(t, x_mix, r), c(-5, 5), maximum = TRUE, tol = 1e-10)$maximum
    ref_map <- optimize(function(t) loglik_ref(t, x_mix, r) + dnorm(t, 0.3, 1.2, log = TRUE),
                        c(-5, 5), maximum = TRUE, tol = 1e-10)$maximum
    if (ml$se.theta[i] != 99.9999) {
      expect_equal(ml$est.theta[i], ref_ml, tolerance = 1e-6)
      expect_equal(ml$se.theta[i], 1 / sqrt(info_ref(ref_ml, x_mix[ok, ])), tolerance = 1e-6)
    }
    expect_equal(map$est.theta[i], ref_map, tolerance = 1e-6)
    expect_equal(map$se.theta[i], 1 / sqrt(info_ref(ref_map, x_mix[ok, ]) + 1 / 1.2^2), tolerance = 1e-6)
  }
})

test_that("EAP and its posterior SD agree with direct quadrature", {
  w <- gen.weight(41, dist = "norm", mu = 0.3, sigma = 1.2)
  eap <- est_score(x_mix, resp_mix, D = 1, method = "EAP", norm.prior = c(0.3, 1.2), nquad = 41)
  for (i in seq_len(nrow(resp_mix))) {
    post <- exp(loglik_ref(w[, 1], x_mix, resp_mix[i, ])) * w[, 2]
    post <- post / sum(post)
    m <- sum(post * w[, 1])
    expect_equal(eap$est.theta[i], m, tolerance = 1e-12)
    expect_equal(eap$se.theta[i], sqrt(sum(post * w[, 1]^2) - m^2), tolerance = 1e-12)
  }
})

test_that("ML returns the range limits with SE 99.9999 for perfect patterns", {
  res <- est_score(x_drm, rbind(rep(1, 10), rep(0, 10)), D = 1, method = "ML", range = c(-4, 4))
  expect_equal(res$est.theta, c(4, -4))
  expect_equal(res$se.theta, c(99.9999, 99.9999))
})

test_that("a missing code gives the same result as NA", {
  d9 <- rbind(c(1, 0, 1, -9, 1, 0, 1, 0, 1, 0))
  dna <- d9
  dna[dna == -9] <- NA
  expect_identical(
    est_score(x_drm, d9, D = 1, missing = -9),
    est_score(x_drm, dna, D = 1)
  )
})

test_that("a single response vector is scored like a one-row matrix", {
  r <- c(1, 0, 1, 0, 1, 0, 1, 0, 1, 0)
  expect_equal(est_score(x_drm, r, D = 1), est_score(x_drm, rbind(r), D = 1))
})

test_that("lwrc() matches enumeration of all response patterns", {
  grid <- as.matrix(expand.grid(lapply(x_mix$cats, function(k) 0:(k - 1))))
  th <- c(-2.3, 0, 1.1)
  ref <- sapply(th, function(t) {
    pr <- apply(grid, 1, function(r) exp(loglik_ref(t, x_mix, r)))
    tapply(pr, rowSums(grid), sum)
  })
  expect_equal(unname(lwrc(x_mix, theta = th, D = 1)), unname(ref), tolerance = 1e-12)

  # a single theta value returns a one-column matrix
  expect_equal(unname(lwrc(x_mix, theta = 0)[, 1]), unname(lwrc(x_mix, theta = c(0, 1))[, 1]))
})

test_that("EAP.SUM matches the brute-force summed-score posterior", {
  grid <- as.matrix(expand.grid(lapply(x_mix$cats, function(k) 0:(k - 1))))
  w <- gen.weight(41, dist = "norm", mu = 0.2, sigma = 1.1)
  lss <- sapply(w[, 1], function(t) {
    pr <- apply(grid, 1, function(r) exp(loglik_ref(t, x_mix, r)))
    tapply(pr, rowSums(grid), sum)
  })
  post <- t(t(lss) * w[, 2])
  post <- post / rowSums(post)
  m <- c(post %*% w[, 1])
  res <- est_score(x_mix, grid[1:3, ], D = 1, method = "EAP.SUM", norm.prior = c(0.2, 1.1))
  expect_equal(res$score.table$est.theta, m, tolerance = 1e-12)
  expect_equal(res$score.table$se.theta, sqrt(c(post %*% w[, 1]^2) - m^2), tolerance = 1e-12)
})

test_that("INV.TCC matches uniroot() on the TCC and interpolates below the guessing sum", {
  tcc <- function(t) {
    sum(sapply(seq_len(nrow(x_mix)), function(j) sum((0:(x_mix$cats[j] - 1)) * cat_prob(t, x_mix[j, ]))))
  }
  res <- est_score(x_mix, rbind(rep(0, 6)), D = 1, method = "INV.TCC", tol = 1e-10)$score.table
  ref <- sapply(1:10, function(s) uniroot(function(t) tcc(t) - s, c(-20, 20), tol = 1e-12)$root)
  expect_equal(res$est.theta[2:11], ref, tolerance = 1e-8)
  expect_equal(res$est.theta[c(1, 12)], c(-7, 7))
})

test_that("llike_score() matches the reference log-likelihood", {
  th <- seq(-2, 2, 1)
  ll <- llike_score(x_mix, resp_mix[1:4, ], th, D = 1, method = "ML")
  ref <- sapply(1:4, function(i) loglik_ref(th, x_mix, resp_mix[i, ]))
  expect_equal(unname(as.matrix(ll)), ref, tolerance = 1e-12)

  # the MAP values add the log prior density
  llp <- llike_score(x_mix, resp_mix[1:4, ], th, D = 1, method = "MAP", norm.prior = c(0.2, 1.1))
  expect_equal(unname(as.matrix(llp)), ref + dnorm(th, 0.2, 1.1, log = TRUE), tolerance = 1e-12)
})

test_that("bisection() finds the root within tol", {
  f <- function(t, p) p - drm(theta = t, a = 1, b = 0.2, D = 1)
  r <- bisection(f, p = 0.2, lb = -10, ub = 10, tol = 1e-8)
  expect_lt(abs(r$root - (0.2 + qlogis(0.2))), 1e-8)
})
