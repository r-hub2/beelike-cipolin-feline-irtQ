# Expected item scores of trace() use only the categories of each item

gpcm_score_ref <- function(theta, a, d, D) {
  lz <- t(apply(outer(theta, c(0, d), "-") * D * a, 1, cumsum))
  e <- exp(lz - apply(lz, 1, max))
  as.vector((e / rowSums(e)) %*% (seq_along(c(0, d)) - 1))
}

grm_score_ref <- function(theta, a, d, D) {
  rowSums(sapply(d, function(dk) 1 / (1 + exp(-D * a * (theta - dk)))))
}

test_that("trace() expected scores are exact when items have fewer categories than the widest item", {
  th <- c(-2, -1, 0, 1, 2)
  x <- shape_df(
    par.drm = list(a = 1.1, b = 0.2, g = 0.2),
    par.prm = list(a = c(1.2, 0.9, 1.3), d = list(c(-1, 0, 1, 2), c(-0.5, 0.5), c(-0.3, 0.8))),
    cats = c(2, 5, 3, 3), model = c("3PLM", "GRM", "GPCM", "GRM")
  )
  x <- irtQ:::confirm_df(x)
  r <- irtQ:::trace(irtQ:::breakdown(x), theta = th, D = 1.702)
  expect_equal(r$icc[, 2], grm_score_ref(th, 1.2, c(-1, 0, 1, 2), 1.702), tolerance = 1e-12)
  expect_equal(r$icc[, 3], gpcm_score_ref(th, 0.9, c(-0.5, 0.5), 1.702), tolerance = 1e-12)
  expect_equal(r$icc[, 4], grm_score_ref(th, 1.3, c(-0.3, 0.8), 1.702), tolerance = 1e-12)
  expect_equal(r$tcc, rowSums(r$icc), tolerance = 1e-12)
})

test_that("trace() works when unused parameter columns are left over", {
  th <- c(-1, 0, 1)
  x <- irtQ:::confirm_df(shape_df(
    par.prm = list(a = c(1.2, 0.9), d = list(c(-1, 0, 1, 2), c(-0.5, 0.5))),
    cats = c(5, 3), model = c("GRM", "GPCM")
  ))
  # removing the five-category item keeps the par.5 column
  r <- irtQ:::trace(irtQ:::breakdown(x[2, ]), theta = th, D = 1)
  expect_equal(as.vector(r$icc), gpcm_score_ref(th, 0.9, c(-0.5, 0.5), 1), tolerance = 1e-12)
  expect_equal(r$tcc, as.vector(r$icc), tolerance = 1e-12)
  expect_equal(unname(r$prob.cats[[1]][, 1:3]),
    unname(irtQ::prm(th, 0.9, c(-0.5, 0.5), 1, "GPCM")),
    tolerance = 1e-9
  )
})
