# item.id replaces the IDs in x in est_mg() and est_item()

M <- make_mg_data(seed = 11)
gn <- c("G1", "G2", "G3")
q <- function(...) est_mg(..., verbose = FALSE, MaxE = 3, se = FALSE)

test_that("est_mg() replaces the IDs in x with item.id", {
  skip_on_cran()
  ids2 <- lapply(M$x, function(z) paste0(z$id, "_n"))
  f <- suppressWarnings(q(x = M$x, data = M$data, group.name = gn, free.group = 2:3, item.id = ids2))
  expect_identical(f$par.est$group$G2$id, ids2[[2]])
  expect_identical(f$par.est$overall$id, unique(unlist(ids2)))
  # the new IDs only relabel the items
  f0 <- suppressWarnings(q(x = M$x, data = M$data, group.name = gn, free.group = 2:3))
  expect_identical(f$par.est$overall[, -1], f0$par.est$overall[, -1])
  # item.id with the IDs of x gives the run without item.id
  id <- lapply(M$x, function(z) z$id)
  f1 <- suppressWarnings(q(x = M$x, data = M$data, group.name = gn, free.group = 2:3, item.id = id))
  expect_identical(f1$par.est, f0$par.est)
  # one ID per item of each group is required
  id_short <- id
  id_short[[2]] <- id_short[[2]][-1]
  expect_error(q(x = M$x, data = M$data, group.name = gn, free.group = 2:3, item.id = id_short), "Check group\\(s\\) 2")
})

test_that("est_item() replaces the IDs in x with item.id", {
  set.seed(7)
  xi <- shape_df(par.drm = list(a = rlnorm(5, 0, .3), b = rnorm(5), g = NA), cats = 2, model = "2PLM")
  th <- rnorm(300)
  d <- simdat(x = xi, theta = th, D = 1)
  f <- est_item(x = xi, data = d, score = th, item.id = paste0("Q", 1:5), verbose = FALSE)
  expect_identical(f$par.est$id, paste0("Q", 1:5))
  f0 <- est_item(x = xi, data = d, score = th, verbose = FALSE)
  expect_identical(f$par.est[, -1], f0$par.est[, -1])
  expect_error(est_item(x = xi, data = d, score = th, item.id = "Q1", verbose = FALSE), "length of 'item.id'")
})
