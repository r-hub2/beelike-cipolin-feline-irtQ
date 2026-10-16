# warnings about items without a maximum likelihood estimate in est_mg() and est_item()

M <- make_mg_data(seed = 11)
gn <- c("G1", "G2", "G3")
q <- function(...) est_mg(..., verbose = FALSE, MaxE = 2, se = FALSE)

test_that("est_mg() flags items with an empty category in the pooled data", {
  skip_on_cran()
  d7 <- M$data
  d7[[1]][, 5] <- 1
  w <- with_warnings(q(x = M$x, data = d7, group.name = gn, free.group = 2:3))$warnings
  expect_true(any(grepl("Item\\(s\\) A1 have a score category with no responses", w)))
  # C1 is constant in G1 only; its pooled responses vary, so it is not flagged
  d8 <- M$data
  d8[[1]][, 1] <- 1
  w <- with_warnings(q(x = M$x, data = d8, group.name = gn, free.group = 2:3))$warnings
  expect_false(any(grepl("score category with no responses", w)))
  # in FIPC only the items to be estimated are checked
  w <- with_warnings(q(x = M$x, data = d7, group.name = gn, free.group = 1:3, fipc = TRUE,
                       fix.loc = list(1:5, 1:4, NULL)))$warnings
  expect_false(any(grepl("score category with no responses", w)))
  w <- with_warnings(q(x = M$x, data = d7, group.name = gn, free.group = 1:3, fipc = TRUE,
                       fix.loc = list(1:4, 1:4, NULL)))$warnings
  expect_true(any(grepl("Item\\(s\\) A1 have a score category with no responses", w)))
})

test_that("est_item() flags items with an empty category", {
  set.seed(7)
  xi <- shape_df(par.drm = list(a = rlnorm(5, 0, .3), b = rnorm(5), g = NA), cats = 2, model = "2PLM")
  th <- rnorm(400)
  d <- simdat(x = xi, theta = th, D = 1)
  d[, 4] <- 1
  expect_warning(est_item(x = xi, data = d, score = th, verbose = FALSE),
                 "Item\\(s\\) V4 have a score category with no responses")
  # items without any response are removed first and reported separately
  d[, 2] <- NA
  w <- with_warnings(est_item(x = xi, data = d, score = th, verbose = FALSE))$warnings
  expect_true(any(grepl("excluded in the item parameter estimation", w)))
  expect_true(any(grepl("Item\\(s\\) V4 have", w)))
})
