# input checks of est_mg()

M <- make_mg_data(seed = 11)
gn <- c("G1", "G2", "G3")
q <- function(...) est_mg(..., verbose = FALSE, MaxE = 2, se = FALSE)

test_that("free.group must name existing groups and leave a reference group without FIPC", {
  skip_on_cran()
  expect_error(q(x = M$x, data = M$data, group.name = gn, free.group = c(2, 5)), "group indices between 1 and 3")
  expect_error(q(x = M$x, data = M$data, group.name = gn, free.group = c("G2", "g3")), "g3 in 'free.group'")
  expect_error(q(x = M$x, data = M$data, group.name = gn, free.group = 1:3), "At least one group")
  # a repeated free group is counted once
  f1 <- suppressWarnings(q(x = M$x, data = M$data, group.name = gn, free.group = 2))
  f2 <- suppressWarnings(q(x = M$x, data = M$data, group.name = gn, free.group = c(2, 2)))
  expect_identical(f2$npar.est, f1$npar.est)
  expect_identical(f2$par.est, f1$par.est)
})

test_that("group.name, x, model and item.id must match the groups in data", {
  expect_error(q(x = M$x, data = M$data, group.name = gn[1:2], free.group = 2), "unique names")
  expect_error(q(x = M$x, data = M$data, group.name = c("G1", "G1", "G3"), free.group = 3), "unique names")
  d2 <- M$data
  d2[[2]] <- d2[[2]][, -1]
  expect_error(q(x = M$x, data = d2, group.name = gn, free.group = 2:3), "differ in group\\(s\\) 2")
  md <- lapply(M$x, function(z) z$model)
  ct <- lapply(M$x, function(z) z$cats)
  id <- lapply(M$x, function(z) z$id)
  expect_error(q(data = M$data, group.name = gn, free.group = 2:3, model = md[1:2], cats = ct, item.id = id), "'model' must be")
  expect_error(q(x = M$x, data = M$data, group.name = gn, free.group = 2:3, item.id = id[1:2]), "'item.id' must be")
})

test_that("items with the same ID must have the same categories and model", {
  x2 <- M$x
  x2[[2]]$model[1] <- "3PLM"
  expect_error(q(x = x2, data = M$data, group.name = gn, free.group = 2:3), "Check item\\(s\\) C1")
  x3 <- M$x
  x3[[2]]$cats[4] <- 4
  expect_error(q(x = x3, data = M$data, group.name = gn, free.group = 2:3), "Check item\\(s\\) C4")
  x4 <- M$x
  x4[[1]]$id[2] <- "C1"
  expect_error(q(x = x4, data = M$data, group.name = gn, free.group = 2:3), "unique within each group")
})

test_that("responses must be scored 0, ..., cats - 1", {
  d3 <- M$data
  d3[[2]][1, 1] <- 2
  expect_error(q(x = M$x, data = d3, group.name = gn, free.group = 2:3), "scored as 0, 1")
  d4 <- M$data
  d4[[1]][1, 1] <- 0.5
  expect_error(q(x = M$x, data = d4, group.name = gn, free.group = 2:3), "scored as 0, 1")
})

test_that("fix.loc, fix.id and fipc.method are checked in MG-FIPC", {
  qf <- function(...) q(x = M$x, data = M$data, group.name = gn, free.group = 1:3, fipc = TRUE, ...)
  expect_error(qf(fix.loc = list(c(1:4, 20), 1:4, NULL)), "Check group\\(s\\) 1")
  expect_error(qf(fix.loc = list(c(1, 1, 2), 1:4, NULL)), "Check group\\(s\\) 1")
  expect_error(qf(fix.loc = list(1:4, 1:4)), "one element for each of the 3 groups")
  # a vector with one value per group is not taken as one position per group
  expect_error(qf(fix.loc = 1:3), "one element for each of the 3 groups")
  expect_error(qf(fix.loc = list(integer(0), NULL, NULL)), "at least one fixed item")
  expect_error(qf(fix.id = c("C1", "ZZ")), "not found in the item IDs: ZZ")
  expect_error(qf(fix.loc = list(1:4, 1:4, NULL), fipc.method = "oem"), "'fipc.method' must be")
})

test_that("missing values are recoded in the group data as well", {
  d5 <- M$data
  d5[[1]][1:5, 2] <- 99
  f <- suppressWarnings(q(x = M$x, data = d5, group.name = gn, free.group = 2:3, missing = 99))
  expect_false(any(f$data$group$G1 == 99, na.rm = TRUE))
  expect_true(all(is.na(f$data$group$G1[1:5, 2])))
  # data frames are recoded as well, with and without FIPC
  d7 <- lapply(d5, as.data.frame)
  f1 <- suppressWarnings(q(x = M$x, data = d7, group.name = gn, free.group = 2:3, missing = 99))
  expect_true(all(is.na(f1$data$group$G1[1:5, 2])))
  f2 <- suppressWarnings(q(x = M$x, data = d5, group.name = gn, free.group = 1:3, fipc = TRUE, fix.id = c("C1", "C2"),
                           missing = 99))
  expect_true(all(is.na(f2$data$group$G1[1:5, 2])))
})

test_that("DRM items are accepted when x = NULL", {
  md <- lapply(M$x, function(z) ifelse(z$cats == 2, "DRM", z$model))
  ct <- lapply(M$x, function(z) z$cats)
  id <- lapply(M$x, function(z) z$id)
  f <- suppressWarnings(q(data = M$data, group.name = gn, free.group = 2:3, model = md, cats = ct, item.id = id))
  expect_true(all(f$par.est$overall$model[f$par.est$overall$cats == 2] == "3PLM"))
  # the warning about the DRM items is given once, not once per group
  r <- with_warnings(q(data = M$data, group.name = gn, free.group = 2:3, model = md, cats = ct, item.id = id))
  expect_equal(sum(grepl("All 'DRM' items are treated as '3PLM' items", r$warnings)), 1L)
  # a single warning is also given when x holds the DRM items
  xd <- lapply(M$x, function(z) {
    z$model[z$cats == 2] <- "DRM"
    z
  })
  r2 <- with_warnings(q(x = xd, data = M$data, group.name = gn, free.group = 2:3))
  expect_equal(sum(grepl("All 'DRM' items are treated as '3PLM' items", r2$warnings)), 1L)
})

test_that("cats = NULL gives two categories when all models are dichotomous", {
  k <- lapply(M$x, function(z) z$cats == 2)
  dd <- lapply(1:2, function(g) M$data[[g]][, k[[g]]])
  id <- lapply(1:2, function(g) M$x[[g]]$id[k[[g]]])
  f0 <- suppressWarnings(q(data = dd, group.name = gn[1:2], free.group = 2, model = list("2PLM", "2PLM"),
                           cats = list(2, 2), item.id = id))
  f1 <- suppressWarnings(q(data = dd, group.name = gn[1:2], free.group = 2, model = list("2PLM", "2PLM"), item.id = id))
  expect_identical(f1$par.est, f0$par.est)
})
