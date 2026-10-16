# messages of est_mg() about groups without common items and about items without responses

M <- make_mg_data(seed = 11)
gn <- c("G1", "G2", "G3")
q <- function(...) est_mg(..., verbose = FALSE, MaxE = 2, se = FALSE)

test_that("the messages name the right groups and items", {
  x4 <- M$x
  x4[[3]]$id <- paste0("Z", 1:6)
  expect_error(q(x = x4, data = M$data, group.name = gn, free.group = 2:3), "G3 group\\(s\\) .* not share")
  expect_error(q(x = x4, data = M$data, group.name = gn, free.group = 1:3, fipc = TRUE, fix.loc = list(1:4, NULL, NULL)),
               "^G3 group\\(s\\)")
  d6 <- M$data
  d6[[3]][, 4] <- NA
  expect_error(q(x = M$x, data = d6, group.name = gn, free.group = 1:3, fipc = TRUE, fix.loc = list(1:4, 1:4, NULL)),
               "^E2 has/have no item response")
})
