# bring.flexmirt(), bring.bilog(), and bring.parscale() read the file type and the model code

flex_file <- system.file("extdata", "3PLM_example-prm.txt", package = "irtQ")
parscale_file <- system.file("extdata", "parscale_sample.PAR", package = "irtQ")

# write a small BILOG-MG parameter file with two items
write_bilog_par <- function() {
  file <- tempfile(fileext = ".PAR")
  field <- function(values) paste(sprintf("%10.5f", values), collapse = "")
  lines <- c(
    "header 1", "header 2", "header 3", "header 4",
    paste0(sprintf("%-8s", "ITEM0001"), sprintf("%-8s", "1"), field(c(0, 0, 1.2, 0, -0.5, 0, 0, 0, 0.2, 0))),
    paste0(sprintf("%-8s", "ITEM0002"), sprintf("%-8s", "1"), field(c(0, 0, 0.8, 0, 0.7, 0, 0, 0, 0.1, 0)))
  )
  writeLines(lines, file)
  file
}

test_that("the file type is 'par' when it is omitted and its case is ignored", {
  expect_identical(bring.flexmirt(flex_file), bring.flexmirt(flex_file, "par"))
  expect_identical(bring.flexmirt(flex_file, "PAR"), bring.flexmirt(flex_file, "par"))
  expect_identical(bring.parscale(parscale_file), bring.parscale(parscale_file, "par"))
  expect_identical(bring.parscale(parscale_file, "Par"), bring.parscale(parscale_file, "par"))
  bilog_file <- write_bilog_par()
  expect_identical(bring.bilog(bilog_file), bring.bilog(bilog_file, "par"))
  expect_identical(bring.bilog(bilog_file, "PAR"), bring.bilog(bilog_file, "par"))
  expect_equal(bring.bilog(bilog_file)$full_df$par.2, c(-0.5, 0.7))
})

test_that("an unknown file type stops with the choices", {
  expect_error(bring.flexmirt(flex_file, type = "prm"), "should be one of")
  expect_error(bring.parscale(parscale_file, type = "prm"), "should be one of")
  expect_error(bring.bilog(write_bilog_par(), type = "prm"), "should be one of")
})

test_that("bring.parscale() stops for the even model codes with polytomous items", {
  lines <- readLines(parscale_file)
  # the third line holds the model code in its fourth field
  expect_match(lines[3], "    3    0    1$")
  for (code in c(2, 4, 6, 8)) {
    even_file <- tempfile(fileext = ".PAR")
    lines[3] <- sub("    [0-9]    0    1$", sprintf("    %d    0    1", code), lines[3])
    writeLines(lines, even_file)
    expect_error(bring.parscale(even_file), "Model codes 2, 4, 6, and 8 of PARSCALE are not supported")
  }
})
