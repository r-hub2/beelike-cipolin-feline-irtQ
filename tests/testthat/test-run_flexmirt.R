# run_flexmirt() needs Windows and an installation of flexMIRT, so only the failure paths are tested

test_that("run_flexmirt() stops on systems other than Windows", {
  skip_on_os("windows")
  expect_error(run_flexmirt("model.flexmirt"), "requires Windows")
  expect_error(run_flexmirt("model.flexmirt", dir.flex = tempdir()), "requires Windows")
})

test_that("run_flexmirt() stops when the folder of flexMIRT does not exist", {
  skip_if_not(.Platform$OS.type == "windows")
  missing_dir <- file.path(tempdir(), "no_such_flexmirt_folder")
  expect_error(run_flexmirt("model.flexmirt", dir.flex = missing_dir), "does not exist")
})
