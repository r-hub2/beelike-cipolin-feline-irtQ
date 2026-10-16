#' Run flexMIRT from Within R
#'
#' This function runs flexMIRT (Cai, 2017) from within R by executing a model
#' specified in a flexMIRT syntax file (i.e., *.flexmirt). flexMIRT runs on
#' Windows only, so this function works only on Windows. To use this function,
#' the flexMIRT software must be installed on your system. This interface is
#' especially useful for conducting simulation studies or automating batch
#' analyses involving flexMIRT.
#'
#' @param file.syntax A single string or character vector specifying the path(s)
#'   to one or more flexMIRT syntax files (with extension *.flexmirt) to be
#'   executed. For example: `"C:/Users/Data/irtmodel.flexmirt"`.
#' @param dir.flex A character string giving the folder where flexMIRT is
#'   installed (e.g., `"C:/Program Files/flexMIRT3"`). If `NULL` (default),
#'   the function looks in `"C:/Program Files"` for a folder whose name
#'   contains "flexMIRT"; give the folder explicitly when there is no such
#'   folder or more than one.
#' @param show.output.on.console Logical. If `TRUE`, the output of the system
#'   command is printed to the R console. Default is `FALSE`. See
#'   [base::system()].
#' @param ... Additional arguments passed to [base::system()].
#'
#' @details When using a version of flexMIRT earlier than 3.6, the directory
#'   specified in `dir.flex` must contain the following six files:
#' \itemize{
#'   \item \code{WinFlexMIRT.exe}
#'   \item \code{FlexMIRT_x64.exe}
#'   \item \code{FlexMIRT_x86.exe}
#'   \item \code{vpg.dll}
#'   \item \code{vpg.licensing.client.dll}
#'   \item \code{vpg.licensing.dll}
#' }
#'
#'   For flexMIRT version 3.6 or later, the directory must include the following
#'   five files:
#' \itemize{
#'   \item \code{WinFlexMIRT.exe}
#'   \item \code{vpg.dll}
#'   \item \code{vpg.licensing.client.dll}
#'   \item \code{vpg.licensing.dll}
#'   \item \code{VPGLicenseClientNet.dll}
#' }
#'   along with a subdirectory named \code{Resources} that contains the
#'   following two files:
#' \itemize{
#'   \item \code{flexMIRT_x64_AVX.exe}
#'   \item \code{flexMIRT_x86_AVX.exe}
#' }
#'
#'   The syntax files "2PLM_example.flexmirt" and "3PLM_example.flexmirt" in
#'   the extdata folder of the package fit the 2PL and 3PL models to the
#'   responses of 2,844 examinees to 12 items in "g341-19.dat". The folder also
#'   contains the resulting "-prm.txt" and "-irt.txt" files.
#'
#' @return Called for its side effect. flexMIRT writes its output files (e.g.,
#'   "-prm.txt", "-irt.txt") to the folder of each syntax file, and they can be
#'   read with [irtQ::bring.flexmirt()]. The value returned to R is not meant
#'   to be used.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @references Cai, L. (2017). flexMIRT 3.5 Flexible multilevel multidimensional
#'   item analysis and test scoring (Computer software). Chapel Hill, NC: Vector
#'   Psychometric Group.
#'
#' @examples
#' # Examples below will run if the flexMIRT software is installed
#' # in the default directory "C:/Program Files/flexMIRT3".
#' # Otherwise, specify the directory where flexMIRT is installed
#' # using the 'dir.flex' argument.
#'
#' \dontrun{
#' # (1) Run a single syntax file
#' # Copy an example syntax file and its data file to a temporary folder
#' ext <- system.file("extdata", package = "irtQ")
#' file.copy(file.path(ext, c("2PLM_example.flexmirt", "g341-19.dat")),
#'           tempdir())
#' file.syntax <- file.path(tempdir(), "2PLM_example.flexmirt")
#'
#' # Run flexMIRT to estimate item parameters for the 2PL model
#' run_flexmirt(file.syntax = file.syntax, dir.flex = NULL,
#'              show.output.on.console = TRUE)
#'
#' # Read the item parameter file written by flexMIRT
#' out.file <- file.path(tempdir(), "2PLM_example-prm.txt")
#' bring.flexmirt(out.file, type = "par")
#'
#' # (2) Run multiple syntax files
#' # Copy two example syntax files and the data file to a temporary folder
#' file.copy(file.path(ext, c("3PLM_example.flexmirt", "g341-19.dat")),
#'           tempdir())
#' file.syntax1 <- file.path(tempdir(), "2PLM_example.flexmirt")
#' file.syntax2 <- file.path(tempdir(), "3PLM_example.flexmirt")
#'
#' # Run flexMIRT to estimate item parameters for both models
#' run_flexmirt(file.syntax = c(file.syntax1, file.syntax2), dir.flex = NULL,
#'              show.output.on.console = FALSE)
#'
#' # Read the item parameter files written by flexMIRT
#' out.file1 <- file.path(tempdir(), "2PLM_example-prm.txt")
#' out.file2 <- file.path(tempdir(), "3PLM_example-prm.txt")
#' bring.flexmirt(out.file1, type = "par")
#' bring.flexmirt(out.file2, type = "par")
#' }
#'
#' @export
run_flexmirt <- function(file.syntax,
                         dir.flex = NULL,
                         show.output.on.console = FALSE,
                         ...) {

  # stop unless R runs on Windows, where flexMIRT runs
  if (.Platform$OS.type != "windows") {
    stop("run_flexmirt() requires Windows, where flexMIRT runs.", call. = FALSE)
  }

  if (is.null(dir.flex)) {
    # find the folders in the default directory whose names contain flexMIRT
    flex_exist <- dir(path = "C:/Program Files", pattern = "flexMIRT", full.names = TRUE)

    # stop when no folder of flexMIRT is found in the default directory
    if (length(flex_exist) == 0) {
      stop("No folder with 'flexMIRT' in its name was found in 'C:/Program Files'. ",
        "Provide the folder where flexMIRT is installed in 'dir.flex'.",
        call. = FALSE
      )
    }

    # stop when several folders of flexMIRT are found in the default directory
    if (length(flex_exist) > 1) {
      stop("More than one folder with 'flexMIRT' in its name was found in 'C:/Program Files': ",
        paste0("'", basename(flex_exist), "'", collapse = ", "), ". ",
        "Provide the folder to use in 'dir.flex'.",
        call. = FALSE
      )
    }

    # a path of flexMIRT execution file
    file.flex <- file.path(flex_exist, "WinFlexMIRT.exe")
  } else {
    # check if the provided directory where flexMIRT is installed exists
    flex_exist <- dir.exists(paths = dir.flex)

    # stop when the provided directory where flexMIRT is installed does not exist
    if (!flex_exist) {
      stop("The provided directory '", dir.flex, "' does not exist. ",
        "Provide the full path of the folder where flexMIRT is installed.",
        call. = FALSE
      )
    }

    # a path of flexMIRT execution file
    file.flex <- file.path(dir.flex, "WinFlexMIRT.exe")
  }

  # run flexMIRT
  # when a single syntax files are provided
  if (length(file.syntax) == 1) {
    system(
      command = paste("\"", file.flex, "\"", " -r ", "\"", file.syntax, "\"", sep = ""),
      show.output.on.console = show.output.on.console, ...
    )
  }

  # when multiple syntax files are provided
  if (length(file.syntax) > 1) {
    # start the progress bar
    pb <- utils::txtProgressBar(min = 0, max = 100, initial = 0, width = 50, style = 3)
    for (r in 1:length(file.syntax)) {
      system(
        command = paste("\"", file.flex, "\"", " -r ", "\"", file.syntax[r], "\"", sep = ""),
        show.output.on.console = show.output.on.console, ...
      )

      # modify the progress bar
      info <- sprintf("Parameter Estimation Progress: %d%%", round((r / length(file.syntax)) * 100))
      utils::setTxtProgressBar(pb, value = r / length(file.syntax) * 100, title = info, label = paste0("Estimation ", r, " is done"))
    }

    # closing the progress Bar
    close(pb)
  }
}
