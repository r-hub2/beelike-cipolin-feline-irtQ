#' Write a "-prm.txt" File for flexMIRT
#'
#' @description This function writes a flexMIRT-compatible "-prm.txt" file (Cai,
#'   2017). It supports unidimensional IRT models only. This function was
#'   developed by modifying `read.flexmirt()` from Pritikin & Falk (2020).
#'
#' @param x A data frame of item metadata (e.g., item parameters, number of
#'   categories, model types) for a single group, or a list of such data frames
#'   for multiple groups. See [irtQ::est_irt()] or [irtQ::simdat()] for item
#'   metadata format. You can also create metadata using [irtQ::shape_df()].
#' @param file A character string giving the path of the file to be written.
#'   It must be supplied. flexMIRT parameter files usually end in "-prm.txt".
#' @param norm.pop A numeric vector of length two giving the mean and the
#'   variance of the normal ability distribution of a single group, or a list
#'   of such vectors, one per group (e.g.,
#'   `norm.pop = list(c(0, 1), c(0, 0.64), c(0.5, 1.44))` for three
#'   groups). The second value is written to the variance field of the file.
#'   If `mgroup = TRUE` and a single vector is given, it is used for all
#'   groups. Default is `c(0, 1)`.
#' @param rePar Logical. Indicates how the parameters in `x` are given. If
#'   `TRUE` (default), `x` is taken to hold flexMIRT slope-intercept values:
#'   slopes, intercepts, and logits of the guessing parameters, and for GPCM
#'   items the scoring function and intercept values of the nominal model;
#'   they are written as they are. If `FALSE`, `x` is taken to hold
#'   \pkg{irtQ} item metadata (slopes, difficulty or step parameters, and
#'   guessing parameters), which is converted to the flexMIRT
#'   parameterization. Item metadata from \pkg{irtQ} (for example, from
#'   [irtQ::shape_df()] or [irtQ::est_irt()]) must be written with
#'   `rePar = FALSE`. Note that the meaning is the reverse of the `rePar`
#'   argument of [irtQ::bring.flexmirt()], where `TRUE` converts the values
#'   of the file to \pkg{irtQ} metadata.
#' @param mgroup A logical value indicating whether the file includes multiple
#'   groups. Default is `FALSE`.
#' @param group.name A character vector of group names. If `NULL`, group names
#'   are automatically generated (e.g., "Group1", "Group2", ...).
#'
#' @details The function writes items of the models `"1PLM"`, `"2PLM"`,
#'   `"3PLM"`, `"DRM"`, `"GRM"`, and `"GPCM"`. Dichotomous items are written
#'   as flexMIRT 3PL (`"3PLM"` and `"DRM"`) or graded (`"1PLM"` and `"2PLM"`)
#'   items, so [irtQ::bring.flexmirt()] reads them back as `"3PLM"` or
#'   `"2PLM"` items. Every item and every group is written as one line, whatever
#'   the width of the console. The parameters are written without rescaling.
#'   flexMIRT uses the logistic metric (D = 1); when `x` was calibrated with
#'   another scaling constant, multiply the slopes by that constant before
#'   writing.
#'
#' @return Called for its side effect of writing the file to `file`. The value
#'   is returned invisibly and is not meant to be used.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @references Cai, L. (2017). flexMIRT 3.5 Flexible multilevel multidimensional
#'   item analysis and test scoring (Computer software). Chapel Hill, NC: Vector
#'   Psychometric Group.
#'
#'   Pritikin, J. N., & Falk, C. F. (2020). OpenMx: A modular research
#'   environment for item response theory method development. *Applied
#'   Psychological Measurement, 44*(7-8), 561-562.
#'   \doi{10.1177/0146621620929431}.
#'
#' @examples
#' \donttest{
#' ## 1. Create a "-prm.txt" file for a single group
#' ##    using the simulated CAT data
#' # 1-(1) Extract the item metadata
#' x <- simCAT_MX$item.prm
#'
#' # 1-(2) Set the name of the "-prm.txt" file
#' temp_prm <- file.path(tempdir(), "single_group_temp-prm.txt")
#'
#' # 1-(3) Write the "-prm.txt" file
#' write.flexmirt(x, file = temp_prm, norm.pop = c(0, 1), rePar = FALSE)
#'
#' # 1-(4) Read the file back and compare the item parameters
#' x_back <- bring.flexmirt(file = temp_prm, type = "par")$Group1$full_df
#' head(x_back)
#' max(abs(as.matrix(x[, 4:7]) - as.matrix(x_back[, 4:7])), na.rm = TRUE)
#'
#' ## 2. Create a "-prm.txt" file for multiple groups
#' ##    using simulated multi-group data
#' # 2-(1) Extract the item metadata
#' x <- simMG$item.prm
#'
#' # Set the name of the "-prm.txt" file
#' temp_prm <- file.path(tempdir(), "mg_group_temp-prm1.txt")
#'
#' # Write the "-prm.txt" file, giving the mean and variance of each group
#' write.flexmirt(x,
#'   file = temp_prm, norm.pop = list(c(0, 1), c(0.5, 0.64), c(-0.3, 1.69)),
#'   rePar = FALSE, mgroup = TRUE, group.name = c("GR1", "GR2", "GR3")
#' )
#'
#' # Or write the "-prm.txt" file so that
#' # all groups share the same ability distribution
#' # and group names are generated automatically
#' temp_prm <- file.path(tempdir(), "mg_group_temp-prm2.txt")
#' write.flexmirt(x,
#'   file = temp_prm, norm.pop = c(0, 1),
#'   rePar = FALSE, mgroup = TRUE, group.name = NULL
#' )
#' }
#'
#' @export
write.flexmirt <- function(x,
                           file = NULL,
                           norm.pop = c(0, 1),
                           rePar = TRUE,
                           mgroup = FALSE,
                           group.name = NULL) {

  if (!mgroup) {
    ngroup <- 1
    x_gr <- list(x)
    norm.pop <- list(norm.pop)
    if (is.null(group.name)) {
      group.name <- "Group1"
    }
  } else {
    x_gr <- x
    ngroup <- length(x)
    if (!is.list(norm.pop)) {
      norm.pop <- purrr::map(.x = 1:ngroup, ~ {
        norm.pop
      })
    }
    if (is.null(group.name)) {
      group.name <- paste0("Group", 1:ngroup)
    }
  }

  # stop when the destination of the file is not given
  if (is.null(file)) {
    stop("You must specify the destination of the file.", call. = FALSE)
  }

  # open a new file and close it when the function exits, even after an error
  prm_file <- file(file, open = "w")
  on.exit(close(prm_file), add = TRUE)

  for (g in 1:ngroup) {
    # give column names
    x <- data.frame(x_gr[[g]])
    colnames(x) <- c("id", "cats", "model", paste0("par.", 1:(ncol(x) - 3)))

    # select parameter columns
    param <- dplyr::select(x, dplyr::starts_with("par."))

    # write prm file one item by one item
    for (i in 1:nrow(x)) {
      # check metadata for each item
      id <- as.character(x$id[i])
      cats <- x$cats[i]
      model <- as.character(x$model[i])
      group <- g
      n.factor <- 1

      # (1) Dichotomous item
      if (model %in% c("3PLM", "DRM")) {
        a <- param[i, 1]
        if (!rePar) {
          c <- -a * param[i, 2]
          if (param[i, 3] == 0L) {
            param[i, 3] <- 1e-100
          }
          logitg <- log(param[i, 3] / (1 - param[i, 3]))
        } else {
          c <- param[i, 2]
          logitg <- param[i, 3]
        }
        a <- format(a, nsmall = 7)
        c <- format(c, nsmall = 7)
        logitg <- format(logitg, nsmall = 7)
        # write the item as one tab-separated line that is never wrapped
        writeLines(paste(c(1, id, group, n.factor, 1, cats, logitg, c, a), collapse = "\t"), con = prm_file)
      }

      if (model %in% c("1PLM", "2PLM")) {
        a <- param[i, 1]
        if (!rePar) {
          c <- -a * param[i, 2]
        } else {
          c <- param[i, 2]
        }
        a <- format(a, nsmall = 7)
        c <- format(c, nsmall = 7)
        # write the item as one tab-separated line that is never wrapped
        writeLines(paste(c(1, id, group, n.factor, 2, cats, c, a), collapse = "\t"), con = prm_file)
      }

      # (2) Polytomous item: GRM
      if (model %in% c("GRM")) {
        a <- param[i, 1]
        if (!rePar) {
          cs <- -a * as.numeric(param[i, 2:cats])
        } else {
          cs <- as.numeric(param[i, 2:cats])
        }
        a <- format(a, nsmall = 7)
        cs <- format(cs, nsmall = 7)
        # write the item as one tab-separated line that is never wrapped
        writeLines(paste(c(1, id, group, n.factor, 2, cats, cs, a), collapse = "\t"), con = prm_file)
      }

      # (3) Polytomous item: PCM and GPCM
      if (model %in% c("PCM", "GPCM")) {
        a <- param[i, 1]
        if (!rePar) {
          alpha <- c(1, rep(0, cats - 2))
          ds.new <- param[i, 2:cats]
          # mean of the step parameters over the cats - 1 steps
          b <- sum(ds.new) / (cats - 1)
          ds <- b - ds.new

          Tmat <- matrix(0, nrow = cats, ncol = (cats - 1))
          Tmat[, 1] <- 0:(cats - 1)
          # sine terms of the categories 2 to cats - 1 (none for two categories)
          for (k in seq_len(cats - 2) + 1) {
            for (j in seq_len(cats - 2) + 1) {
              Tmat[k, j] <- sin(pi * (j - 1) * (k - 1) / (cats - 1))
            }
          }

          c.vec <- rep(0, cats)
          c.vec[cats] <- -a * (cats - 1) * b
          # intercepts of the categories cats - 1 down to 2 (none for two categories)
          for (k in rev(seq_len(cats - 2) + 1)) {
            c.vec[k] <- c.vec[k + 1] - a * (as.numeric(ds[k]) - b)
          }

          # keep the matrix form when there is a single step
          gam <- as.numeric(solve(Tmat[-1, , drop = FALSE]) %*% c.vec[-1])
        } else {
          # all cats - 1 scoring function values
          alpha <- as.numeric(param[i, 2:cats])
          gam <- as.numeric(param[i, (1 + cats - 2 + 1 + 1):(1 + cats - 2 + 1 + 1 + cats - 2)])
        }
        a <- format(a, nsmall = 7)
        alpha <- format(alpha, nsmall = 7)
        gam <- format(gam, nsmall = 7)
        # write the item as one tab-separated line that is never wrapped
        writeLines(paste(c(1, id, group, n.factor, 3, cats, 0, alpha, a, 0, gam), collapse = "\t"), con = prm_file)
      }
    }

    # group information
    # write the group as one tab-separated line that is never wrapped
    group_line <- c(
      0, group.name[g], g, n.factor, 0,
      format(norm.pop[[g]][1], nsmall = 7), format(norm.pop[[g]][2], nsmall = 7)
    )
    writeLines(paste(group_line, collapse = "\t"), con = prm_file)
  }

  # return nothing; the file is closed by the exit handler
  invisible(NULL)
}
