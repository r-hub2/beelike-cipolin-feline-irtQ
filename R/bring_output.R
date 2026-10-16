#' Import Item and Ability Parameters from IRT Software
#'
#' @description These functions import item and/or ability parameter estimates
#'   from the output files of BILOG-MG 3, PARSCALE 4, and flexMIRT.
#'   `bring.mirt()` imports the item parameters, together with the latent mean
#'   and variance, of a model fitted with the \pkg{mirt} R package.
#'
#' @param file The path to the output file that contains the item or ability
#'   parameter estimates.
#' @param type A character string giving the type of output file: `"par"` for
#'   an item parameter file or `"sco"` for an ability (score) file. The case is
#'   ignored. Default is `"par"`. Not used by `bring.mirt()`.
#' @param rePar Logical. Used only by `bring.flexmirt()` with `type = "par"`.
#'   If `TRUE` (default), the intercepts of dichotomous and graded response
#'   model (GRM) items are converted to difficulty (threshold) parameters and
#'   the logits of the guessing parameters to guessing parameters, so that
#'   `full_df` follows the \pkg{irtQ} item metadata format. This requires a
#'   unidimensional model. If `FALSE`, the slope-intercept values in the file
#'   are kept, and `full_df` cannot then be used directly by other
#'   \pkg{irtQ} functions.
#' @param rePar.gpc Logical. Used only by `bring.flexmirt()` with
#'   `type = "par"`. If `TRUE` (default), the nominal model parameters of
#'   generalized partial credit model (GPCM) items are converted to a slope
#'   and step parameters. The function stops when the file contains a nominal
#'   model item or a multidimensional GPCM item; use `rePar.gpc = FALSE` to
#'   read such a file, in which case the scoring function and intercept values
#'   are kept and nominal items get the model name `"NOMM"`.
#' @param n.factor An integer giving the number of latent traits (factors).
#'   Used only by `bring.flexmirt()` with `type = "sco"`. Default is 1.
#' @param x A fitted model object returned by [mirt::mirt()] for a single
#'   group and a unidimensional model. The \pkg{mirt} package, listed in
#'   Suggests, must be installed to use `bring.mirt()`.
#'
#' @details The [irtQ::bring.flexmirt()] function was developed by modifying
#'   the `read.flexmirt()` function (Pritikin & Falk, 2020). Similarly,
#'   [irtQ::bring.bilog()] and [irtQ::bring.parscale()] were based on modified
#'   versions of the `read.bilog()` and `read.parscale()` functions
#'   (Weeks, 2010), respectively.
#'
#'   The functions read the layout of the file, not its name. The usual files
#'   are the ".PAR" and ".SCO" files of BILOG-MG and PARSCALE and the
#'   "-prm.txt" and "-sco.txt" files of flexMIRT; other flexMIRT output files
#'   (e.g., "-irt.txt") cannot be read. For \pkg{mirt}, the fitted model
#'   object is passed directly.
#'
#'   [irtQ::bring.flexmirt()] can also read the files of a multidimensional
#'   model: item parameter files with `rePar = FALSE` (and
#'   `rePar.gpc = FALSE` when GPCM items are present), and ability files with
#'   `n.factor` set to the number of factors. The other functions of
#'   \pkg{irtQ} work only with unidimensional models.
#'
#'   For polytomous items, [irtQ::bring.flexmirt()],
#'   [irtQ::bring.parscale()], and [irtQ::bring.mirt()] import the parameters
#'   of the graded response model (GRM) and the generalized partial credit
#'   model (GPCM); a polytomous \pkg{mirt} item of type `"Rasch"` is imported
#'   as a GPCM item with a slope of 1. [irtQ::bring.bilog()] reads
#'   dichotomous items only. [irtQ::bring.parscale()] does not support the
#'   model codes 2, 4, 6, and 8 of PARSCALE for polytomous items and stops
#'   with an error when the file has them.
#'
#'   [irtQ::bring.mirt()] accepts the item types `"Rasch"`, `"2PL"`, and
#'   `"3PL"` for dichotomous items and `"Rasch"`, `"graded"`, `"gpcm"`, and
#'   `"gpcmIRT"` for polytomous items. A two-category item of a mixed-format
#'   test may also be fitted with `"graded"`, `"gpcm"`, or `"gpcmIRT"`. A
#'   model fitted to multiple groups is not supported, and the function stops
#'   with an error; extract one group first with `mirt::extract.group()`.
#'
#'   In `full_df`, [irtQ::bring.flexmirt()] labels dichotomous items `"3PLM"`
#'   or `"2PLM"`, while [irtQ::bring.bilog()], [irtQ::bring.parscale()], and
#'   [irtQ::bring.mirt()] label them `"DRM"`, with a guessing parameter of 0
#'   for items without guessing. For [irtQ::bring.mirt()], this includes
#'   two-category items fitted with the `"graded"`, `"gpcm"`, or `"gpcmIRT"`
#'   model. Change the model column if another model is wanted in
#'   [irtQ::est_irt()], which fits `"DRM"` items as 3PLM items. The parameters
#'   are kept on the scale of the program, and the files do not record the
#'   scaling constant; pass the value of `D` used in the calibration (e.g., 1
#'   for flexMIRT and \pkg{mirt}) to the other \pkg{irtQ} functions.
#'
#' @return
#'   - `bring.flexmirt()`: a named list with one element per group (`Group1`,
#'     `Group2`, ...). With `type = "par"`, each element is a list with `id`,
#'     `factor` (number of dimensions), `cats`, `model`, `par` (a matrix of
#'     item parameters), `full_df` (the item metadata), `mean` and `cov` (the
#'     latent mean vector and covariance matrix), and `emphist` (a data frame
#'     with `theta` and `weight`, or `NULL` when no empirical histogram was
#'     estimated). With `type = "sco"`, each element is a list with `theta`
#'     and `se`, plus `iter` (the number of iterations) or `imputation` (the
#'     imputation number) when the file contains that column. `se` is `NULL`
#'     for multidimensional models and for imputed values.
#'   - `bring.bilog()`: with `type = "par"`, a list with `par` and `full_df`;
#'     with `type = "sco"`, a list with one element per group code in the
#'     file, each a list with `theta`, `se`, and `full_df` (a data frame with
#'     `group`, `id`, `nitem`, `ncorrect`, `p`, `theta`, and `se`).
#'   - `bring.parscale()`: with `type = "par"`, a list with `id`, `factor`,
#'     `cats`, `model`, `par_parscale` (the parameters and standard errors as
#'     reported by PARSCALE), `par`, and `full_df`; with `type = "sco"`, a
#'     data frame with `id`, `theta`, and `se` (a value of 999 becomes `NA`).
#'   - `bring.mirt()`: a list with `id`, `factor`, `cats`, `model` (a list
#'     with the \pkg{mirt} item types, `old`, and the \pkg{irtQ} model names,
#'     `new`), `emphist` (a matrix of quadrature points and weights),
#'     `par_mirt` (the \pkg{mirt} parameters), `par`, `full_df`, and `ability`
#'     (the latent mean and variance).
#'
#' @note For an item parameter file, the element `full_df` of the result is
#'   the item metadata (item IDs, numbers of score categories, model names,
#'   and item parameters) used by the other functions of \pkg{irtQ}. See
#'   [irtQ::info()] or [irtQ::simdat()] for the metadata format.
#'
#'   For GPCM items, `full_df` contains the step parameters
#'   \eqn{b_v = \beta - \tau_v}, computed from the overall location
#'   \eqn{\beta} and the category parameters \eqn{\tau_v} reported by the
#'   program (PARSCALE) or from the nominal model intercepts (flexMIRT). See
#'   the **IRT Models** section in [irtQ::irtQ-package] for the
#'   parameterization of the GPCM.
#'
#' @section Sample Output Files:
#'
#'   Two sample item parameter files are included in the package: the
#'   PARSCALE 4 file "parscale_sample.PAR" and the flexMIRT file
#'   "flexmirt_sample-prm.txt". Both contain the estimates for the same
#'   mixed-format test of 55 items: 50 dichotomous items fit with the 3PL
#'   model and five items with five score categories fit with the GRM (items
#'   39, 40, and 53 to 55). The examples below import these files.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [irtQ::irtQ-package]
#'
#' @references Cai, L. (2017). flexMIRT 3.5 Flexible multilevel multidimensional
#'   item analysis and test scoring (Computer software). Chapel Hill, NC: Vector
#'   Psychometric Group.
#'
#'   Chalmers, R. P. (2012). mirt: A multidimensional item response theory
#'   package for the R environment. *Journal of Statistical Software, 48*(6),
#'   1-29. \doi{10.18637/jss.v048.i06}
#'
#'   Muraki, E., & Bock, R. D. (2003). PARSCALE 4: IRT item analysis and test
#'   scoring for rating scale data (Computer software). Chicago, IL: Scientific
#'   Software International. URL http://www.ssicentral.com
#'
#'   Pritikin, J. (2018). *rpf: Response Probability Functions*. R package
#'   version 0.59. https://CRAN.R-project.org/package=rpf.
#'
#'   Pritikin, J. N., & Falk, C. F. (2020). OpenMx: A modular research
#'   environment for item response theory method development. *Applied
#'   Psychological Measurement, 44*(7-8), 561-562.
#'   \doi{10.1177/0146621620929431}.
#'
#'   Weeks, J. P. (2010). plink: An R Package for Linking Mixed-Format Tests
#'   Using IRT-Based Methods. *Journal of Statistical Software, 35*(12), 1-33.
#'   \doi{10.18637/jss.v035.i12}.
#'
#'   Zimowski, M. F., Muraki, E., Mislevy, R. J., & Bock, R. D. (2003). BILOG-MG
#'   3: Multiple-group IRT analysis and test maintenance for binary items
#'   (Computer software). Chicago, IL: Scientific Software International. URL
#'   http://www.ssicentral.com
#'
#' @examples
#' ## Example 1
#' # Import the "-prm.txt" output file from flexMIRT
#' flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt",
#'                         package = "irtQ")
#'
#' # Read item parameters and convert them to item metadata
#' bring.flexmirt(file = flex_sam, "par")$Group1$full_df
#'
#' ## Example 2
#' # Import the ".PAR" output file from PARSCALE
#' pscale_sam <- system.file("extdata", "parscale_sample.PAR", package = "irtQ")
#'
#' # Read item parameters and convert them to item metadata
#' bring.parscale(file = pscale_sam, "par")$full_df
#'
#' ## Example 3
#' # Import the item parameters of a model fitted with the mirt package
#' \donttest{
#' if (requireNamespace("mirt", quietly = TRUE)) {
#'   fit <- mirt::mirt(as.data.frame(LSAT6), 1, itemtype = "2PL",
#'                     verbose = FALSE)
#'   bring.mirt(fit)$full_df
#' }
#' }
#'
#' @export
bring.flexmirt <- function(file,
                           type = c("par", "sco"),
                           rePar = TRUE,
                           rePar.gpc = TRUE,
                           n.factor = 1) {


  # match the file type, which is "par" when omitted (the case is ignored)
  type <- match.arg(tolower(type), c("par", "sco"))

  type <- toupper(type)

  # When "par" file is inserted
  if (type == "PAR") {
    results <- bring.flexPrm(file = file, rePar = rePar, rePar.gpc = rePar.gpc)
  }

  # When "sco" file is inserted
  if (type == "SCO") {
    results <- bring.flexSco(file = file, n.factor = n.factor)
  }
  results
}

# "bring.flexPrm" is a modified version of "read.flexmirt" in R package of "rpf" (Pritikin, Weeks, Cai, Houts, and Chalmers, 2016)
# Read a flexMIRT PRM file
# Load the item parameters from a flexMIRT PRM file.
bring.flexPrm <- function(file, rePar = TRUE, rePar.gpc = TRUE) {
  groups <- list()

  ncol <- max(utils::count.fields(file, sep = "\t"))
  prm <- utils::read.delim(file, header = FALSE, col.names = 1:ncol)

  # determine the number of groups
  ng <- length(unique(prm[, 3]))

  for (g in 1:ng) {
    # select all items and this group
    index <- (prm[, 3] == g)
    thisGroup <- prm[index, ]

    g.dist <- list()
    g.label <- list()
    g.category <- list()
    g.dim <- list()
    g.param <- list()
    g.model <- list()

    for (i in 1:nrow(thisGroup)) {
      param <- c()
      dims <- thisGroup[i, 4]

      if (thisGroup[i, 1] == 1) { # item
        catg <- thisGroup[i, 6]

        if (thisGroup[i, 5] == 1) { # 3PL
          # grab item parameters
          logitg <- thisGroup[i, 7]
          c <- thisGroup[i, 8]
          a <- thisGroup[i, 9:(9 + dims - 1)]
          if (rePar) {
            if (!dims == 1) stop("Re-parameterization works only for one-dimensional model", call. = FALSE)
            param <- c(unlist(a), -c / unlist(a), exp(logitg) / (1 + exp(logitg)))
          } else {
            param <- c(unlist(a), c, logitg)
          }
          model <- "3PLM"
        }

        if (thisGroup[i, 5] == 2) { # graded
          # grab item parameters
          c <- rep(0, catg - 1)
          for (k in 1:(catg - 1)) {
            c[k] <- thisGroup[i, 6 + k]
          }
          a <- thisGroup[i, (6 + catg):(6 + catg + dims - 1)]
          if (rePar) {
            if (!dims == 1) stop("Re-parameterization works only for one-dimensional model", call. = FALSE)
            param <- c(unlist(a), -c / unlist(a))
          } else {
            param <- c(unlist(a), c)
          }
          model <- ifelse(catg == 2, "2PLM", "GRM")
        }

        if (thisGroup[i, 5] == 3) { # gpc
          offset <- 7
          # grab item parameters
          Tmat <- thisGroup[i, offset]
          offset <- offset + 1

          # scoring fn
          alf <- rep(0, catg - 1)
          for (k in 1:(catg - 1)) {
            alf[k] <- thisGroup[i, offset]
            offset <- offset + 1
          }
          alfGpc <- c(1, rep(0, catg - 2))
          isGpc <- all(alf == alfGpc)

          # slope
          a <- thisGroup[i, offset:(offset + dims - 1)]
          offset <- offset + dims

          # intercept (skip the type of the intercept matrix)
          offset <- offset + 1
          gam <- rep(0, catg - 1)
          for (k in 1:(catg - 1)) {
            gam[k] <- thisGroup[i, offset]
            offset <- offset + 1
          }

          # For GPC parameterization (only when dimension is one)
          if (dims == 1 & isGpc) {
            Tmat <- matrix(0, nrow = catg, ncol = (catg - 1))
            Tmat[, 1] <- 0:(catg - 1)
            if (catg > 2) {
              for (k in 2:(catg - 1)) {
                for (j in 2:(catg - 1)) {
                  Tmat[k, j] <- sin(pi * (j - 1) * (k - 1) / (catg - 1))
                }
              }
            }
            c.Vec <- Tmat %*% gam
            b <- -gam[1] / a
            d <- rep(0, length(c.Vec))
            for (k in 2:length(d)) {
              d[k] <- ((c.Vec[k] - c.Vec[k - 1]) / a) + b
            }
            d.new <- b - d[-1]
          }

          if (rePar.gpc) {
            if (!dims == 1 | !isGpc) stop("Re-parameterization works only for one-dimensional GPC model", call. = FALSE)
            if (dims == 1 & isGpc) param <- c(a, d.new)
          } else {
            param <- c(unlist(a), alf, gam)
          }
          if (catg == 2) model <- "2PLM"
          if (catg > 2 & isGpc) model <- "GPCM"
          if (catg > 2 & !isGpc) model <- "NOMM"
        }

        g.label[[i]] <- as.character(thisGroup[i, 2])
        g.category[[i]] <- catg
        g.dim[[i]] <- dims
        g.param[[i]] <- param
        g.model[[i]] <- model
      } else { # group
        if (thisGroup[i, 5] == 0) { # EmpHis = No (default)
          Cov <- matrix(NA, nrow = dims, ncol = dims)
          col <- 6 + dims
          for (rx in 1:dims) {
            for (cx in 1:rx) {
              Cov[rx, cx] <- thisGroup[i, col]
              if (rx != cx) Cov[cx, rx] <- thisGroup[i, col]
              col <- col + 1
            }
          }
          Mean <- unlist(thisGroup[i, 6:(6 + dims - 1)])
          names(Mean) <- paste0("Mu.", 1:length(Mean))
          colnames(Cov) <- paste0("Theta.", 1:ncol(Cov))
          rownames(Cov) <- paste0("Theta.", 1:nrow(Cov))
          g.dist <- list(Mean = Mean, Cov = Cov)
        }
        if (thisGroup[i, 5] == 1) { # EmpHis = Yes
          nquad <- thisGroup[i, 6]
          quads <- seq(from = -thisGroup[i, 7], to = thisGroup[i, 7], length.out = nquad)
          Cov <- matrix(NA, nrow = dims, ncol = dims)
          col <- 8 + dims
          for (rx in 1:dims) {
            for (cx in 1:rx) {
              Cov[rx, cx] <- thisGroup[i, col]
              if (rx != cx) Cov[cx, rx] <- thisGroup[i, col]
              col <- col + 1
            }
          }
          Mean <- unlist(thisGroup[i, 8:(8 + dims - 1)])
          names(Mean) <- paste0("Mu.", 1:length(Mean))
          colnames(Cov) <- paste0("Theta.", 1:ncol(Cov))
          rownames(Cov) <- paste0("Theta.", 1:nrow(Cov))
          weight <- as.numeric(thisGroup[i, col:(col + nquad - 1)])
          Emphist <- data.frame(theta = quads, weight = weight)
          g.dist <- list(Mean = Mean, Cov = Cov, Emphist = Emphist)
        }
      }
    } # for every item

    pmat <- matrix(NA, nrow = length(g.param), ncol = max(vapply(g.param, length, 0)))
    for (i in 1:length(g.param)) {
      v <- g.param[[i]]
      pmat[i, 1:length(v)] <- v
    }
    # DF <- data.frame(ID=unlist(g.label), FACTOR=unlist(g.dim), CATEGORY=unlist(g.category),
    # 				MODEL=unlist(g.model), PARAM=pmat)
    DF <- data.frame(
      id = unlist(g.label), cats = unlist(g.category),
      model = unlist(g.model), par = pmat, stringsAsFactors = FALSE
    )
    rownames(pmat) <- g.label
    groups[[g]] <- list(
      id = unlist(g.label), factor = unlist(g.dim), cats = unlist(g.category),
      model = unlist(g.model), par = pmat, full_df = DF, mean = g.dist$Mean,
      cov = g.dist$Cov, emphist = g.dist$Emphist
    )
    names(groups)[[g]] <- paste0("Group", g)
  } # for every group

  groups
}


# "bring.flexSco" function
# This function reads a scoring output of flexMIRT, "-sco.txt",
# and returns a list
#' @importFrom utils count.fields read.delim read.fwf read.table
bring.flexSco <- function(file, n.factor = 1) {

  groups <- list()

  ncol <- max(utils::count.fields(file, sep = ""))

  # read "-sco.txt" file
  sco <- utils::read.delim(file, header = FALSE, sep = "", col.names = 1:ncol)

  # (1) unidimensional model
  if (n.factor == 1) {
    if (ncol == 4) {
      if (typeof(sco[1, 3]) == "integer") {
        # determine the number of groups
        ng <- length(unique(sco[, 2]))
        for (g in 1:ng) {
          index <- (sco[, 2] == g)

          g.imputation <- sco[index, 1]
          g.ability <- sco[index, 4]
          g.se <- NULL

          groups[[g]] <- list(theta = g.ability, se = g.se, imputation = g.imputation)
          names(groups)[[g]] <- paste0("Group", g)
        }
      } else {
        # determine the number of groups
        ng <- length(unique(sco[, 1]))
        for (g in 1:ng) {
          index <- (sco[, 1] == g)

          g.ability <- sco[index, 3]
          g.se <- sco[index, 4]

          groups[[g]] <- list(theta = g.ability, se = g.se)
          names(groups)[[g]] <- paste0("Group", g)
        }
      }
    }

    if (ncol == 5) {
      # determine the number of groups
      ng <- length(unique(sco[, 1]))
      for (g in 1:ng) {
        index <- (sco[, 1] == g)

        g.iter <- sco[index, 3]
        g.ability <- sco[index, 4]
        g.se <- sco[index, 5]

        groups[[g]] <- list(theta = g.ability, se = g.se, iter = g.iter)
        names(groups)[[g]] <- paste0("Group", g)
      }
    }
  }

  # (2) multidimensional model
  if (n.factor > 1) {
    if (typeof(sco[1, 3]) == "integer") {
      # determine the number of groups
      ng <- length(unique(sco[, 2]))
      for (g in 1:ng) {
        index <- (sco[, 2] == g)

        g.imputation <- sco[index, 1]
        g.ability <- sco[index, 4:(4 + (n.factor - 1))]
        g.se <- NULL

        groups[[g]] <- list(theta = g.ability, se = g.se, imputation = g.imputation)
        names(groups)[[g]] <- paste0("Group", g)
      }
    } else {
      # determine the number of groups
      ng <- length(unique(sco[, 1]))
      for (g in 1:ng) {
        index <- (sco[, 1] == g)

        g.ability <- sco[index, 3:(3 + (n.factor - 1))]
        g.se <- NULL # will be updated

        groups[[g]] <- list(theta = g.ability, se = g.se)
        names(groups)[[g]] <- paste0("Group", g)
      }
    }
  }

  groups
}

######################################################################################################################################
# "bring.bilog" function
# This function is to read "par" or "sco" output files from Bilog-MG

#' @rdname bring.flexmirt
#' @importFrom utils count.fields read.delim read.fwf read.table
#' @export
bring.bilog <- function(file, type = c("par", "sco")) {

  # match the file type, which is "par" when omitted (the case is ignored)
  type <- match.arg(tolower(type), c("par", "sco"))

  type <- toupper(type)

  # When "par" file is inserted
  if (type == "PAR") {
    df <- utils::read.fwf(file, widths = c(8, 8, rep(10, 10)), skip = 4)
    # item IDs without the padding spaces
    df[, 1] <- trimws(as.character(df[, 1]))
    prm <- df[c(1, 5, 7, 11)]
    prm.df <- data.frame(
      id = prm[, 1], cats = rep(2, nrow(prm)), model = "DRM",
      par.1 = prm[, 2], par.2 = prm[, 3], par.3 = prm[, 4], stringsAsFactors = FALSE
    )
    results <- list(par = prm[-1], full_df = prm.df)
  }

  # When "sco" file is inserted
  if (type == "SCO") {
    # df <- utils::read.table(file, skip=2, fill=TRUE)
    df1 <- utils::read.fwf(file, widths = c(7, 8, 5, 5, 10, 12, 12), skip = 2)
    df2 <- utils::read.fwf(file, widths = c(3, 2, 30), skip = 2)
    num1 <- seq(2, nrow(df1), 2)
    num2 <- seq(1, (nrow(df2) - 1), 2)

    # Delete rows including examinees' id
    sco.df1 <- df1[num1, 3:7]
    sco.df2 <- df2[num2, c(1, 3)]
    sco.df <- cbind(sco.df2, sco.df1)
    colnames(sco.df) <- c("group", "id", "nitem", "ncorrect", "p", "theta", "se")

    # Splite dataframe based on a group factor
    sco.df$group <- as.factor(sco.df$group)

    # Make a list including all info per a group
    scoList <- list()
    ng <- length(levels(sco.df$group))
    dfList <- split(sco.df, sco.df$group)

    for (g in 1:ng) {
      g.name <- levels(sco.df$group)[g]
      g.ability <- dfList[[g]]$theta
      g.se <- dfList[[g]]$se
      g.df <- dfList[[g]]

      scoList[[g]] <- list(theta = g.ability, se = g.se, full_df = g.df)
      names(scoList)[[g]] <- as.character(g.name)
    }

    results <- scoList
  }
  results
}



######################################################################################################################################
# "bring.parscale" function
# This function is to read "par" or "sco" output files from Parscale
# "bring.parscale" is a modified version of "read.parscale" in R package of "plink" (Jonathan P. Weeks, 2017)

#' @rdname bring.flexmirt
#' @importFrom utils count.fields read.delim read.fwf read.table
#' @export
bring.parscale <- function(file, type = c("par", "sco")) {

  # match the file type, which is "par" when omitted (the case is ignored)
  type <- match.arg(tolower(type), c("par", "sco"))

  type <- toupper(type)

  ##   Import item parameters
  if (type == "PAR") {
    ##   Read the third line of the .PAR file which identifies the number
    ##   of parameter blocks, the total number of items, and the model code
    nums <- as.numeric(utils::read.fwf(file, c(8, rep(5, 5)), skip = 2, n = 1)[-1])

    ##   Read the fourth line in the .PAR file which identifies
    ##   the number of items in each block
    line_quot <- nums[1] %/% 30 # quotient
    numLines <- 2 * line_quot - 1
    if (nums[1] %% 30 != 0) {
      numLines <- numLines + 2
    }
    block_str1 <- readLines(file, n = (3 + numLines))[-c(1, 2, 3)] # strings containing block information
    block_str2 <- strsplit(block_str1, "\\s+")
    block_str3 <- as.numeric(do.call("c", block_str2))
    block <- block_str3[!is.na(block_str3)]

    ##   Initialize an object to store the item parameters
    pars <- NULL

    ##   Initialize a list to store the deviation threshold/step parameters for each block
    pars1 <- vector("list", length(block))

    ##   Initialize an object to store the number of response
    ##   categories for each item
    cat <- NULL

    ##   Identify the starting row for reading in the
    ##   parameters for a given block
    skip <- 4 + numLines

    ##   Identify the starting row for reading in the
    ##   category parameters for a given block
    skip1 <- skip + block[1]

    ##   Loop through each block of parameters
    for (i in 1:length(block)) {
      if (i > 1) {
        ##   Update the starting rows for the given block
        skip <- (4 + numLines) + sum(block[1:(i - 1)]) + (i - 1) * 2
        skip1 <- (4 + numLines) + sum(block[1:i]) + (i - 1) * 2
      }

      ##   Read the item parameters for the given block
      tmp <- utils::read.fwf(file, c(8, 5, 4, rep(10, 6)), skip = skip, n = block[i], colClasses = c("character", "numeric", "character", rep("numeric", 6)))

      ##   Add this block of parameters to the full set of parameters
      pars <- rbind(pars, tmp)

      ##   If there is more than one item in the block
      if (nrow(tmp) > 1) {
        cat <- c(cat, as.numeric(tmp[, 2]))

        ##   If there is only one item in the block
      } else {
        cat <- c(cat, as.numeric(tmp[2]))
      }

      ##   Read in the deviation threshold/step parameters parameters for the given block
      tmp1 <- utils::read.fwf(file, rep(10, cat[length(cat)]), skip = skip1, n = 2)

      ##   For the graded response models
      ##   Eliminate the last category parameter
      if (nums[3] < 5) {
        tmp1 <- tmp1[, -ncol(tmp1)]

        ##   For the  partial credit models
        ##   Eliminate the first category parameter
      } else {
        tmp1 <- tmp1[, -1]
      }
      pars1[[i]] <- tmp1
    }
    colnames(pars) <- c("block.name", "cat", "item", "slope", "slope.se", "location", "loc.se", "asymptote", "asymp.se")

    ##   For polytomous items
    if (max(cat) > 2) {
      ##   Create a matrix for the threshold/step parameters
      step <- matrix(NA, nrow(pars), 2 * (max(cat) - 1))

      ##   Because multiple items can all have the same
      ##   deviation threshold/step parameters, expand {pars1}
      ##   so that each item has a corresponding list element
      ##   containing threshold/step parameters
      pars1 <- rep(pars1, block)

      ##   Loop through all of the list elements for {pars1}
      ##   and insert the deviation threshold/step parameters
      ##   into the matrix {step}
      for (i in 1:length(pars1)) {
        ##   We only need deviation threshold/step parameters
        ##   for polytomous items. In the list, these will be
        ##   formatted as a data.frame whereas the parameters
        ##   for dichotomous items will be formatted as a
        ##   vector (of zeros)
        if (is.data.frame(pars1[[i]])) {
          step[i, ] <- unlist(pars1[[i]])
        }
      }

      if (nums[3] < 5) prefix <- "thresh" else prefix <- "step"

      ##   Create column names for the matrix of deviation threshold/step parameters
      step.names <- character()
      for (j in 1:(max(cat) - 1)) {
        step.names <- c(step.names, paste(prefix, j, sep = ""), paste(prefix, j, ".se", sep = ""))
      }
      colnames(step) <- step.names

      ##   Combine the matrix of item parameters with the matrix
      ##   of deviation threshold/step parameters
      pars <- cbind(pars, step)
    }
    pars_loc <- pars <- data.frame(pars)

    ##   For polytomous items, reformulate the parameters to
    ##   exclude the location parameter (if applicable)
    if (max(cat) > 2) {
      # stop for the model codes whose location parameter is not handled
      if (nums[3] %in% c(2, 4, 6, 8)) {
        stop("Model codes 2, 4, 6, and 8 of PARSCALE are not supported for polytomous items.", call. = FALSE)
      }
      if (nums[3] %in% c(1, 3, 5, 7)) {
        pars[cat > 2, seq(10, ncol(pars), 2)] <- pars[cat > 2, 6] - pars[cat > 2, seq(10, ncol(pars), 2)]
        pars[cat > 2, 6:7] <- 0
      }
    }

    ##   Get rid of all the columns of containing standard errors
    ##   the item names, and the number of categories
    pars_nose <- pars[, seq(4, ncol(pars), 2)]

    ##   reverse the order of threshold/step parameters for polytomous items
    # if(max(cat)>2) {
    #   pars_nose <- cbind(pars_nose[, 1:3], rev(pars_nose[, 4:ncol(pars_nose)]))
    # }

    ##   Make sure that the matrix of parameters is formatted properly
    ##   The threshold/step parameters in the matrix need to be shifted
    ##   one or two columns to the left
    if (max(cat) > 2) {
      pars_nose[cat > 2, 2:(ncol(pars_nose) - 2)] <- pars_nose[cat > 2, 4:ncol(pars_nose)]
      pars_nose <- pars_nose[, 1:max(cat)]
    }
    colnames(pars_nose) <- paste0("par.", 1:ncol(pars_nose))

    ##   Number of items
    n <- nrow(pars)

    ##   Model of items
    if (nums[3] < 5) {
      mod <- ifelse(cat == 2, "DRM", "GRM")
    } else {
      mod <- ifelse(cat == 2, "DRM", "GPCM")
    }

    ##   create a full data frame
    # pars_full <- data.frame(ID=pars$item, FACTOR=rep(1, n), CATEGORY=cat, MODEL=mod, pars_nose)
    pars_full <- data.frame(id = pars$item, cats = cat, model = mod, pars_nose, stringsAsFactors = FALSE)

    ##   a list to be returned
    rr <- list(
      id = pars$item, factor = rep(1, n), cats = cat, model = mod,
      par_parscale = pars_loc, par = pars_nose, full_df = pars_full
    )
  }

  ##   Import ability estimates
  if (type == "SCO") {
    pars <- read.fwf(file, c(21, 32, 12, 12))
    id <- as.character(pars[seq(1, nrow(pars), 2), 1])
    ability <- pars[seq(2, nrow(pars), 2), 3:4]
    ability[ability == 999] <- NA
    pars <- cbind(id, ability)
    colnames(pars) <- c("id", "theta", "se")
    rr <- pars
  }

  return(rr)
}


######################################################################################################################################
# "bring.mirt" function
# This function is to read parameter estimates of an object obtained from MIRT

#' @rdname bring.flexmirt
#' @importFrom utils count.fields read.delim read.fwf read.table
#' @export
bring.mirt <- function(x) {
  # bring.mirt() needs the mirt package, which is only suggested
  if (!requireNamespace("mirt", quietly = TRUE)) {
    stop("Package 'mirt' is required for bring.mirt(). Please install it.", call. = FALSE)
  }

  # stop when the object is a model of multiple groups
  if (inherits(x, "MultipleGroupClass")) {
    stop(
      "bring.mirt() does not support multiple-group mirt objects; ",
      "extract one group with mirt::extract.group() first.",
      call. = FALSE
    )
  }

  # read parameter estimates from an object of mirt
  prm_all <- mirt::coef(x, simplify = TRUE, IRTpars = TRUE)
  cats <- x@Data$K
  any_py <- any(cats > 2)

  # extract other information from an object
  wts <- as.numeric(x@Internals$Prior[[1]])
  nods <- as.numeric(x@Model$Theta)
  EMPHIST <- cbind(nods = nods, weights = wts)
  mods <- x@Model$itemtype
  MODEL <- mods
  MODEL[MODEL %in% c("gpcm", "gpcmIRT")] <- "GPCM"
  MODEL[MODEL == "Rasch" & cats > 2] <- "GPCM"
  MODEL[MODEL == "graded"] <- "GRM"
  MODEL[MODEL == "2PL" | MODEL == "3PL"] <- "DRM"
  MODEL[MODEL == "Rasch" & cats == 2L] <- "DRM"
  FACTOR <- mirt::extract.mirt(x, what = "nfact")
  ID <- mirt::extract.mirt(x, what = "itemnames")
  nitem <- mirt::extract.mirt(x, what = "nitems")
  MEAN <- as.numeric(prm_all$mean)
  COV <- as.numeric(prm_all$cov)

  # check dimensionality
  if (FACTOR > 1) stop("This function only works for an unidimensional model", call. = FALSE)

  if (!any_py) {
    # check models for dichotomous data
    if (!all(mods %in% c("Rasch", "2PL", "3PL"))) stop("For dichotomous data, only Rasch, 2PL, 3PL models are available", call. = FALSE)
  }

  if (any_py) {
    mods_py <- mods[cats > 2]

    # check models for polytomous data
    if (!all(mods_py %in% c("Rasch", "graded", "gpcm", "gpcmIRT"))) stop("For polytomous data, only Rasch, graded, gpcm, gpcmIRT models are available", call. = FALSE)

    # check models for the two-category items of mixed-format data
    if (!all(mods[cats == 2L] %in% c("Rasch", "2PL", "3PL", "graded", "gpcm", "gpcmIRT"))) {
      stop("For dichotomous items, only Rasch, 2PL, 3PL, graded, gpcm, gpcmIRT models are available", call. = FALSE)
    }
  }

  # the coefficient matrix, whose columns follow the order of the first item of each model
  prm_coef <- prm_all$items

  # names of the threshold columns sorted by their numbers
  b_cols <- grep("^b[0-9]+$", colnames(prm_coef), value = TRUE)
  b_cols <- b_cols[order(as.numeric(sub("^b", "", b_cols)))]

  # a matrix of item parameters, filled item by item with the column names
  prm_item <- matrix(NA_real_,
    nrow = nitem, ncol = ifelse(any_py, max(cats), 3L),
    dimnames = list(rownames(prm_coef), NULL)
  )
  for (i in seq_len(nitem)) {
    if (mods[i] %in% c("Rasch", "2PL", "3PL") && cats[i] == 2L) {
      # slope, difficulty, and guessing parameters of a dichotomous item
      prm_item[i, 1:3] <- prm_coef[i, c("a", "b", "g")]
    } else {
      if (mods[i] == "gpcmIRT") {
        # the intercept c of the gpcmIRT model is moved into the thresholds (b_v - c / a)
        slope <- prm_coef[i, "a1"]
        thr <- prm_coef[i, b_cols[seq_len(cats[i] - 1)]] - prm_coef[i, "c"] / slope
      } else {
        # slope and thresholds of a graded or partial credit item
        slope <- prm_coef[i, "a"]
        thr <- prm_coef[i, b_cols[seq_len(cats[i] - 1)]]
      }
      if (cats[i] == 2L) {
        # a two-category item fitted with a polytomous model is a dichotomous item without guessing
        prm_item[i, 1:3] <- c(slope, thr, 0)
      } else {
        # slope and thresholds of a polytomous item
        prm_item[i, seq_len(cats[i])] <- c(slope, thr)
      }
    }
  }

  # a two-category item is a dichotomous item in irtQ
  MODEL[cats == 2L] <- "DRM"

  # give colnames for a prm matrix
  colnames(prm_item) <- paste0("par.", 1:ncol(prm_item))

  # original item parameters matrix obtained from mirt
  prm_mirt <- mirt::coef(x, simplify = TRUE, IRTpars = FALSE)$items

  # full data.frame of item parameters
  # prm_full <- data.frame(ID=ID, FACTOR=rep(FACTOR, nitem), CATEGORY=cats, MODEL=MODEL, prm_item)
  prm_full <- data.frame(id = ID, cats = cats, model = MODEL, prm_item, stringsAsFactors = FALSE)

  rr <- list(
    id = ID, factor = rep(FACTOR, nitem), cats = cats, model = list(old = mods, new = MODEL),
    emphist = EMPHIST, par_mirt = prm_mirt, par = prm_item, full_df = prm_full,
    ability = c(MEAN = MEAN, COV = COV)
  )

  rr
}
