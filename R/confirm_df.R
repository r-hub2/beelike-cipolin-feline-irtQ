# Check that the item metadata is a correctly formatted data frame with all
# components required for IRT analyses, and standardize it
confirm_df <- function(x, g2na = FALSE) {
  # change all factor variables into character variables
  x <- data.frame(x, stringsAsFactors = FALSE)
  x <- purrr::modify_if(x, is.factor, as.character)
  x[, 3] <- toupper(x[, 3])

  # convert a cats column stored as text or factor labels into numbers
  if (!is.numeric(x[, 2])) {
    cats.num <- suppressWarnings(as.numeric(x[, 2]))
    # stop when a non-missing value cannot be read as a number
    if (anyNA(cats.num[!is.na(x[, 2])])) {
      stop("The cats column in 'x' must be numeric.", call. = FALSE)
    }
    x[, 2] <- cats.num
  }
  modelGood <- all(x[, 3] %in% c("1PLM", "2PLM", "3PLM", "DRM", "GRM", "GPCM"))
  # require a known number of at least two score categories for every item
  catsGood <- !anyNA(x[, 2]) && all(x[, 2] >= 2)
  if (!modelGood) {
    stop(paste0(
      "At least one model name is mis-specified in the model column.\n",
      "Available model names are 1PLM, 2PLM, 3PLM, DRM, GRM, and GPCM"
    ), call. = FALSE)
  }
  if (!catsGood) {
    stop(paste0(
      "At least one score category count in the cats column is missing or less than 2.\n",
      "Each item must have at least two score categories."
    ), call. = FALSE)
  }

  # add a par.3 (guessing parameter) column when it is missing and all items are 1PLM or 2PLM
  if (ncol(x[, -c(1, 2, 3)]) == 2) {
    if (all(x[, 3] %in% c("1PLM", "2PLM"))) {
      x <- data.frame(x, par.3 = NA)
    } else {
      stop("Add a par.3 column to the item metadata 'x' unless all items are 1PLM or 2PLM.", call. = FALSE)
    }
  }

  # remove the parameter columns with all NAs except par.1
  col.na.lg <- Rfast::colAll(is.na(x[, -c(1:4)]))
  if (!all(col.na.lg)) {
    x <- x[, c(!logical(4), !col.na.lg)]
  }

  # re-assign column names
  colnames(x) <- c("id", "cats", "model", paste0("par.", 1:(ncol(x) - 3)))

  # check that every parameter column is numeric or logical (all NA)
  par.num <- vapply(
    x[, -(1:3), drop = FALSE],
    function(v) is.numeric(v) || is.logical(v), logical(1)
  )
  # stop on text columns, which breakdown() would turn into factor codes
  if (!all(par.num)) {
    stop("All item parameter columns in 'x' must be numeric.", call. = FALSE)
  }

  # stop when a dichotomous model is given a number of score categories other than two
  if (any(x$model %in% c("1PLM", "2PLM", "3PLM", "DRM") & x$cats != 2)) {
    stop("Dichotomous items (1PLM, 2PLM, 3PLM, DRM) must have cats = 2.", call. = FALSE)
  }

  # threshold columns that follow the slope column
  thr <- as.matrix(x[, -(1:4), drop = FALSE])
  if (ncol(thr) > 0L) {
    # flag thresholds placed beyond the cats - 1 thresholds of a polytomous item
    thr.extra <- (x$model %in% c("GRM", "GPCM")) & (col(thr) > (x$cats - 1)) & !is.na(thr)
    # stop instead of using the extra thresholds in the category probabilities
    if (any(thr.extra)) {
      stop("A polytomous item has more threshold parameters than cats - 1.", call. = FALSE)
    }
  }

  # handle the g parameters for the 1PLM and 2PLM
  if (g2na) {
    # assign NAs to the par.3 column for the 1PLM, 2PLM items
    x[x$model %in% c("1PLM", "2PLM"), "par.3"] <- NA_real_
  } else {
    # assign 0 values to the par.3 column for the 1PLM, 2PLM items
    x[x$model %in% c("1PLM", "2PLM"), "par.3"] <- 0
  }

  # consider DRM as 3PLM
  if ("DRM" %in% x$model) {
    x$model[x$model == "DRM"] <- "3PLM"
    memo <- "All 'DRM' items are treated as '3PLM' items.\n"
    warning(memo, call. = FALSE)
  }

  # return the results
  x
}
