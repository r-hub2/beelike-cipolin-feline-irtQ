#' Create a Data Frame of Item Metadata
#'
#' This function creates a data frame of item metadata (item parameters,
#' numbers of score categories, and IRT models) for use in IRT analyses with
#' the \pkg{irtQ} package.
#'
#' @param par.drm A list containing three numeric vectors for dichotomous item
#'   parameters: item discrimination (`a`), item difficulty (`b`), and guessing
#'   parameters (`g`), each with one value per item with `cats = 2`. A single
#'   value of `g` is used for every dichotomous item. If `g` is omitted or
#'   `NULL`, the guessing parameters are set to 0.
#' @param par.prm A list containing polytomous item parameters. The list must
#'   include a numeric vector `a` of slope parameters and a list `d` of numeric
#'   vectors of threshold parameters, with one slope and one vector of
#'   `cats - 1` thresholds per item with `cats > 2`. See the **Details** section
#'   for more information.
#' @param item.id A character vector of item IDs. If `NULL`, default IDs (e.g.,
#'   "V1", "V2", ...) are assigned automatically.
#' @param cats A numeric vector indicating the number of score categories for
#'   each item.
#' @param model A character vector specifying the IRT model for each item.
#'   Available options are `"1PLM"`, `"2PLM"`, `"3PLM"`, and `"DRM"` for
#'   dichotomous items, and `"GRM"` and `"GPCM"` for polytomous items. The label
#'   `"DRM"` serves as a general category that encompasses all dichotomous
#'   models (`"1PLM"`, `"2PLM"`, and `"3PLM"`), while `"GRM"` and `"GPCM"` refer
#'   to the graded response model and (generalized) partial credit model,
#'   respectively.
#' @param default.par Logical. If `TRUE`, default item parameters are generated
#'   based on the specified `cats` and `model`. In this case, the slope
#'   parameter is set to 1, all difficulty (or threshold) parameters are set to
#'   0, and the guessing parameter is set to 0.2 for `"3PLM"` or `"DRM"` items.
#'   The default is `FALSE`.
#'
#' @details For any item where `"1PLM"` or `"2PLM"` is specified in `model`, the
#'   guessing parameter is set to `NA`. If `cats` or `model` has length 1, it is
#'   recycled across all items.
#'
#'   As in the [irtQ::simdat()] function, when constructing a mixed-format test
#'   form, it is important to specify the `cats` argument to reflect the correct
#'   number of score categories for each item, in the exact order that the items
#'   appear. See [irtQ::simdat()] for further guidance on how to specify `cats`.
#'
#'   When specifying item parameters using `par.drm` and/or `par.prm`, the
#'   internal structure and ordering of elements must be followed.
#'   - `par.drm` should be a list with three components:
#'     - `a`: a numeric vector of slope parameters
#'     - `b`: a numeric vector of difficulty parameters
#'     - `g`: a numeric vector of guessing parameters
#'   - `par.prm` should be a list with two components:
#'     - `a`: a numeric vector of slope parameters for polytomous items
#'     - `d`: a list of numeric vectors specifying the threshold parameters
#'       for each polytomous item
#'
#'   For GPCM items, the threshold parameters are \eqn{b_v = \beta - \tau_v},
#'   the overall item location minus the threshold of each score category. An
#'   item with *K* score categories requires *K - 1* threshold parameters,
#'   because the term for the lowest score category is fixed at 0 and is not
#'   supplied. For GRM items, the thresholds must be in increasing order.
#'
#' @return A data frame of item metadata with columns `id`, `cats`, `model`,
#'   and `par.1`, `par.2`, ..., one row per item. The number of parameter
#'   columns is the larger of 3 and the maximum of `cats`; unused cells are
#'   `NA`. This data frame can be used as input for other functions in the
#'   \pkg{irtQ} package, such as [irtQ::est_irt()] or [irtQ::simdat()].
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [irtQ::est_irt()], [irtQ::simdat()], [irtQ::shape_df_fipc()]
#'
#' @examples
#' ## A mixed-format test form
#' ## containing five dichotomous items and two polytomous items
#' # Create a list of dichotomous item parameters
#' par.drm <- list(
#'   a = c(1.1, 1.2, 0.9, 1.8, 1.4),
#'   b = c(0.1, -1.6, -0.2, 1.0, 1.2),
#'   g = rep(0.2, 5)
#' )
#'
#' # Create a list of polytomous item parameters
#' par.prm <- list(
#'   a = c(1.4, 0.6),
#'   d = list(
#'     c(-1.9, 0.0, 1.2),
#'     c(0.4, -1.1, 1.5, 0.2)
#'   )
#' )
#'
#' # Create a numeric vector indicating the number of score categories for each item
#' cats <- c(2, 4, 2, 2, 5, 2, 2)
#'
#' # Create a character vector specifying the IRT model for each item
#' model <- c("DRM", "GRM", "DRM", "DRM", "GPCM", "DRM", "DRM")
#'
#' # Generate an item metadata set using the specified parameters
#' shape_df(par.drm = par.drm, par.prm = par.prm, cats = cats, model = model)
#'
#' ## An item metadata frame with default parameters for five dichotomous and two polytomous items
#' # Create a numeric vector indicating the number of score categories for each item
#' cats <- c(2, 4, 3, 2, 5, 2, 2)
#'
#' # Create a character vector specifying the IRT model for each item
#' model <- c("1PLM", "GRM", "GRM", "2PLM", "GPCM", "DRM", "3PLM")
#'
#' # Generate an item metadata frame with default parameters
#' shape_df(cats = cats, model = model, default.par = TRUE)
#'
#' ## A single-format test form consisting of five dichotomous items
#' # Generate the item metadata
#' shape_df(par.drm = par.drm, cats = rep(2, 5), model = "DRM")
#'
#' @export
#'
shape_df <- function(par.drm = list(a = NULL, b = NULL, g = NULL),
                     par.prm = list(a = NULL, d = NULL),
                     item.id = NULL,
                     cats,
                     model,
                     default.par = FALSE) {

  # convert the model names to uppercase
  model <- toupper(model)

  # check model names
  if (!all(model %in% c("1PLM", "2PLM", "3PLM", "DRM", "GRM", "GPCM"))) {
    stop(paste0(
      "At least one model name is mis-specified in the model argument.\n",
      "Available model names are 1PLM, 2PLM, 3PLM, DRM, GRM, and GPCM"
    ), call. = FALSE)
  }

  # generate default item parameters when default.par = TRUE
  if (default.par) {
    if (missing(cats) | missing(model)) {
      stop("The number of score categories and IRT models must be specified.", call. = FALSE)
    }

    # count the items from the longest of cats, model, and item.id
    n_item_default <- max(length(cats), length(model), length(item.id))

    # replicate a single cats value across all items
    if (length(cats) == 1) cats <- rep(cats, n_item_default)

    # replicate a single model name across all items
    if (length(model) == 1) model <- rep(model, n_item_default)

    # stop when cats, model, or item.id cannot be matched to the number of items
    if (length(cats) != n_item_default || length(model) != n_item_default ||
        (!is.null(item.id) && length(item.id) != n_item_default)) {
      stop(
        "`cats` and `model` must have length 1 or one value per item, ",
        "and `item.id` must have one value per item.",
        call. = FALSE
      )
    }

    # find the index of drm items
    idx.drm <- which(cats == 2)
    if (sum(idx.drm) == 0) idx.drm <- NULL

    # find the index of prm items
    idx.prm <- which(cats > 2)
    if (sum(idx.prm) == 0) idx.prm <- NULL

    # assign default values to the item parameters
    if (!is.null(idx.drm)) {
      par.drm <- list(
        a = rep(1, length(idx.drm)),
        b = rep(0, length(idx.drm)), g = rep(NA_real_, length(idx.drm))
      )
      par.drm$g[model[idx.drm] == "3PLM" | model[idx.drm] == "DRM"] <- 0.2
    } else {
      par.drm <- list(a = NULL, b = NULL, g = NULL)
    }
    if (!is.null(idx.prm)) {
      par.prm <- list(a = rep(1, length(idx.prm)), d = vector("list", length(idx.prm)))
      for (i in 1:length(idx.prm)) {
        par.prm$d[[i]] <- rep(0, cats[idx.prm[i]] - 1)
      }
    } else {
      par.prm <- list(a = NULL, d = NULL)
    }
  }

  # number of the items
  nitem <- length(par.drm$b) + length(par.prm$d)

  # max score categories
  max.cat <- max(cats)

  # create a vector of item ids when item.id = NULL
  if (is.null(item.id)) item.id <- paste0("V", 1:nitem)

  # create a vector of score categories when length(cats) = 1
  if (length(cats) == 1) cats <- rep(cats, nitem)

  # create a vector of model names when length(model) = 1
  if (length(model) == 1) model <- rep(model, nitem)

  # stop when cats, model, or item.id do not match the number of items
  if (length(cats) != nitem || length(model) != nitem || length(item.id) != nitem) {
    stop(
      "`cats` and `model` must have length 1 or one value per item, and ",
      "`item.id` one value per item (", nitem, " items given in `par.drm` and `par.prm`).",
      call. = FALSE
    )
  }

  # find the indices of DRM and PRM items when default.par = FALSE
  if (!default.par) {
    # find the index of drm items
    idx.drm <- which(cats == 2)
    if (sum(idx.drm) == 0) idx.drm <- NULL

    # find the index of prm items
    idx.prm <- which(cats > 2)
    if (sum(idx.prm) == 0) idx.prm <- NULL

    # repeat a single guessing value for every dichotomous item
    if (length(par.drm) >= 3L && length(par.drm[[3]]) == 1L && length(idx.drm) > 1L) {
      par.drm[[3]] <- rep(par.drm[[3]], length(idx.drm))
    }

    # count the guessing values given (zero when the g element is absent)
    n.g <- if (length(par.drm) < 3L) 0L else length(par.drm[[3]])

    # stop when the dichotomous parameter vectors do not match the dichotomous items
    if (length(par.drm[[1]]) != length(idx.drm) || length(par.drm[[2]]) != length(idx.drm) ||
        !(n.g %in% c(0L, length(idx.drm)))) {
      stop("`par.drm` must give a, b (and g) for every item with cats = 2.", call. = FALSE)
    }

    # stop when the polytomous slopes or threshold lists do not match the polytomous items
    if (length(par.prm[[1]]) != length(idx.prm) || length(par.prm[[2]]) != length(idx.prm) ||
        any(lengths(par.prm[[2]]) != cats[idx.prm] - 1)) {
      stop("`par.prm` must give a and cats - 1 thresholds for every item with cats > 2.", call. = FALSE)
    }
  }

  # create an empty matrix to contain item parameters
  if (is.null(idx.prm)) {
    par_mat <- array(NA, c(nitem, 3))
  } else {
    par_mat <- array(NA, c(nitem, max.cat))
  }

  # if drm items exist
  if (!is.null(idx.drm)) {
    # fill zero guessing values when g is absent or NULL
    if (length(par.drm) < 3L || is.null(par.drm[[3]])) par.drm[[3]] <- rep(0, length(idx.drm))
    par_mat[idx.drm, 1:3] <- bind.fill(par.drm, type = "cbind")
  }

  # if prm items exist
  if (!is.null(idx.prm)) {
    par_mat[idx.prm, 1] <- par.prm[[1]]
    par_mat[idx.prm, 2:max.cat] <- bind.fill(par.prm[[2]], type = "rbind")
  }

  # create an item metadata
  x <- data.frame(id = item.id, cats = cats, model = model, par_mat, stringsAsFactors = FALSE)

  # re-assign column names
  colnames(x) <- c("id", "cats", "model", paste0("par.", 1:(ncol(x) - 3)))

  # assign NAs to the par.3 column for the 1PLM and 2PLM items
  x[x$model %in% c("1PLM", "2PLM"), "par.3"] <- NA_real_

  # assign 0s to the par.3 column for the 3PLM when par.3 = NA
  x[x$model == "3PLM" & is.na(x$par.3), "par.3"] <- 0

  # last check
  if (any(x[x$cats == 2, 3] %in% c("GRM", "GPCM"))) {
    stop("Dichotomous items must have models among '1PLM', '2PLM', '3PLM', and 'DRM'.", call. = FALSE)
  }
  if (any(x[x$cats > 2, 3] %in% c("1PLM", "2PLM", "3PLM", "DRM"))) {
    stop("Polytomous items must have models among 'GRM' and 'GPCM'.", call. = FALSE)
  }

  # return the results
  x
}


# Create item metadata containing the starting values
startval_df <- function(cats, model, item.id = NULL) {
  # convert the model names to uppercase
  model <- toupper(model)

  # check model names
  if (!all(model %in% c("1PLM", "2PLM", "3PLM", "DRM", "GRM", "GPCM"))) {
    stop(paste0(
      "At least one model name is mis-specified in the model argument.\n",
      "Available model names are 1PLM, 2PLM, 3PLM, DRM, GRM, and GPCM"
    ), call. = FALSE)
  }

  # count the items from the longest of cats, model, and item.id
  n_item_start <- max(length(cats), length(model), length(item.id))

  # replicate a single cats value across all items
  if (length(cats) == 1) cats <- rep(cats, n_item_start)

  # replicate a single model name across all items
  if (length(model) == 1) model <- rep(model, n_item_start)

  # find the index of drm items
  idx.drm <- which(cats == 2)
  if (sum(idx.drm) == 0) idx.drm <- NULL

  # find the index of prm items
  idx.prm <- which(cats > 2)
  if (sum(idx.prm) == 0) idx.prm <- NULL

  # assign default values to the item parameters
  if (!is.null(idx.drm)) {
    par.drm <- list(
      a = rep(1, length(idx.drm)),
      b = rep(0, length(idx.drm)), g = rep(0, length(idx.drm))
    )
    par.drm$g[model[idx.drm] == "3PLM" | model[idx.drm] == "DRM"] <- 0.2
  } else {
    par.drm <- list(a = NULL, b = NULL, g = NULL)
  }
  if (!is.null(idx.prm)) {
    par.prm <- list(a = rep(1, length(idx.prm)), d = vector("list", length(idx.prm)))
    for (i in 1:length(idx.prm)) {
      par.prm$d[[i]] <- seq(-1.0, 1.0, length.out = (cats[idx.prm[i]] - 1))
    }
  } else {
    par.prm <- list(a = NULL, d = NULL)
  }

  # number of the items
  nitem <- length(par.drm$b) + length(par.prm$d)

  # max score categories
  max.cat <- max(cats)

  # create a vector of item ids when item.id = NULL
  if (is.null(item.id)) item.id <- paste0("V", 1:nitem)

  # create a vector of score categories when length(cats) = 1
  if (length(cats) == 1) cats <- rep(cats, nitem)

  # create a vector of model names when length(model) = 1
  if (length(model) == 1) model <- rep(model, nitem)

  # create an empty matrix to contain item parameters
  if (is.null(idx.prm)) {
    par_mat <- array(NA, c(nitem, 3))
  } else {
    par_mat <- array(NA, c(nitem, max.cat))
  }

  # if drm items exist
  if (!is.null(idx.drm)) {
    if (is.null(par.drm[[3]])) par.drm[[3]] <- rep(0, length(idx.drm))
    par_mat[idx.drm, 1:3] <- bind.fill(par.drm, type = "cbind")
  }

  # if prm items exist
  if (!is.null(idx.prm)) {
    par_mat[idx.prm, 1] <- par.prm[[1]]
    par_mat[idx.prm, 2:max.cat] <- bind.fill(par.prm[[2]], type = "rbind")
  }

  # create an item metadata
  x <- data.frame(id = item.id, cats = cats, model = model, par_mat, stringsAsFactors = FALSE)

  # re-assign column names
  colnames(x) <- c("id", "cats", "model", paste0("par.", 1:(ncol(x) - 3)))

  # assign 0s to the par.3 column for the 3PLM when par.3 = NA
  x[x$model == "3PLM" & is.na(x$par.3), "par.3"] <- 0

  # last check
  if (any(x[x$cats == 2, 3] %in% c("GRM", "GPCM"))) {
    stop("Dichotomous items must have models among '1PLM', '2PLM', '3PLM', and 'DRM'.", call. = FALSE)
  }
  if (any(x[x$cats > 2, 3] %in% c("1PLM", "2PLM", "3PLM", "DRM"))) {
    stop("Polytomous items must have models among 'GRM' and 'GPCM'.", call. = FALSE)
  }

  # return the results
  x
}
