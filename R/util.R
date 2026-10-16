#' Bind Fill
#'
#' This function creates a matrix using either row-wise (`rbind`) or column-wise
#' (`cbind`) binding of a list of numeric vectors with varying lengths. Shorter
#' vectors are padded with a specified fill value.
#'
#' @param List A list containing numeric vectors of possibly different lengths.
#' @param type A character string indicating the type of binding to perform.
#'   Options are `"rbind"` or `"cbind"`. Default is `"rbind"`.
#' @param fill A value used to fill missing elements when aligning the vectors.
#'   For `type = "cbind"`, this fills missing rows in shorter columns; for `type
#'   = "rbind"`, this fills missing columns in shorter rows. Accepts any R
#'   object (e.g., numeric, character, logical). Default is `NA`.
#'
#' @return A matrix formed by binding the elements of the list either row-wise
#'   or column-wise, with shorter vectors padded by the specified `fill` value.
#'   The matrix has no row or column names.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @examples
#' # Sample list
#' score_list <- list(item1 = 0:3, item2 = 0:2, item3 = 0:5, item4 = 0:4)
#'
#' # 1) Create a row-bound matrix (rbind)
#' bind.fill(score_list, type = "rbind")
#'
#' # 2) Create a column-bound matrix (cbind)
#' bind.fill(score_list, type = "cbind")
#'
#' # 3) Create a column-bound matrix and fill missing values with 0
#' bind.fill(score_list, type = "cbind", fill = 0L)
#'
#' @export
#' @import dplyr
bind.fill <- function(List, type=c("rbind", "cbind"), fill = NA){
  type <- tolower(type)
  type <- match.arg(type)
  nm <- List
  nm <- purrr::map(nm, as.matrix)
  names(nm) <- 1:length(nm)
  n <- max(purrr::map_dbl(nm, nrow))
  df <-
    purrr::map_dfc(nm, function(x) {rbind(x, matrix(fill, n - nrow(x), ncol(x)))}) %>%
    as.matrix()
  switch(type,
         cbind = unname(df),
         rbind = unname(t(df))
  )

}

# Check the positions of the items to be skipped in an IPD analysis and return them as
# sorted unique integers (NULL when no item is skipped)
check_item_skip <- function(item.skip, nitem) {
  # nothing to check when no item is skipped
  if (is.null(item.skip)) {
    return(NULL)
  }

  # read a logical vector as the positions of TRUE
  if (is.logical(item.skip)) {
    if (length(item.skip) != nitem || anyNA(item.skip)) {
      stop("A logical 'item.skip' must have one non-missing value per item.", call. = FALSE)
    }
    item.skip <- which(item.skip)
  }

  # stop when the values are not whole item positions within the items
  if (!is.numeric(item.skip) || anyNA(item.skip) || any(item.skip != round(item.skip)) ||
      any(item.skip < 1 | item.skip > nitem)) {
    stop("'item.skip' must contain item positions (row numbers of 'x') from 1 to ", nitem, ".",
         call. = FALSE)
  }

  # return the unique positions in ascending order
  item.skip <- sort(unique(as.integer(item.skip)))
  if (length(item.skip) == 0L) item.skip <- NULL
  item.skip
}

# Convert the response data of the model-data fit functions to a numeric matrix
# and check that every observed response is a whole number within the score
# categories of its item
resp_to_matrix <- function(data, cats, id = seq_along(cats)) {
  # convert the responses to a matrix
  resp <- as.matrix(data)

  # stop when the data do not have one column per item
  if (ncol(resp) != length(cats)) {
    stop(
      "The number of columns in 'data' (", ncol(resp), ") must equal ",
      "the number of items in 'x' (", length(cats), ").",
      call. = FALSE
    )
  }

  # read character or factor responses as numbers
  if (!is.numeric(resp) && !is.logical(resp)) {
    resp_num <- suppressWarnings(as.numeric(resp))

    # stop when an observed response cannot be read as a number
    if (any(is.na(resp_num) & !is.na(resp))) {
      stop("Responses must be numeric scores.", call. = FALSE)
    }

    # keep the numeric responses with the dimensions and names of the input
    resp <- matrix(resp_num, nrow = nrow(resp), dimnames = dimnames(resp))
  }

  # read logical responses as 0 and 1
  if (is.logical(resp)) {
    resp <- matrix(as.numeric(resp), nrow = nrow(resp), dimnames = dimnames(resp))
  }

  # flag the items that have an observed response outside their score categories
  max_cat <- matrix(cats - 1, nrow = nrow(resp), ncol = ncol(resp), byrow = TRUE)
  bad_item <- which(colSums(!is.na(resp) & (resp < 0 | resp > max_cat | resp != round(resp)),
    na.rm = TRUE) > 0)

  # stop when any observed response is invalid
  if (length(bad_item) > 0L) {
    stop(
      "Responses outside the score categories 0, ..., (cats - 1) are found for item(s) ",
      paste0(id[bad_item], " (column ", bad_item, ")", collapse = ", "), ".",
      call. = FALSE
    )
  }

  resp
}

# Find the indices of the DRM and PRM items in the output of breakdown()
idxfinder <- function(x) {

  # find the index of drm items
  idx.drm <- which(x$cats == 2)
  if(sum(idx.drm) == 0) idx.drm <- NULL

  # find the index of prm items
  idx.prm <- which(x$cats > 2)
  if(sum(idx.prm) == 0) idx.prm <- NULL

  # return the results
  list(idx.drm = idx.drm, idx.prm = idx.prm)

}

# Compute the mean and variance of a discrete distribution given its nodes and weights
cal_moment <- function(node, weight) {
  mu <- sum(node * weight)
  sigma2 <- sum(node^2 * weight) - mu^2
  rst <- c(mu=mu, sigma2=sigma2)
  rst
}


# Build the per-item frequency-category list used by divide_data() and
# info_xpd() in the EM pipeline.  For each item k the result is a
# `nstd x cats[k]` integer matrix where entry (i, j) = 1 iff examinee i
# selected response category (j - 1), and NA responses contribute
# all-zero rows (no category indicator is on for a missing response).
#
# Args:
#   data : numeric matrix of responses (nstd x nitem); responses are
#          coded 0..cats[k]-1, with NA for missing
#   cats : integer vector of length nitem giving the number of score
#          categories per item
#
# Returns:
#   list of length nitem; element [[k]] is an nstd x cats[k] integer
#   matrix as described above.
build_freqcat <- function(data, cats) {

  # number of examinees and items in the response matrix
  nstd <- nrow(data)
  nitem <- ncol(data)

  # preallocate the output list to avoid repeated growth
  freq <- vector("list", nitem)

  # build each item's one-hot indicator matrix independently because
  # cats[k] varies per item, so a single vectorized call is not feasible
  for (k in seq_len(nitem)) {

    # response vector for item k (length = nstd, may contain NAs)
    resp_k <- data[, k]

    # destination: one row per examinee, one column per score category;
    # 0L initializer also doubles as the value for NA-response rows
    m <- matrix(0L, nrow = nstd, ncol = cats[k])

    # row indices of examinees with an observed (non-NA) response
    not_na <- which(!is.na(resp_k))

    # set m[i, resp_k[i] + 1L] = 1 for each non-NA examinee.  Using a
    # two-column index matrix lets matrix-indexing assign all positions
    # in one shot without an inner loop.  +1L converts the 0-based
    # response value to a 1-based column index.
    if (length(not_na) > 0L) {
      m[cbind(not_na, as.integer(resp_k[not_na]) + 1L)] <- 1L
    }

    freq[[k]] <- m
  }

  # name the list by the column names of data, or X1..Xn when there are none
  cn <- colnames(data)
  if (is.null(cn)) cn <- paste0("X", seq_len(nitem))
  names(freq) <- cn

  freq
}


# Split the response data into sparse matrices of DRM correct responses,
# DRM incorrect responses, PRM category indicators, and all category indicators.
#' @importFrom Matrix Matrix
divide_data <- function(data, idx.item, freq.cat) {

  # divide the response data set into DRM and PRM parts
  if(!is.null(idx.item$idx.drm)) {
    data_drm_p <- Matrix::Matrix(data[, idx.item$idx.drm], sparse = TRUE)
    data_drm_q <- Matrix::Matrix(1 - data_drm_p, sparse = TRUE)
    data_drm_p[is.na(data_drm_p)] <- 0
    data_drm_q[is.na(data_drm_q)] <- 0
  } else {
    data_drm_p <- NULL
    data_drm_q <- NULL
  }
  if(!is.null(idx.item$idx.prm)) {
    data_prm <-
      Matrix::Matrix(do.call(what='cbind',
                             freq.cat[idx.item$idx.prm]), sparse = TRUE)
  } else {
    data_prm <- NULL
  }

  # create a response data matrix including all incorrect + correct responses
  data_all <- Matrix::Matrix(do.call(what = "cbind", freq.cat), sparse = TRUE)

  # return the results
  list(data_drm_p = data_drm_p, data_drm_q = data_drm_q,
       data_prm = data_prm, data_all = data_all)

}

# This function returns the column numbers of the frequency
# matrix corresponding to all items
cols4item <- function(nitem, cats, loc_1p_const=NULL) {

  cols.all <- vector('list', nitem)
  for(i in 1:nitem) {
    if(i == 1) {
      cat.st <- 1
      cat.ed <- cats[i]
    } else {
      cat.st <- cat.ed + 1
      cat.ed <- cat.st + (cats[i] - 1)
    }
    cols.tmp <- cat.st:cat.ed
    cols.all[[i]] <- cols.tmp
  }

  if(!is.null(loc_1p_const)) {
    cols.1pl <- unlist(cols.all[loc_1p_const])
  } else {
    cols.1pl <- NULL
  }

  # return
  rst <- list(cols.all = cols.all, cols.1pl=cols.1pl)
  rst

}


#' @export
coef.est_irt <- function(object, ...) {
  object$estimates
}

#' @export
coef.est_mg <- function(object, ...) {
  object$estimates
}

#' @export
coef.est_item <- function(object, ...) {
  object$estimates
}

#' @export
logLik.est_irt <- function(object, ...) {
  object$loglikelihood
}

#' @export
logLik.est_mg <- function(object, ...) {
  object$loglikelihood
}

#' @export
logLik.est_item <- function(object, ...) {
  object$loglikelihood
}

#' @export
vcov.est_irt <- function(object, ...) {
  object$covariance
}

#' @export
vcov.est_mg <- function(object, ...) {
  object$covariance
}

#' @export
vcov.est_item <- function(object, ...) {
  object$covariance
}

