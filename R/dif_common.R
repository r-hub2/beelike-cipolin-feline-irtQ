# This function updates the ability estimates after an item is removed during the purification
# of the DIF and IPD functions. With the pattern scoring methods, the examinees are scored
# independently, so after the first iteration only the examinees who responded to the removed
# item are rescored when a single core is used
dif_rescore <- function(x, data, score, loc_resp, first, D, method, range,
                        norm.prior, nquad, weights, ncore, quiet_na = FALSE,
                        quiet_ncore = FALSE, se_prev = NULL, return_se = FALSE, ...) {
  # estimate the abilities and their standard errors, and mute the warning about the examinees
  # without any response when it is reported separately and the warning about the number of
  # cores when it was already given
  score_est <- function(data_est) {
    withCallingHandlers(
      est_score(
        x = x, data = data_est, D = D, method = method, range = range,
        norm.prior = norm.prior, nquad = nquad, weights = weights, ncore = ncore, ...),
      warning = function(w) {
        if ((quiet_na && grepl("^NA values are returned", conditionMessage(w))) ||
            (quiet_ncore && grepl("^ncore > 1 is not recommended", conditionMessage(w)))) {
          invokeRestart("muffleWarning")
        }
      }
    )
  }

  if (!first && ncore == 1 && method %in% c("ML", "MLF", "WL", "MAP", "EAP")) {
    # rescore the examinees who responded to the removed item
    if (length(loc_resp) > 0L) {
      score_rst <- score_est(data[loc_resp, , drop = FALSE])
      score[loc_resp] <- score_rst$est.theta
      se_prev[loc_resp] <- score_rst$se.theta
    }
  } else {
    # rescore all examinees
    score_rst <- score_est(data)
    score <- score_rst$est.theta
    se_prev <- score_rst$se.theta
  }

  # return the updated ability estimates, together with their standard errors when requested
  if (return_se) {
    list(score = score, se = se_prev)
  } else {
    score
  }
}


# This function stops when the group vector or the focal group is invalid for a two-group
# analysis of the DIF and IPD functions
check_group <- function(group, focal.name, nrow_data) {
  # stop when the group vector does not match the rows of the response data
  if (length(group) != nrow_data) {
    stop("The length of 'group' must equal the number of rows in 'data'.", call. = FALSE)
  }

  # stop when the focal group is not a single value found in the group vector
  if (length(focal.name) != 1L || !any(group == focal.name, na.rm = TRUE)) {
    stop("'focal.name' must be a single value found in 'group'.", call. = FALSE)
  }

  # stop when there is no examinee in the reference group
  if (!any(group != focal.name, na.rm = TRUE)) {
    stop("'group' must contain at least one examinee of the reference group.", call. = FALSE)
  }
}


# This function stops when the group vector or the focal groups are invalid for a multiple-group
# analysis, in which every value of 'group' not in 'focal.name' is the reference group
check_group_mg <- function(group, focal.name, nrow_data) {
  # stop when the group vector does not match the rows of the response data
  if (length(group) != nrow_data) {
    stop("The length of 'group' must equal the number of rows in 'data'.", call. = FALSE)
  }

  # stop when the focal groups are not distinct values found in the group vector
  if (length(focal.name) == 0L || anyNA(focal.name) || anyDuplicated(focal.name) > 0L ||
      !all(focal.name %in% group)) {
    stop("'focal.name' must contain distinct values, each found in 'group'.", call. = FALSE)
  }

  # find the groups that are not focal groups
  ref_name <- unique(group[!is.na(group) & !(group %in% focal.name)])

  # stop when there is no examinee in the reference group
  if (length(ref_name) == 0L) {
    stop("'group' must contain at least one examinee of the reference group.", call. = FALSE)
  }

  # stop when more than one group is left as the reference group
  if (length(ref_name) > 1L) {
    stop("'group' must contain exactly one reference group, which is the group not included ",
         "in 'focal.name'.", call. = FALSE)
  }
}


# This function drops the item parameter columns that no remaining item uses, which are the
# columns beyond the larger of three and the maximum number of score categories
trim_par_cols <- function(x) {
  # find the parameter columns that are not used by any item
  n_par <- max(3, max(x$cats))
  par_col <- grep("^par\\.", names(x))
  drop_col <- par_col[seq_along(par_col) > n_par]

  # remove the unused columns
  if (length(drop_col) > 0L) {
    x <- x[, -drop_col, drop = FALSE]
  }
  x
}


# This function checks the inputs of rdif(), crdif(), and grdif(), recodes the missing responses, and
# returns the response matrix, the ability estimates, and the positions of the items to be skipped
dif_prepare <- function(x, data, score, group, focal.name, item.skip, D, alpha, missing,
                        purify, max.iter, min.resp, method, range, norm.prior, nquad,
                        weights, ncore, fn, multi_focal = FALSE, ...) {
  # transform the response data to a matrix form
  data <- as.matrix(data)

  # re-code missing values
  if (!is.na(missing)) {
    data[data == missing] <- NA
  }

  # transform the response data to a numeric matrix and check the responses
  data <- resp_to_matrix(data, x$cats, x$id)

  if (multi_focal) {
    # check the group vector and the focal groups of a multiple-group analysis
    check_group_mg(group, focal.name, nrow(data))
  } else {
    # check the group vector and the focal group
    check_group(group, focal.name, nrow(data))

    # stop when the group vector has more than two groups
    if (length(unique(group[!is.na(group)])) > 2L) {
      stop("'group' must contain exactly two groups (the reference and focal groups). ",
           "Use grdif() for more than two groups.", call. = FALSE)
    }
  }

  # stop when the significance level is not a single value between 0 and 1
  check_alpha(alpha)

  # stop when the maximum number of iterations is not a positive whole number
  check_max_iter(max.iter, purify)

  # stop when the minimum number of responses is not a single nonnegative number
  check_min_resp(min.resp)

  # stop when the abilities are estimated or updated with an unsupported scoring method
  if (is.null(score) || purify) {
    check_pattern_method(method, fn = fn)
  }

  # check the positions of the items to be skipped
  item.skip <- check_item_skip(item.skip, nrow(x))

  # if min.resp is not NULL, find the examinees who have the number of responses
  # less than specified value (e.g., 5). Then, replace their all responses with NA
  loc_less <- integer(0)
  if (!is.null(min.resp)) {
    n_resp <- rowSums(!is.na(data))
    loc_less <- which(n_resp < min.resp & n_resp > 0)

    # warn that the examinees with too few responses are excluded
    if (length(loc_less) > 0L) {
      warning(length(loc_less), " examinee(s) with fewer than ", min.resp,
              " responses were excluded.", call. = FALSE)
    }
    data[loc_less, ] <- NA
  }

  # compute the score if score = NULL
  if (!is.null(score)) {
    # transform the scores to a vector form and check them
    score <- check_score(score, nrow(data))

    # treat the ability estimates of the examinees excluded by min.resp as missing
    score[loc_less] <- NA
  } else {
    # estimate the abilities; the examinees without a response are excluded without a warning
    score <- dif_rescore(
      x = x, data = data, score = NULL, loc_resp = NULL, first = TRUE, D = D, method = method,
      range = range, norm.prior = norm.prior, nquad = nquad, weights = weights, ncore = ncore,
      quiet_na = TRUE, ...
    )
  }

  # stop when no examinee has an ability estimate
  if (all(is.na(score[!is.na(group)]))) {
    stop("No examinee has an ability estimate.", call. = FALSE)
  }

  # exclude the examinees with responses but without an ability estimate, with a warning
  na_score <- is.na(score) & rowSums(!is.na(data)) > 0
  if (any(na_score)) {
    warning(sum(na_score), " examinee(s) with a missing ability estimate were excluded.",
            call. = FALSE)
    data[na_score, ] <- NA
  }

  # return the checked inputs
  list(data = data, score = score, item.skip = item.skip)
}


# This function chooses the item to be removed in a purification iteration: among the flagged
# items, the one with the smallest p-value, which is compared on the log scale to avoid ties at zero
pick_flagged_item <- function(flag_loc, log_p) {
  flag_loc[which.min(log_p[flag_loc])]
}


# This function warns once about the examinees who lost their ability estimates during the
# purification of the DIF and IPD functions
warn_purify_excluded <- function(n_excluded) {
  if (n_excluded > 0L) {
    warning(n_excluded, " examinee(s) were excluded during the purification because they had ",
            "too few responses to the remaining items to estimate the ability.", call. = FALSE)
  }
}


# This function stops when the minimum number of responses is not NULL or a single nonnegative
# number
check_min_resp <- function(min.resp) {
  if (!is.null(min.resp) && !(is.numeric(min.resp) && length(min.resp) == 1L &&
                              is.finite(min.resp) && min.resp >= 0)) {
    stop("'min.resp' must be NULL or a single nonnegative number.", call. = FALSE)
  }
}


# This function transforms the supplied ability estimates to a vector and stops when they are not
# numbers or do not match the rows of the response data
check_score <- function(score, nrow_data) {
  # transform the scores to a vector form
  if (is.matrix(score) | is.data.frame(score)) {
    score <- as.numeric(data.matrix(score))
  }

  # stop when the scores are not numbers
  if (!is.numeric(score)) {
    stop("'score' must be a numeric vector.", call. = FALSE)
  }

  # stop when the scores do not match the rows of the response data
  if (length(score) != nrow_data) {
    stop("The length of 'score' must equal the number of rows in 'data'.", call. = FALSE)
  }

  # return the scores
  score
}


# This function stops when the significance level is not a single value between 0 and 1
check_alpha <- function(alpha) {
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single number between 0 and 1.", call. = FALSE)
  }
}


# This function stops when the purification is requested and the maximum number of iterations is
# not a positive whole number
check_max_iter <- function(max.iter, purify) {
  if (purify && !(is.numeric(max.iter) && length(max.iter) == 1L && is.finite(max.iter) &&
                  max.iter >= 1 && max.iter == round(max.iter))) {
    stop("'max.iter' must be a single positive whole number when purify = TRUE.", call. = FALSE)
  }
}


# This function stops when the critical value is not NULL or a single nonnegative number
check_crit_val <- function(crit.val) {
  if (!is.null(crit.val) && !(is.numeric(crit.val) && length(crit.val) == 1L &&
                              is.finite(crit.val) && crit.val >= 0)) {
    stop("'crit.val' must be NULL or a single nonnegative number.", call. = FALSE)
  }
}


# This function stops when the scoring method does not give one ability estimate per response
# pattern, which the DIF and IPD functions need to estimate or update the abilities
check_pattern_method <- function(method, fn) {
  if (length(method) != 1L || !method %in% c("ML", "MLF", "WL", "MAP", "EAP")) {
    stop("The scoring method ", paste(method, collapse = ", "), " is not supported by ", fn, ". ",
         "Use one of \"ML\", \"MLF\", \"WL\", \"MAP\", or \"EAP\".", call. = FALSE)
  }
}


# This function stops when the response data, the scaling factor D, or (for an est_item object)
# the ability estimates are passed through '...' to a method for an estimation object, which
# takes them from the object
check_obj_dots <- function(dots, obj_score = FALSE) {
  if (any(c("data", "D") %in% names(dots))) {
    stop("The response data and the scaling factor D are taken from the object 'x'.",
         call. = FALSE)
  }
  if (obj_score && "score" %in% names(dots)) {
    stop("The response data, the ability estimates, and the scaling factor D are taken from ",
         "the object 'x'.", call. = FALSE)
  }
}


# This function computes the quadratic form of a deviation vector in the inverse of a covariance
# matrix. When the matrix is singular or nearly singular, it uses the Moore-Penrose generalized
# inverse and returns the rank of the matrix as the degrees of freedom
quad_form <- function(cov_mat, dev_vec, df) {
  # return a missing statistic when the covariance matrix has non-finite values
  if (!all(is.finite(cov_mat))) {
    return(list(stat = NA_real_, df = df, reduced = FALSE))
  }

  # compute the reciprocal condition number of the covariance matrix
  rc_cov <- tryCatch(rcond(cov_mat), error = function(e) 0)

  # use the inverse of the covariance matrix when it is well conditioned
  if (rc_cov > 1e-10) {
    stat <- as.numeric(t(dev_vec) %*% solve(cov_mat) %*% dev_vec)
    return(list(stat = stat, df = df, reduced = FALSE))
  }

  # otherwise use the generalized inverse and the rank of the covariance matrix
  eig <- eigen(cov_mat, symmetric = TRUE)
  keep <- eig$values > max(eig$values) * 1e-10
  if (!any(keep)) {
    return(list(stat = NA_real_, df = df, reduced = FALSE))
  }
  proj <- crossprod(eig$vectors[, keep, drop = FALSE], dev_vec)
  list(stat = sum(proj^2 / eig$values[keep]), df = sum(keep), reduced = sum(keep) < df)
}
