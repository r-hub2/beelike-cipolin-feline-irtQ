#' Log-Likelihood of Ability Parameters
#'
#' This function computes the log-likelihood values for a set of ability
#' values, given item parameters and response data.
#'
#' @inheritParams rdif
#' @inheritParams est_score
#' @param x A data frame containing item metadata (e.g., item parameters, number
#'   of categories, IRT model types, etc.). See [irtQ::est_irt()] or
#'   [irtQ::simdat()] for more details about the item metadata. This data frame
#'   can be easily created using the [irtQ::shape_df()] function.
#' @param D A scaling constant used in IRT models to make the logistic function
#'   closely approximate the normal ogive function. A value of 1.702 is commonly
#'   used for this purpose. Default is 1.
#' @param theta A numeric vector of ability values at which to evaluate the
#'   log-likelihood function.
#' @param method A character string specifying the estimation method. Available
#'   options are:
#'   - `"ML"`: Maximum likelihood estimation
#'   - `"MLF"`: Maximum likelihood estimation with fences
#'   - `"MAP"`: Maximum a posteriori estimation
#'
#'   Default is `"ML"`.
#' @param norm.prior A numeric vector of length two specifying the mean and
#'   standard deviation of the normal prior distribution. Used only when
#'   `method = "MAP"`. Default is `c(0, 1)`.
#' @param fence.b A numeric vector of length two specifying the difficulty (*b*)
#'   parameters of the lower and upper fence items. If `NULL`, `c(-5, 5)` is
#'   used. Used only when `method = "MLF"`. Default is `NULL`.
#'
#' @details This function evaluates the log-likelihood at each value of
#' `theta` for one or more examinees, based on item parameters (`x`) and item
#' response data (`data`). For `method = "MAP"`, the log density of the normal
#' prior is added, so the values are log-posterior values.
#'
#' If `method = "MLF"`, two fence items with slope `fence.a` and difficulties
#' `fence.b` are appended, with the lower fence answered correctly and the
#' upper fence answered incorrectly. The fences give the likelihood a finite
#' maximum for all-correct and all-incorrect patterns. See [irtQ::est_score()]
#' for details.
#'
#' For example, to compute the log-likelihood curves of two examinees' responses
#' to the same test items, supply a 2-row matrix to `data` and a vector of
#' ability values to `theta`.
#'
#' @return A data frame of log-likelihood values (log-posterior values for
#'   `method = "MAP"`).
#' - Each **row** corresponds to a value of `theta`.
#' - Each **column** corresponds to an examinee and is named `Resp.1`,
#'   `Resp.2`, and so on. Examinees without any observed response have `NA` in
#'   their column.
#'
#' @examples
#' ## Import the "-prm.txt" output file from flexMIRT
#' flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
#'
#' # Read item parameters and convert them to item metadata
#' x <- bring.flexmirt(file = flex_sam, "par")$Group1$full_df
#'
#' # Generate ability values from N(0, 1)
#' set.seed(10)
#' score <- rnorm(5, mean = 0, sd = 1)
#'
#' # Simulate response data
#' data <- simdat(x = x, theta = score, D = 1)
#'
#' # Specify ability values for log-likelihood evaluation
#' theta <- seq(-3, 3, 0.5)
#'
#' # Compute log-likelihood values (ML)
#' llike_score(x = x, data = data, theta = theta, D = 1, method = "ML")
#'
#' @export
#' @importFrom reshape2 melt
#' @import dplyr
llike_score <- function(x,
                        data,
                        theta,
                        D = 1,
                        method = "ML",
                        norm.prior = c(0, 1),
                        fence.a = 3.0,
                        fence.b = NULL,
                        missing = NA) {

  # convert a single examinee's response vector to a one-row matrix
  if (is.vector(data)) {
    data <- rbind(data)
  }

  # drop row names so that duplicated names cannot merge examinees
  rownames(data) <- NULL

  # re-code missing values
  if (!is.na(missing)) {
    data[data == missing] <- NA
  }

  # check the number of examinees
  nstd <- nrow(data)

  # confirm and correct all item metadata information
  x <- confirm_df(x)

  # add two fence items and their responses when MLF is used
  if (method == "MLF") {
    # when fence.b = NULL, use the default ability range of est_score() as the fences
    if (is.null(fence.b)) {
      fence.b <- c(-5, 5)
    }

    # flag examinees without any observed response before adding the fences
    allmiss <- rowSums(!is.na(data)) == 0L

    # add two more response columns for the two fence items
    data <- cbind(data, f.lower = 1, f.upper = 0)

    # keep the fence responses missing for those examinees so that they get NA
    data[allmiss, ncol(data) - 1:0] <- NA

    # create item metadata for the two fence items
    x.fence <- shape_df(
      par.drm = list(a = rep(fence.a, 2), b = fence.b, g = rep(0, 2)),
      item.id = c("fence.lower", "fence.upper"), cats = 2,
      model = "3PLM"
    )

    # create the new item metadata by adding two fence items
    x <- dplyr::bind_rows(x, x.fence)
  }

  # check the maximum score category across all items
  max.cats <- max(x$cats)

  # count the number of columns of the item metadata
  max.col <- ncol(x)

  # create an empty elm_item object
  elm_item <- list(pars = NULL, model = NULL, cats = NULL)

  # reshape the data
  data <-
    data.frame(x, t(data), check.names = FALSE) %>%
    reshape2::melt(
      id.vars = 1:max.col,
      variable.name = "std",
      value.name = "resp",
      na.rm = TRUE,
      factorsAsStrings = FALSE
    )
  data$resp <- factor(data$resp, levels = (seq_len(max.cats) - 1))
  data$std <- as.numeric(data$std)

  # compute the log-likelihood values for all the discrete theta values
  lls <-
    purrr::map(
      .x = 1:nstd,
      .f = function(x) {
        exam_dat <-
          data %>%
          dplyr::filter(.data$std == x)
        llike_score_one(
          exam_dat = exam_dat, theta = theta,
          elm_item = elm_item, max.col = max.col, D = D,
          method = method, norm.prior = norm.prior
        )
      }
    ) %>%
    bind.fill(type = "cbind") %>%
    data.frame() %>%
    dplyr::rename_all(.f = ~ {
      paste0("Resp.", 1:nstd)
    })

  # return results
  lls
}


# Compute the log-likelihood values for one examinee
llike_score_one <- function(exam_dat, theta, elm_item, max.col, D = 1, method = "ML",
                            norm.prior = c(0, 1)) {
  # return NA when the examinee has no observed response
  if (nrow(exam_dat) == 0L) {
    return(rep(NA_real_, length(theta)))
  }

  # extract the required objects from the individual exam data
  elm_item$pars <- data.matrix(exam_dat[, 4:max.col])
  elm_item$model <- exam_dat$model
  elm_item$cats <- exam_dat$cats
  resp <- exam_dat$resp
  n.resp <- length(resp)

  # classify the items into DRM and PRM item groups
  idx.item <- idxfinder(elm_item)
  idx.drm <- idx.item$idx.drm
  idx.prm <- idx.item$idx.prm

  # calculate the score categories
  tmp.id <- 1:n.resp
  freq.cat <-
    matrix(
      stats::xtabs(~ tmp.id + resp,
        na.action = stats::na.pass, addNA = FALSE
      ),
      nrow = n.resp
    )

  # compute the negative log-likelihood values for all the discrete theta values
  ll_val <- ll_score(
    theta = theta, elm_item = elm_item, freq.cat = freq.cat, method = method,
    idx.drm = idx.drm, idx.prm = idx.prm, D = D, norm.prior = norm.prior,
    logL = TRUE
  )

  # return the log-likelihood values
  return(-ll_val)
}
