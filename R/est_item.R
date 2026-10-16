#' Fixed ability parameter calibration
#'
#' This function performs fixed ability parameter calibration (FAPC), often
#' called Stocking's (1988) Method A, which estimates item parameters by
#' maximum likelihood given ability estimates (Baker & Kim, 2004; Ban et al.,
#' 2001; Stocking, 1988). When prior distributions are specified for item
#' parameters (by default, for the guessing parameters of 3PLM items), the
#' estimates are the modes of the posterior distribution (Bayes modal
#' estimates). It can also be considered a special case of joint maximum
#' likelihood estimation in which only one cycle of item parameter estimation
#' is conducted, conditioned on the given ability estimates (Birnbaum, 1968).
#' FAPC is a potentially useful method for calibrating pretest
#' (or newly developed) items in computerized adaptive testing (CAT), as it
#' enables placing their parameter estimates on the same scale as operational
#' items. In addition, it can be used to recalibrate operational items in the
#' item bank to evaluate potential parameter drift (Chen & Wang, 2016; Stocking,
#' 1988).
#'
#' @inheritParams est_irt
#' @param x A data frame containing item metadata. This metadata is required to
#'   retrieve essential information for each item (e.g., number of score
#'   categories, IRT model type, etc.) necessary for calibration. You can create
#'   an empty item metadata frame using the function [irtQ::shape_df()].
#'
#'   When `use.startval = TRUE`, the item parameters specified in the metadata
#'   will be used as starting values for parameter estimation; otherwise they
#'   are ignored. If `x = NULL`, the models and numbers of score categories are
#'   taken from the `model` and `cats` arguments. See [irtQ::est_irt()]
#'   or [irtQ::simdat()] for more details about the item metadata. Default is
#'   `NULL`.
#' @param score A numeric vector of examinees' ability estimates (theta values),
#'   which are treated as known. It must have one finite value per row of
#'   `data`, in the same order as the rows.
#' @param model A character vector specifying the IRT model to fit each item.
#'   Available values are:
#'   - `"1PLM"`, `"2PLM"`, `"3PLM"`, `"DRM"` for dichotomous items
#'   - `"GRM"`, `"GPCM"` for polytomous items
#'
#'   Here, `"GRM"` denotes the graded response model and `"GPCM"` the
#'   (generalized) partial credit model. Note that `"DRM"` is a general label
#'   for dichotomous items; items labeled `"DRM"` are fitted with the 3PLM, and
#'   a warning is issued. If a single model name is provided, it is recycled for
#'   all items. This argument is only used when `x = NULL`. Default is `NULL`.
#' @param cats Numeric vector specifying the number of score categories per
#'   item. For dichotomous items, this should be 2. If a single value is
#'   supplied, it will be recycled across all items. When `cats = NULL` and all
#'   models specified in the `model` argument are dichotomous (`"1PLM"`,
#'   `"2PLM"`, `"3PLM"`, or `"DRM"`), the function defaults to 2 categories per
#'   item. This argument is used only when `x = NULL`. Default is `NULL`.
#' @param item.id Character vector of item identifiers. If `NULL`, the IDs in
#'   `x` are used, or IDs are generated automatically when `x = NULL`. A
#'   provided `item.id` replaces any IDs in `x` and must have one value per
#'   item. Default is `NULL`.
#' @param control A named list of control options for [stats::nlminb()], which
#'   maximizes the objective of each item (see Details). The default is
#'   `list(eval.max = 500, iter.max = 200, x.tol = 1e-4)`, where
#'   - `eval.max` is the maximum number of objective function evaluations,
#'   - `iter.max` is the maximum number of optimizer iterations, and
#'   - `x.tol` is the X-convergence tolerance of [stats::nlminb()].
#'
#'   A partial list may be supplied; options that are not supplied keep their
#'   default values. Other options of [stats::nlminb()], such as `rel.tol` or
#'   `trace`, can also be supplied.
#' @param verbose Logical. If `FALSE`, all progress messages are suppressed.
#'   Default is `TRUE`.
#'
#' @details The parameters of each item are estimated separately by maximizing
#'   the log-likelihood of its responses given the abilities in `score` (plus
#'   the log prior densities when priors are used), except that the 1PLM items
#'   share one slope, estimated jointly with their difficulties, when
#'   `fix.a.1pl = FALSE`. Missing responses are omitted. Items without any
#'   response are excluded from the estimation with a warning.
#'
#'   In most cases, the function [irtQ::est_item()] returns successfully
#'   converged item parameter estimates using its default internal starting
#'   values. However, if convergence issues arise during calibration, one
#'   possible solution is to use alternative starting values. If item parameter
#'   values are already specified in the item metadata (i.e., the `x` argument),
#'   they can be used as starting values for item parameter calibration by
#'   setting `use.startval = TRUE`.
#'
#'   A maximum likelihood estimate does not exist for an item with a score
#'   category that has no responses, which includes an item whose responses are
#'   all the same. The estimates of such an item can be extreme or remain near
#'   their starting values, and their standard errors are unreliable (a standard
#'   error of 99999 means that it could not be obtained). After the items
#'   without any response are excluded, one warning lists the IDs of the
#'   remaining items of this kind. A parameter with a prior distribution, such
#'   as the guessing parameter of a 3PLM item with `use.gprior = TRUE`, can
#'   still be stable.
#'
#'   When prior distributions are used, the item parameters are Bayes modal
#'   estimates, and `loglikelihood` is the log-likelihood at these estimates
#'   without the log prior densities.
#'
#'   The print method of this class has an argument `digits` that is not used.
#'   It prints the log-likelihood with the number of significant digits given
#'   by `getOption("digits")` (7 by default).
#'
#' @return This function returns an object of class `est_item`. The returned
#'   object contains the following components:
#'
#'   \item{estimates}{A data frame containing both the item parameter estimates
#'   and their corresponding standard errors.}
#'
#'   \item{par.est}{A data frame of item parameter estimates, structured
#'   according to the item metadata format.}
#'
#'   \item{se.est}{A data frame of standard errors for the item parameter
#'   estimates, computed from the inverse of the observed information matrix
#'   (the Hessian of the negative log-likelihood, plus the negative second
#'   derivatives of the log prior densities when priors are used). The abilities
#'   in `score` are treated as known, so their estimation error is not reflected
#'   in the standard errors.}
#'
#'   \item{pos.par}{A data frame of the positions of the estimated item
#'   parameters in `covariance`. The rows of the excluded items are `NA`, and a
#'   slope shared by 1PLM items appears only in the row of the first 1PLM
#'   item.}
#'
#'   \item{covariance}{A variance-covariance matrix of the item parameter
#'   estimates. Items are estimated separately, so the covariances between the
#'   parameters of different items are 0, except among the 1PLM items that
#'   share a common slope.}
#'
#'   \item{loglikelihood}{The log-likelihood of the observed responses given
#'   the abilities in `score`, summed over the items, at the final estimates.
#'   The log prior densities of the item parameters are not included. Printed
#'   results report -2 times this value as "-2loglikelihood".}
#'
#'   \item{group.par}{A named numeric vector with the mean (`mu`) and the
#'   standard deviation (`sigma`) of the abilities in `score`.}
#'
#'   \item{data}{A data frame of the response data used in the estimation,
#'   with missing responses as `NA` and without the excluded items (see
#'   `deleted.item`).}
#'
#'   \item{score}{A vector of examinees' ability estimates used as fixed values
#'   during item parameter estimation.}
#'
#'   \item{scale.D}{The scaling factor used in the IRT model.}
#'
#'   \item{convergence}{A message indicating whether item parameter estimation
#'   successfully converged.}
#'
#'   \item{nitem}{The number of items in `data` as supplied, including the
#'   excluded items.}
#'
#'   \item{deleted.item}{A numeric vector of the positions of the items without
#'   any response, which are excluded from the estimation (empty if there are
#'   none). Their rows in `par.est` and `se.est` are `NA`.}
#'
#'   \item{npar.est}{The number of estimated item parameters.}
#'
#'   \item{n.response}{A numeric vector of the numbers of non-missing responses
#'   to the items, with 0 for the excluded items.}
#'
#'   \item{TotalTime}{Total computation time in seconds.}
#'
#'   \item{call}{The matched call.}
#'
#'   Note that you can easily extract components from the output using the
#'   [irtQ::getirt()] function.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [irtQ::est_irt()], [irtQ::shape_df()], [irtQ::getirt()]
#'
#' @references Baker, F. B., & Kim, S. H. (2004). *Item response theory:
#'   Parameter estimation techniques.* CRC Press.
#'
#'   Ban, J. C., Hanson, B. A., Wang, T., Yi, Q., & Harris, D. J. (2001). A
#'   comparative study of on-line pretest item calibration/scaling methods in
#'   computerized adaptive testing. *Journal of Educational Measurement, 38*(3),
#'   191-212. \doi{10.1111/j.1745-3984.2001.tb01123.x}.
#'
#'   Birnbaum, A. (1968). Some latent trait models and their use in inferring an
#'   examinee's ability. In F. M. Lord & M. R. Novick (Eds.),
#'   *Statistical theories of mental test scores* (pp. 397-479). Reading, MA: Addison-Wesley.
#'
#'   Chen, P., & Wang, C. (2016). A new online calibration method for
#'   multidimensional computerized adaptive testing.
#'   *Psychometrika, 81*(3), 674-701.
#'
#'   Stocking, M. L. (1988). *Scale drift in on-line calibration* (Research Rep.
#'   88-28). Princeton, NJ: ETS.
#'
#' @examples
#' ## Import the "-prm.txt" output file from flexMIRT
#' flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
#'
#' # Extract the item metadata
#' x <- bring.flexmirt(file = flex_sam, "par")$Group1$full_df
#'
#' # Modify the item metadata so that some items follow 1PLM, 2PLM, and GPCM
#' x[c(1:3, 5), 3] <- "1PLM"
#' x[c(1:3, 5), 4] <- 1
#' x[c(1:3, 5), 6] <- 0
#' x[c(4, 8:12), 3] <- "2PLM"
#' x[c(4, 8:12), 6] <- 0
#' x[54:55, 3] <- "GPCM"
#'
#' # Generate examinees' abilities from N(0, 1)
#' set.seed(23)
#' score <- rnorm(500, mean = 0, sd = 1)
#'
#' # Simulate response data based on the item metadata and ability values
#' data <- simdat(x = x, theta = score, D = 1)
#'
#' \donttest{
#' # 1) Estimate item parameters: constrain the slope parameters of 1PLM items
#' #    to be equal
#' (mod1 <- est_item(x, data, score,
#'   D = 1, fix.a.1pl = FALSE, use.gprior = TRUE,
#'   gprior = list(dist = "beta", params = c(5, 16)), use.startval = FALSE
#' ))
#' summary(mod1)
#'
#' # Extract the item parameter estimates
#' getirt(mod1, what = "par.est")
#'
#' # 2) Estimate item parameters: fix the slope parameters of 1PLM items to 1
#' (mod2 <- est_item(x, data, score,
#'   D = 1, fix.a.1pl = TRUE, a.val.1pl = 1, use.gprior = TRUE,
#'   gprior = list(dist = "beta", params = c(5, 16)), use.startval = FALSE
#' ))
#' summary(mod2)
#'
#' # Extract the standard error estimates
#' getirt(mod2, what = "se.est")
#'
#' # 3) Estimate item parameters: fix the guessing parameters of 3PLM items to 0.2
#' (mod3 <- est_item(x, data, score,
#'   D = 1, fix.a.1pl = TRUE, fix.g = TRUE, a.val.1pl = 1, g.val = .2,
#'   use.startval = FALSE
#' ))
#' summary(mod3)
#'
#' # Extract both item parameter and standard error estimates
#' getirt(mod3, what = "estimates")
#'
#' }
#'
#' @importFrom Matrix sparseMatrix
#' @importFrom utils modifyList
#' @import dplyr
#' @export
#'
est_item <- function(x = NULL,
                     data,
                     score,
                     D = 1,
                     model = NULL,
                     cats = NULL,
                     item.id = NULL,
                     fix.a.1pl = FALSE,
                     fix.a.gpcm = FALSE,
                     fix.g = FALSE,
                     a.val.1pl = 1,
                     a.val.gpcm = 1,
                     g.val = .2,
                     use.aprior = FALSE,
                     use.bprior = FALSE,
                     use.gprior = TRUE,
                     aprior = list(dist = "lnorm", params = c(0, 0.5)),
                     bprior = list(dist = "norm", params = c(0.0, 1.0)),
                     gprior = list(dist = "beta", params = c(5, 16)),
                     missing = NA,
                     use.startval = FALSE,
                     control = list(eval.max = 500, iter.max = 200, x.tol = 1e-4),
                     verbose = TRUE) {

  # check start time
  start.time <- Sys.time()

  # match.call
  cl <- match.call()

  # Merge user-supplied control list with defaults, enabling partial specification
  # (e.g., control = list(iter.max = 500) keeps eval.max and x.tol at defaults)
  default_control <- list(eval.max = 500, iter.max = 200, x.tol = 1e-4)
  control <- modifyList(default_control, control)

  ## -------------------------------------------------------------------------------------------------------
  ## 1. preparation of item parameter estimation
  if (verbose) {
    cat("Starting...", "\n")
  }

  # check if the starting values are available
  if (use.startval & is.null(x)) {
    stop(paste0(
      "To use starting values for item parameter estimation, \n",
      "the item metadata must be specified in the argument 'x'."
    ), call. = FALSE)
  }

  # transform a data set to data.frame
  if (is.vector(data)) {
    data <- data.frame(X1 = data, stringsAsFactors = FALSE)
  } else {
    data <- data.frame(data, stringsAsFactors = FALSE)
  }

  # extract information about the number of score categories and models
  if (verbose) {
    cat("Parsing input...", "\n")
  }
  if (!is.null(x)) {
    # confirm and correct all item metadata information
    x <- confirm_df(x)

    # override the item IDs in x when item.id is supplied
    if (!is.null(item.id)) {
      # stop when item.id does not have one value per item of x
      if (length(item.id) != nrow(x)) {
        stop("The length of 'item.id' must be equal to the number of items in 'x'.", call. = FALSE)
      }
      x$id <- item.id
    }

    # extract score categories for all items
    cats <- x$cats

    # extract model names for all items
    model <- x$model
  } else {
    # make the model names as upper cases
    model <- toupper(model)

    # when a character string scalar is provided in the model argument
    # create a vector of the same model names
    if (length(model) == 1) {
      model <- rep(model, ncol(data))
    }

    # check if the score category information is provided in the cats argument
    if (is.null(cats)) {
      if (all(model %in% c("1PLM", "2PLM", "3PLM", "DRM"))) {
        cats <- rep(2, ncol(data))
      } else {
        stop("The number of score categories for the items should be specified in the argument 'cats'.", call. = FALSE)
      }
    }

    # when an integer scalar is provided in the cats argument
    # create a vector of the same category numbers
    if (length(cats) == 1) {
      cats <- rep(cats, ncol(data))
    }

    # create an empty item metadata data frame
    x <- shape_df(item.id = item.id, cats = cats, model = model, default.par = TRUE)

    # confirm and correct all item metadata information
    x <- confirm_df(x)

    # extract score categories for all items
    cats <- x$cats

    # extract model names for all items
    model <- x$model
  }

  # check the total number of item in the response data set
  nitem <- ncol(data)

  # check the total number of examinees
  nstd <- nrow(data)

  # check whether included data are correct
  if (nrow(x) != nitem) stop("The number of items included in 'x' and 'data' must be the same.", call. = FALSE)

  # recode missing values
  if (!is.na(missing)) {
    data[data == missing] <- NA
  }

  # stop when a response is not one of the score categories 0, ..., cats - 1
  # (character and factor columns are checked by their values)
  check_resp_range(data = resp_as_numeric(data), cats = cats, id = x$id)

  # transform scores to a vector if the score has other formats
  if (is.matrix(score) | is.data.frame(score)) {
    score <- as.numeric(data.matrix(score))
  }

  # stop when the scores do not give one finite ability value per examinee
  if (length(score) != nstd || !all(is.finite(score))) {
    stop("'score' must contain one finite ability value (no missing value) for each examinee (row) of 'data'.",
         call. = FALSE)
  }

  # check the number of item responses across all items
  n.resp <- colSums(!is.na(data))

  # check the items which have all missing responses
  loc_allmiss <- which(n.resp == 0L)

  # delete items which have all missing responses from the item metadata set
  if (length(loc_allmiss) > 0L) {
    loc_nomiss <- which(n.resp > 0L)
    x_pre <- x
    x <- x[-loc_allmiss, ]
    model <- x[, 3]
    cats <- x[, 2]
    # keep the data a data frame even when one item is left
    data <- data[, loc_nomiss, drop = FALSE]
    memo2 <- paste0(
      paste0("item ", loc_allmiss, collapse = ", "),
      " is/are excluded in the item parameter estimation because the item(s) has/have no item response data."
    )
    warning(memo2, call. = FALSE)
  } else {
    loc_allmiss <- NULL
  }

  # warn about the items whose maximum likelihood estimates do not exist
  check_item_mle(data = resp_as_numeric(data), cats = cats, id = x$id)

  # compute the mean and standard deviation of the abilities
  group.par <- c(mu = mean(score), sigma = stats::sd(score))

  # copy scores to theta values
  theta <- score

  # find the location of 1PLM items in which slope parameters should be constrained to be equal
  # also, find the location of items with other models
  if ("1PLM" %in% model & !fix.a.1pl) {
    loc_1p_const <- which(model == "1PLM")
    loc_else <- which(model != "1PLM")

    # count the number of 1PLM items to be constrained
    n.1PLM <- length(loc_1p_const)
  } else {
    loc_1p_const <- NULL
    loc_else <- 1:nrow(x)
    n.1PLM <- NULL
  }

  # record the original location of item parameters to be estimated, and
  # the relocated position of item parameters when computing
  # the variance-covariance matrix of item parameter estimates
  param_loc <- parloc(
    x = x, loc_1p_const = loc_1p_const, loc_else = loc_else,
    fix.a.1pl = fix.a.1pl, fix.a.gpcm = fix.a.gpcm, fix.g = fix.g
  )

  # build the per-item one-hot frequency-category list used downstream
  # in the FAPC scoring loop.  See build_freqcat() (R/util.R) for the
  # output structure.
  freq.cat <- build_freqcat(data, cats)

  ## -------------------------------------------------------------------------------------------------------
  ## 2. item parameter estimation
  ## ---------------------------------------------------------------
  # Item parameter estimation
  # create empty vectors to contain results
  est_par <- NULL
  est_se <- NULL
  loc_items <- NULL
  convergence <- NULL
  noconv_items <- NULL
  objective <- NULL
  cov_list <- NULL

  # break down the item metadata into several elements
  elm_item <- breakdown(x)

  # classify the items into DRM and PRM item groups
  idx.item <- idxfinder(elm_item)
  idx.drm <- idx.item$idx.drm
  idx.prm <- idx.item$idx.prm

  # then, divide items into three groups
  # (1) DRM 1PL items with the constrained slope
  # (2) other DRM items
  # (3) PRM items
  # (no item other than the constrained 1PLM items when loc_else is empty)
  if (length(loc_else) == 0) {
    drm.else <- NULL
  } else {
    drm.else <- loc_else
  }
  idx4est <- list(
    drm.slc = loc_1p_const,
    drm.else = drm.else,
    prm = idx.prm
  )

  # create the lower and upper bounds of the item parameters
  parbd <- lubound(model, cats, n.1PLM, idx4est, fix.a.1pl, fix.g, fix.a.gpcm)

  # start estimation
  if (verbose) {
    cat("Estimating item parameters...", "\n")
  }

  # (1) estimation of DRM items: 1PLM with the constrained slope value
  if (!is.null(loc_1p_const)) {
    # prepare input files to estimate the 1PLM item parameters.
    # extract the per-item (incorrect, correct) columns from freq.cat
    # via array preallocation (so the resulting matrices have no
    # dimnames) and compute f_i = s_i + r_i in one elementwise add.
    fc1pl <- freq.cat[loc_1p_const]
    s_i <- r_i <- array(0, c(nstd, n.1PLM))
    for (k in seq_len(n.1PLM)) {
      s_i[, k] <- fc1pl[[k]][, 1]
      r_i[, k] <- fc1pl[[k]][, 2]
    }
    f_i <- s_i + r_i

    # set the starting values
    if (use.startval) {
      startval <-
        set_startval(
          pars = elm_item$pars, item = loc_1p_const,
          use.startval = TRUE, mod = "1PLM",
          score.cat = 2, fix.a.1pl = FALSE, fix.g = fix.g,
          fix.a.gpcm = fix.a.gpcm, n.1PLM = n.1PLM
        )
    } else {
      startval <-
        set_startval(
          pars = elm_item$pars, item = loc_1p_const,
          use.startval = FALSE, mod = "1PLM",
          score.cat = 2, fix.a.1pl = FALSE, fix.g = fix.g,
          fix.a.gpcm = fix.a.gpcm, n.1PLM = n.1PLM
        )
    }

    # bounds of the item parameters
    lower <- parbd$drm.slc$lower
    upper <- parbd$drm.slc$upper

    # parameter estimation
    est <- estimation1(
      f_i = f_i, r_i = r_i, s_i, theta = theta, mod = "1PLM", D = D, nstd = nstd,
      fix.a.1pl = FALSE, n.1PLM = n.1PLM, aprior = aprior, bprior = bprior,
      use.aprior = use.aprior, use.bprior = use.bprior,
      control = control, startval = startval, lower = lower, upper = upper
    )

    # extract the results
    # item parameter estimates
    a <- est$pars[1]
    b <- est$par[-1]
    pars <- purrr::map(1:n.1PLM, .f = function(x) c(a, b[x], NA))
    est_par <- c(est_par, pars)

    # standard errors
    a.se <- est$se[1]
    b.se <- est$se[-1]
    pars.se <- purrr::map(1:n.1PLM, .f = function(x) c(NA, b.se[x], NA))
    pars.se[[1]][1] <- a.se
    est_se <- c(est_se, pars.se)

    # covariance matrix
    cov_list <- c(cov_list, list(est$covariance))

    # convergence indicator
    convergence <- c(convergence, est$convergence)
    if (est$convergence > 0L) noconv_items <- c(noconv_items, loc_1p_const)

    # negative log-likelihood value without the log prior densities
    objective <- c(objective, est$objective.data)
  }

  # all other items
  if (length(loc_else) >= 1) {
    for (i in 1:length(loc_else)) {
      # prepare information to estimate item parameters
      mod <- model[loc_else][i]
      score.cat <- cats[loc_else][i]

      # in case of a DRM item
      if (score.cat == 2) {
        # extract (incorrect, correct) columns from the raw freq.cat
        # matrix for this item; cast to double so that s_i, r_i, and f_i
        # are double vectors for downstream estimation1()
        fc_i <- freq.cat[loc_else][[i]]
        s_i <- as.double(fc_i[, 1])
        r_i <- as.double(fc_i[, 2])
        # f_i (total response indicator: 1 if examinee answered the
        # item, 0 otherwise)
        f_i <- s_i + r_i

        # set the starting values
        if (use.startval) {
          startval <-
            set_startval(
              pars = elm_item$pars, item = loc_else[i],
              use.startval = TRUE, mod = mod,
              score.cat = score.cat, fix.a.1pl = TRUE, fix.g = fix.g,
              fix.a.gpcm = fix.a.gpcm, n.1PLM = NULL
            )
        } else {
          startval <-
            set_startval(
              pars = elm_item$pars, item = loc_else[i],
              use.startval = FALSE, mod = mod,
              score.cat = score.cat, fix.a.1pl = TRUE, fix.g = fix.g,
              fix.a.gpcm = fix.a.gpcm, n.1PLM = NULL
            )
        }

        # bounds of the item parameters
        loc.tmp <- which(idx4est$drm.else == loc_else[i])
        lower <- parbd$drm.else[[loc.tmp]]$lower
        upper <- parbd$drm.else[[loc.tmp]]$upper

        # parameter estimation
        est <- estimation1(
          f_i = f_i, r_i = r_i, s_i = s_i, theta = theta, mod = mod, D = D,
          nstd = nstd, fix.a.1pl = ifelse(mod == "1PLM", TRUE, FALSE),
          fix.g = fix.g, a.val.1pl = a.val.1pl, g.val = g.val, n.1PLM = NULL,
          aprior = aprior, bprior = bprior, gprior = gprior,
          use.aprior = use.aprior, use.bprior = use.bprior, use.gprior = use.gprior,
          control = control, startval = startval, lower = lower, upper = upper
        )

        # extract the results
        # item parameter estimates
        a <- ifelse(mod == "1PLM", a.val.1pl, est$pars[1])
        b <- ifelse(mod == "1PLM", est$pars[1], est$pars[2])
        g <- ifelse(mod == "3PLM", ifelse(fix.g, g.val, est$pars[3]), NA)
        pars <- c(a, b, g)
        est_par <- c(est_par, list(pars))

        # standard errors
        a.se <- ifelse(mod == "1PLM", NA, est$se[1])
        b.se <- ifelse(mod == "1PLM", est$se[1], est$se[2])
        g.se <- ifelse(mod == "3PLM", ifelse(fix.g, NA, est$se[3]), NA)
        pars.se <- c(a.se, b.se, g.se)
        est_se <- c(est_se, list(pars.se))

        # covariance matrix
        cov_list <- c(cov_list, list(est$covariance))

        # convergence indicator
        convergence <- c(convergence, est$convergence)
        if (est$convergence > 0L) noconv_items <- c(noconv_items, loc_else[i])

        # negative log-likelihood value without the log prior densities
        objective <- c(objective, est$objective.data)
      }

      # in case of a PRM item
      if (score.cat > 2) {
        r_i <- freq.cat[loc_else][[i]]

        # set the starting values
        if (use.startval) {
          startval <-
            set_startval(
              pars = elm_item$pars, item = loc_else[i],
              use.startval = TRUE, mod = mod,
              score.cat = score.cat, fix.a.1pl = fix.a.1pl, fix.g = fix.g,
              fix.a.gpcm = fix.a.gpcm, n.1PLM = NULL
            )
        } else {
          startval <-
            set_startval(
              pars = elm_item$pars, item = loc_else[i],
              use.startval = FALSE, mod = mod,
              score.cat = score.cat, fix.a.1pl = fix.a.1pl, fix.g = fix.g,
              fix.a.gpcm = fix.a.gpcm, n.1PLM = NULL
            )
        }

        # bounds of the item parameters
        loc.tmp <- which(idx4est$prm == loc_else[i])
        lower <- parbd$prm[[loc.tmp]]$lower
        upper <- parbd$prm[[loc.tmp]]$upper

        # parameter estimation
        est <- estimation1(
          r_i = r_i, theta = theta, mod = mod, score.cat = score.cat, D = D,
          nstd = nstd, fix.a.gpcm = ifelse(mod == "GPCM", fix.a.gpcm, FALSE),
          a.val.gpcm = a.val.gpcm, n.1PLM = NULL, aprior = aprior, bprior = bprior,
          use.aprior = use.aprior, use.bprior = use.bprior,
          control = control, startval = startval, lower = lower, upper = upper
        )

        # extract the results
        # item parameter estimates
        a <- ifelse(mod == "GRM", est$pars[1], ifelse(fix.a.gpcm, a.val.gpcm, est$pars[1]))
        if (mod == "GRM") {
          bs <- est$pars[-1]
        } else {
          if (fix.a.gpcm) {
            bs <- est$pars
          } else {
            bs <- est$pars[-1]
          }
        }
        pars <- c(a, bs)
        est_par <- c(est_par, list(pars))

        # standard errors
        a.se <- ifelse(mod == "GRM", est$se[1], ifelse(fix.a.gpcm, NA, est$se[1]))
        if (mod == "GRM") {
          bs.se <- est$se[-1]
        } else {
          if (fix.a.gpcm) {
            bs.se <- est$se
          } else {
            bs.se <- est$se[-1]
          }
        }
        pars.se <- c(a.se, bs.se)
        est_se <- c(est_se, list(pars.se))

        # covariance matrix
        cov_list <- c(cov_list, list(est$covariance))

        # convergence indicator
        convergence <- c(convergence, est$convergence)
        if (est$convergence > 0L) noconv_items <- c(noconv_items, loc_else[i])

        # negative log-likelihood value without the log prior densities
        objective <- c(objective, est$objective.data)
      }
    }
  }

  ## ---------------------------------------------------------------
  # check the convergence of parameter estimation
  if (sum(convergence) == 0L) {
    note <- "All item parameters were successfully converged."
  } else {
    # report the positions in the original data when items with no responses were removed
    if (length(loc_allmiss) > 0L) noconv_items <- loc_nomiss[noconv_items]
    note <- paste0(paste0("Item ", sort(noconv_items), collapse = ", "), " was(were) not successfully converged.")
  }

  ## ---------------------------------------------------------------
  # compute the sum of the log-likelihood values without the log prior densities
  llike <- -sum(objective)

  ## ---------------------------------------------------------------
  # bind per-item estimates and permute rows from estimation order
  # [loc_1p_const, loc_else] back to natural item order
  ord    <- order(c(loc_1p_const, loc_else))
  par_df <- data.frame(bind.fill(est_par, type = "rbind")[ord, , drop = FALSE])
  se_df  <- data.frame(bind.fill(est_se,  type = "rbind")[ord, , drop = FALSE])

  # combine all covariance matrices
  cov_mat <- as.matrix(Matrix::bdiag(cov_list))

  # relocate the matrix elements according to the parameter location
  pos_df <- param_loc$loc.par
  reloc.par <- param_loc$reloc.par
  cov_mat2 <-
    Matrix::sparseMatrix(
      i = c(row(cov_mat)[reloc.par, ]), j = c(col(cov_mat)[, reloc.par]), x = c(cov_mat)
    ) %>%
    as.matrix()

  # arrange the estimated item parameters and standard errors
  # when there is any item which have all missing response data
  if (length(loc_allmiss) > 0L) {
    tmp <-
      x_pre[, -c(1:3)] %>%
      dplyr::mutate_all(.funs = function(x) NA) %>%
      purrr::map(.x = 1:nrow(x_pre), .f = function(i, .data) as.numeric(.data[i, ]))
    par_tmp <- se_tmp <- pos_tmp <- tmp
    par_tmp[loc_nomiss] <-
      par_df %>%
      purrr::map(.x = 1:nrow(par_df), .f = function(i, .data) as.numeric(.data[i, ]))
    se_tmp[loc_nomiss] <-
      se_df %>%
      purrr::map(.x = 1:nrow(se_df), .f = function(i, .data) as.numeric(.data[i, ]))
    pos_tmp[loc_nomiss] <-
      pos_df %>%
      purrr::map(.x = 1:nrow(pos_df), .f = function(i, .data) as.numeric(.data[i, ]))
    par_df <- data.frame(bind.fill(par_tmp, type = "rbind"))
    se_df <- data.frame(bind.fill(se_tmp, type = "rbind"))
    pos_df <- data.frame(bind.fill(pos_tmp, type = "rbind"))
  }

  # create a full data.frame for the item parameter estimates
  if (length(loc_allmiss) > 0L) {
    x <- x_pre
  }
  full_par_df <-
    data.frame(x[, 1:3], par_df) %>%
    confirm_df(g2na = TRUE)

  # create a full data.frame for the standard error estimates
  full_se_df <-
    data.frame(x[, 1:3], se_df) %>%
    confirm_df(g2na = TRUE)

  # create a full data.frame containing the position of item parameter estimates
  # this data.frame is useful when interpreting the variance-covariance matrix of item parameter estimates
  full_pos_df <-
    data.frame(full_par_df[, 1:3], pos_df) %>%
    confirm_df(g2na = TRUE)

  # create a full data.frame including both the item parameter estimates and standard error estimates
  all_df <- data.frame(matrix(NA, nrow = nrow(x), ncol = 2 * ncol(par_df)))
  all_df[, seq(1, 2 * ncol(par_df), 2)] <- par_df
  all_df[, seq(2, 2 * ncol(par_df), 2)] <- se_df
  col.names <- rep(NA, 2 * ncol(par_df))
  col.names[seq(1, 2 * ncol(par_df), 2)] <- paste0("par.", 1:ncol(par_df))
  col.names[seq(2, 2 * ncol(par_df), 2)] <- paste0("se.", 1:ncol(se_df))
  colnames(all_df) <- col.names
  full_all_df <- data.frame(full_par_df[, 1:3], all_df)

  ## ---------------------------------------------------------------
  # check end time
  end.time <- Sys.time()

  # record total computation time
  est_time <- round(as.numeric(difftime(end.time, start.time, units = "secs")), 2)

  # return results
  rst <- structure(
    list(
      estimates = full_all_df, par.est = full_par_df, se.est = full_se_df, pos.par = full_pos_df, covariance = cov_mat2,
      loglikelihood = llike, group.par = group.par, data = data, score = score, scale.D = D, convergence = note, nitem = nitem,
      deleted.item = as.numeric(loc_allmiss), npar.est = length(param_loc$reloc.par),
      n.response = as.numeric(n.resp), TotalTime = est_time
    ),
    class = "est_item"
  )
  rst$call <- cl

  if (verbose) {
    cat("Estimation is finished.", "\n")
  }
  return(rst)
}
