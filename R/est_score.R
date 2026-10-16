#' Estimate examinees' ability (proficiency) parameters
#'
#' This function estimates examinees' latent ability parameters. Available
#' scoring methods include maximum likelihood estimation (ML), maximum
#' likelihood estimation with fences (MLF; Han, 2016), weighted likelihood
#' estimation (WL; Warm, 1989), maximum a posteriori estimation (MAP; Hambleton
#' et al., 1991), expected a posteriori estimation (EAP; Bock & Mislevy, 1982),
#' EAP summed scoring (Thissen et al., 1995; Thissen & Orlando, 2001), and
#' inverse test characteristic curve (TCC) scoring (e.g., Kolen & Brennan, 2004;
#' Kolen & Tong, 2010; Stocking, 1996).
#'
#' @inheritParams est_irt
#' @param x A data frame containing item metadata (e.g., item parameters, number
#'   of categories, IRT model types, etc.), or an object of class `est_irt` or
#'   `est_item` obtained from [irtQ::est_irt()] or [irtQ::est_item()].
#'
#'   See [irtQ::est_irt()] or [irtQ::simdat()] for more details about the item
#'   metadata. This data frame can be easily created using the
#'   [irtQ::shape_df()] function.
#' @param method A character string indicating the scoring method to use.
#'   Available options are:
#'   - `"ML"`: Maximum likelihood estimation
#'   - `"MLF"`: Maximum likelihood estimation with fences (Han, 2016)
#'   - `"WL"`: Weighted likelihood estimation (Warm, 1989)
#'   - `"MAP"`: Maximum a posteriori estimation (Hambleton et al., 1991)
#'   - `"EAP"`: Expected a posteriori estimation (Bock & Mislevy, 1982)
#'   - `"EAP.SUM"`: Expected a posteriori summed scoring (Thissen et al., 1995;
#'   Thissen & Orlando, 2001)
#'   - `"INV.TCC"`: Inverse test characteristic curve scoring
#'   (e.g., Kolen & Brennan, 2004; Kolen & Tong, 2010; Stocking, 1996)
#'
#'   Default is `"ML"`.
#' @param range A numeric vector of length two specifying the lower and upper
#'   bounds of the ability scale for the `"ML"`, `"MLF"`, `"WL"`, and `"MAP"`
#'   methods. Estimates outside these bounds are set to the nearest bound. The
#'   bounds also define the search grid for `stval.opt = 1` and the default
#'   fence locations for `"MLF"`. Default is `c(-5, 5)`.
#' @param norm.prior A numeric vector of length two specifying the mean and
#'   standard deviation of the normal prior distribution. For `"MAP"`, it
#'   defines the prior. For `"EAP"` and `"EAP.SUM"`, it is used to generate the
#'   Gaussian quadrature points and weights when `weights = NULL`. Ignored for
#'   `"ML"`, `"MLF"`, `"WL"`, and `"INV.TCC"`. Default is `c(0, 1)`.
#' @param nquad An integer indicating the number of Gaussian quadrature points
#'   to be generated from the normal prior distribution. Used only when `method`
#'   is `"EAP"` or `"EAP.SUM"` and `weights = NULL`. Default is `41`.
#' @param weights A two-column matrix or data frame containing the quadrature
#'   points (in the first column) and their corresponding weights (in the second
#'   column) for the latent variable prior distribution. The weights and points
#'   can be conveniently generated using the function [irtQ::gen.weight()].
#'
#'   If `NULL` and `method` is either `"EAP"` or `"EAP.SUM"`, default quadrature
#'   values are generated based on the `norm.prior` and `nquad` arguments.
#'   Ignored if `method` is `"ML"`, `"MLF"`, `"WL"`, `"MAP"`, or `"INV.TCC"`.
#' @param fence.a A numeric value specifying the slope (*a*) parameter of the
#'   two fence items used in MLF. See **Details** below. Default is 3.0.
#' @param fence.b A numeric vector of length two specifying the difficulty (*b*)
#'   parameters of the lower and upper fence items used in MLF. If `NULL`, the
#'   values in `range` are used. Default is `NULL`.
#' @param tol A positive numeric value specifying the convergence tolerance. For
#'   `"ML"`, `"MLF"`, `"WL"`, and `"MAP"`, the Fisher scoring iterations stop
#'   when the absolute update is less than `tol`, and `tol` is also the
#'   tolerance of the root search used when the iterations do not converge. For
#'   `"INV.TCC"`, the bisection search stops when the search interval is no
#'   wider than `tol`. Default is 1e-4.
#' @param max.iter A positive integer specifying the maximum number of Fisher
#'   scoring iterations for `"ML"`, `"MLF"`, `"WL"`, and `"MAP"`. Default is
#'   100. See **Details** for the case where the iterations do not converge.
#' @param stval.opt A positive integer specifying the starting value option for
#'   the ML, MLF, WL, and MAP scoring methods. Available options are:
#'   - 1: Grid search (default)
#'   - 2: Based on observed sum scores
#'   - 3: Fixed at 0
#'
#'   See **Details** below for more information.
#' @param se Logical. If `TRUE`, standard errors of ability estimates are
#'   computed. If `FALSE`, the standard errors are returned as `NA`. For
#'   `"EAP.SUM"` and `"INV.TCC"`, standard errors are always computed. Default
#'   is `TRUE`.
#' @param intpol Logical. If `TRUE` and `method = "INV.TCC"`, ability estimates
#'   are assigned to sum scores that cannot be mapped through the TCC (sum
#'   scores less than or equal to the sum of the guessing parameters, and the
#'   maximum possible sum score). Default is `TRUE`. See **Details** below.
#' @param range.tcc A numeric vector of length two giving the ability estimates
#'   assigned to the lowest and the maximum possible sum scores when
#'   `method = "INV.TCC"` and `intpol = TRUE`. Default is `c(-7, 7)`.
#' @param missing A value indicating missing responses in the data set. Default
#'   is `NA`. See **Details** below.
#' @param ncore An integer specifying the number of logical CPU cores used for
#'   parallel scoring with `"ML"`, `"MLF"`, `"WL"`, `"MAP"`, and `"EAP"`.
#'   Default is 1. See **Details** below.
#' @param ... Additional arguments passed to [parallel::makeCluster()] when
#'   `ncore > 1`.
#'
#' @details For `"MAP"`, the prior is the normal distribution given by
#'   `norm.prior`.
#'
#'   Each response must be an integer from 0 to the number of categories of the
#'   item minus 1, and `data` must have one column per item in `x`. Missing
#'   responses must be coded as `NA` or declared with the `missing`
#'   argument. For `"ML"`, `"MLF"`, `"WL"`, `"MAP"`, and `"EAP"`, missing
#'   responses are excluded from the likelihood, and examinees with all
#'   responses missing receive `NA` with a warning. For `"EAP.SUM"` and
#'   `"INV.TCC"`, missing responses are recoded as 0 with a warning.
#'
#'   In maximum likelihood estimation with fences (MLF; Han, 2016), two
#'   imaginary 2PL items are added to every response pattern: a lower fence
#'   item answered correctly and an upper fence item answered incorrectly. The
#'   lower fence should have a *b*-parameter below, and the upper fence a
#'   *b*-parameter above, all item difficulties in the test. Both should have
#'   steep slopes (large *a*-parameters). If `fence.b = NULL`, the fences are
#'   placed at the bounds given in `range`. See Han (2016) for details.
#'
#'   For `"INV.TCC"`, the ability estimate for a sum score X is the root of
#'   TCC(theta) = X, found by the bisection method (Howard, 2017). No root
#'   exists for the maximum possible sum score or for sum scores less than or
#'   equal to the sum of the guessing parameters, which is only the score 0 when
#'   there are no guessing parameters. With
#'   `intpol = TRUE`, these scores receive estimates as follows. Let
#'   \eqn{\theta_{min}} and \eqn{\theta_{max}} be the first and second values
#'   of `range.tcc`, and let \eqn{\theta_{X}} be the estimate for the smallest
#'   sum score X that is greater than the sum of the guessing parameters.
#'   Scores below X are mapped onto the line through
#'   \eqn{(x = \theta_{min}, y = 0)} and \eqn{(x = \theta_{X}, y = X)}, and the
#'   maximum possible sum score receives \eqn{\theta_{max}}. If
#'   \eqn{\theta_{min}} is above \eqn{\theta_{X}} or \eqn{\theta_{max}} is
#'   below the largest root, a warning is issued and no values are
#'   assigned. With `intpol = FALSE`, these scores receive `NA`.
#'
#'   For `"INV.TCC"`, the standard error for a sum score is the standard
#'   deviation of the inverse TCC estimates over the conditional sum score
#'   distribution at the estimate for that score, computed with the
#'   Lord-Wingersky recursion (Lim et al., 2021; see [irtQ::lwrc()]). When some
#'   sum scores have no estimate (`intpol = FALSE`, or the warning described
#'   above), the distribution is restricted to the sum scores with an estimate
#'   and rescaled to sum to one. The implementation is based on a modified
#'   version of `SNSequate::irt.eq.tse()` from the \pkg{SNSequate} package
#'   (Gonzalez, 2014).
#'
#'   For the ML, MLF, WL, and MAP scoring methods, different strategies can be
#'   used to determine the starting value for ability estimation based on the
#'   `stval.opt` argument:
#'
#'   - When `stval.opt = 1` (default), the log-likelihood (log-posterior for
#'   `"MAP"`) is evaluated on a grid from `range[1]` to `range[2]` in steps of
#'   0.1. The grid point at the highest local maximum is the starting value. If
#'   there is no local maximum on the grid, the starting value is 0.
#'
#'   - When `stval.opt = 2`, the starting value is the log-odds of the
#'   observed sum score `obs.score` relative to the maximum possible score
#'   `max.score`, both computed over the non-missing items:
#'   `log(obs.score / (max.score - obs.score))`.
#'     - If `obs.score = 0`, the starting value is `log(1 / max.score)`.
#'     - If `obs.score = max.score`, the starting value is `log(max.score)`.
#'
#'   - When `stval.opt = 3`, the starting value is fixed at 0.
#'
#'   For the ML, MLF, WL, and MAP scoring methods, the estimate is found by
#'   Fisher scoring. This is a Newton-type iteration that divides the first
#'   derivative of the log-likelihood (of the log-posterior for `"MAP"`, or
#'   Warm's weighted score function for `"WL"`) by the expected Fisher
#'   information instead of the observed second derivative, and a single update
#'   is capped at 1 in absolute value. The iterations stop when the absolute
#'   update is less than `tol`. If they have not converged after `max.iter`
#'   iterations (for example, when the iterations alternate between two values
#'   for a low-scoring examinee on a test with 3PLM items), the estimate is
#'   found as follows. Among the values visited during the iterations, the
#'   function takes the closest pair in which the score function (the
#'   derivative above) is positive at the lower value and negative at the
#'   upper value, so that the pair encloses a local maximum, and finds the root
#'   of the score function between them with [stats::uniroot()] using the
#'   tolerance `tol`. If there is no such pair, the last value of the
#'   iterations is used. This search is not carried out when the iterations
#'   converge. In all cases, the estimate is truncated to `range`, and the
#'   standard error is computed at the final estimate.
#'
#'   For `"ML"`, `"MLF"`, `"WL"`, `"MAP"`, and `"EAP"`, examinees can be scored
#'   in parallel by setting `ncore` greater than 1. A warning is issued when
#'   `ncore > 1` and there are fewer than 5,000 examinees, because the parallel
#'   overhead then exceeds the computation time.
#'
#'   For `"ML"`, `"MLF"`, `"WL"`, and `"MAP"`, the standard error is the
#'   inverse square root of the expected Fisher information at the estimate.
#'   The information includes the two fence items for `"MLF"` and the prior
#'   information \eqn{1/\sigma^2} for `"MAP"`. For `"EAP"` and `"EAP.SUM"`,
#'   the standard error is the posterior standard deviation. For `"EAP.SUM"`,
#'   the posterior for each sum score is computed with the Lord-Wingersky
#'   recursion (see [irtQ::lwrc()]).
#'
#'   For the implementation of the WL method, the function references the
#'   `catR::Pi()`, `catR::Ji()`, and `catR::Ii()` functions from the \pkg{catR}
#'   package (Magis & Barrada, 2017).
#'
#' @return For `method = "ML"`, `"MLF"`, `"WL"`, `"MAP"`, or `"EAP"`, a data
#'   frame with one row per examinee and two columns:
#'   - `est.theta`: Ability estimates.
#'   - `se.theta`: Standard errors of the ability estimates (`NA` when
#'     `se = FALSE`).
#'
#'   For `"ML"`, `"MLF"`, `"WL"`, and `"MAP"`, the standard error is set to
#'   99.9999 when the estimate equals a bound of `range`. Examinees with all
#'   responses missing receive `NA` in both columns.
#'
#'   For `method = "EAP.SUM"` or `"INV.TCC"`, a list with two data frames:
#'   - `est.par`: One row per examinee with the columns `sum.score`,
#'     `est.theta`, and `se.theta`.
#'   - `score.table`: One row per possible sum score, from 0 to the maximum
#'     possible score, with the same three columns.
#'
#'   For `"INV.TCC"`, sum scores without an estimate (see **Details**) have
#'   `NA` in `est.theta` and `se.theta`.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [irtQ::est_irt()], [irtQ::simdat()], [irtQ::shape_df()],
#'   [irtQ::gen.weight()]
#'
#' @references Bock, R. D., & Mislevy, R. J. (1982). Adaptive EAP estimation of
#'   ability in a microcomputer environment. *Applied Psychological Measurement,
#'   6*(4), 431-444. \doi{10.1177/014662168200600405}.
#'
#'   Gonzalez, J. (2014). SNSequate: Standard and nonstandard statistical models
#'   and methods for test equating. *Journal of Statistical Software, 59*(7),
#'   1-30. \doi{10.18637/jss.v059.i07}.
#'
#'   Hambleton, R. K., Swaminathan, H., & Rogers, H. J. (1991). *Fundamentals of
#'   item response theory*. Newbury Park, CA: Sage.
#'
#'   Han, K. T. (2016). Maximum likelihood score estimation method with fences
#'   for short-length tests and computerized adaptive tests. *Applied
#'   Psychological Measurement, 40*(4), 289-301. \doi{10.1177/0146621616631317}.
#'
#'   Howard, J. P. (2017). *Computational methods for numerical analysis with
#'   R*. New York: Chapman and Hall/CRC.
#'
#'   Kolen, M. J., & Brennan, R. L. (2004). *Test equating, scaling, and
#'   linking* (2nd ed.). Springer.
#'
#'   Kolen, M. J., & Tong, Y. (2010). Psychometric properties of IRT proficiency
#'   estimates. *Educational Measurement: Issues and Practice, 29*(3), 8-14.
#'   \doi{10.1111/j.1745-3992.2010.00179.x}.
#'
#'   Lim, H., Davey, T., & Wells, C. S. (2021). A recursion-based analytical
#'   approach to evaluate the performance of MST. *Journal of Educational
#'   Measurement, 58*(2), 154-178. \doi{10.1111/jedm.12276}.
#'
#'   Magis, D., & Barrada, J. R. (2017). Computerized adaptive testing with R:
#'   Recent updates of the package catR. *Journal of Statistical Software,
#'   76*(Code Snippet 1), 1-19. \doi{10.18637/jss.v076.c01}.
#'
#'   Stocking, M. L. (1996). An alternative method for scoring adaptive tests.
#'   *Journal of Educational and Behavioral Statistics, 21*(4), 365-389.
#'   \doi{10.3102/10769986021004365}.
#'
#'   Thissen, D., & Orlando, M. (2001). Item response theory for items scored in
#'   two categories. In D. Thissen & H. Wainer (Eds.), *Test scoring* (pp.
#'   73-140). Mahwah, NJ: Lawrence Erlbaum.
#'
#'   Thissen, D., Pommerich, M., Billeaud, K., & Williams, V. S. L. (1995). Item
#'   response theory for scores on tests including polytomous items with ordered
#'   responses. *Applied Psychological Measurement, 19*(1), 39-49.
#'   \doi{10.1177/014662169501900105}.
#'
#'   Warm, T. A. (1989). Weighted likelihood estimation of ability in item
#'   response theory. *Psychometrika, 54*(3), 427-450. \doi{10.1007/BF02294627}.
#'
#' @examples
#' ## Import the "-prm.txt" output file from flexMIRT
#' flex_prm <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
#'
#' # Read item parameters and convert them into item metadata
#' x <- bring.flexmirt(file = flex_prm, "par")$Group1$full_df
#'
#' # Generate examinee ability values
#' set.seed(12)
#' theta <- rnorm(10)
#'
#' # Simulate item response data based on the item metadata and abilities
#' data <- simdat(x, theta, D = 1)
#'
#' \donttest{
#' # Estimate abilities using maximum likelihood (ML)
#' est_score(x, data, D = 1, method = "ML", range = c(-4, 4), se = TRUE)
#'
#' # Estimate abilities using weighted likelihood (WL)
#' est_score(x, data, D = 1, method = "WL", range = c(-4, 4), se = TRUE)
#'
#' # Estimate abilities using MLF with default fences
#' # based on the `range` argument
#' est_score(x, data,
#'   D = 1, method = "MLF",
#'   fence.a = 3.0, fence.b = NULL, se = TRUE
#' )
#'
#' # Estimate abilities using MLF with user-specified fences
#' est_score(x, data,
#'   D = 1, method = "MLF", fence.a = 3.0,
#'   fence.b = c(-7, 7), se = TRUE
#' )
#'
#' # Estimate abilities using maximum a posteriori (MAP)
#' est_score(x, data,
#'   D = 1, method = "MAP", norm.prior = c(0, 1),
#'   se = TRUE
#' )
#'
#' # Estimate abilities using expected a posteriori (EAP)
#' est_score(x, data,
#'   D = 1, method = "EAP", norm.prior = c(0, 1),
#'   nquad = 30, se = TRUE
#' )
#'
#' # Estimate abilities using EAP summed scoring
#' est_score(x, data,
#'   D = 1, method = "EAP.SUM", norm.prior = c(0, 1),
#'   nquad = 30
#' )
#'
#' # Estimate abilities using inverse TCC scoring
#' est_score(x, data,
#'   D = 1, method = "INV.TCC", intpol = TRUE,
#'   range.tcc = c(-7, 7)
#' )
#' }
#'
#' @export
est_score <- function(x, ...) UseMethod("est_score")

#' @describeIn est_score Default method to estimate examinees' latent ability
#'  parameters using a data frame `x` containing the item metadata.
#' @importFrom dplyr bind_rows
#' @export
est_score.default <- function(x,
                              data,
                              D = 1,
                              method = "ML",
                              range = c(-5, 5),
                              norm.prior = c(0, 1),
                              nquad = 41,
                              weights = NULL,
                              fence.a = 3.0,
                              fence.b = NULL,
                              tol = 1e-4,
                              max.iter = 100,
                              se = TRUE,
                              stval.opt = 1,
                              intpol = TRUE,
                              range.tcc = c(-7, 7),
                              missing = NA,
                              ncore = 1,
                              ...) {
  # stop when the scoring method is not one of the available options
  if (!(is.character(method) && length(method) == 1L &&
    method %in% c("ML", "MLF", "WL", "MAP", "EAP", "EAP.SUM", "INV.TCC"))) {
    stop(
      "'method' must be one of \"ML\", \"MLF\", \"WL\", \"MAP\", \"EAP\", ",
      "\"EAP.SUM\", or \"INV.TCC\".",
      call. = FALSE
    )
  }

  # stop when the starting value option is not 1, 2, or 3
  if (!(length(stval.opt) == 1L && stval.opt %in% 1:3)) {
    stop("'stval.opt' must be 1, 2, or 3.", call. = FALSE)
  }

  # stop when the convergence tolerance is not a single positive number
  if (method %in% c("ML", "MLF", "WL", "MAP", "INV.TCC") &&
    !(is.numeric(tol) && length(tol) == 1L && is.finite(tol) && tol > 0)) {
    stop("'tol' must be a single positive number.", call. = FALSE)
  }

  # stop when the maximum number of iterations is not a positive whole number
  if (method %in% c("ML", "MLF", "WL", "MAP") &&
    !(is.numeric(max.iter) && length(max.iter) == 1L && is.finite(max.iter) &&
      max.iter >= 1 && max.iter == round(max.iter))) {
    stop("'max.iter' must be a single positive whole number.", call. = FALSE)
  }

  # convert a single examinee's response vector to a one-row matrix
  if (is.vector(data)) {
    data <- rbind(data)
  }

  # re-code missing values
  if (!is.na(missing)) {
    data[data == missing] <- NA
  }

  # check the number of examinees
  nstd <- nrow(data)

  # confirm and correct all item metadata information
  x <- confirm_df(x)

  # stop when the response data do not have one column per item
  if (ncol(data) != nrow(x)) {
    stop(
      "The number of columns in 'data' (", ncol(data), ") must equal ",
      "the number of items in 'x' (", nrow(x), ").",
      call. = FALSE
    )
  }

  # convert the responses to a matrix for validation
  resp_chk <- as.matrix(data)

  # read character or factor responses as numbers, as the scoring does
  if (!is.numeric(resp_chk) && !is.logical(resp_chk)) {
    resp_num <- suppressWarnings(as.numeric(resp_chk))

    # stop when an observed response cannot be read as a number
    if (any(is.na(resp_num) & !is.na(resp_chk))) {
      stop(
        "Responses must be numeric scores; ",
        "check the 'missing' argument for missing-value codes.",
        call. = FALSE
      )
    }

    # keep the numeric responses with the dimensions and names of the input
    resp_chk <- matrix(resp_num, nrow = nrow(resp_chk), dimnames = dimnames(resp_chk))

    # score the numeric responses
    data <- resp_chk
  }

  # flag observed responses that are not integer scores within each item's categories
  bad_resp <- !is.na(resp_chk) &
    (resp_chk != round(resp_chk) | resp_chk < 0 | t(t(resp_chk) > (x$cats - 1)))

  # stop when any observed response is invalid
  if (any(bad_resp)) {
    stop(
      "Responses must be integers from 0 to (number of categories - 1) ",
      "of each item; check the 'missing' argument for missing-value codes.",
      call. = FALSE
    )
  }

  # scoring of ML, WL, MLF, MAP, and EAP
  if (method %in% c("ML", "MAP", "WL", "EAP", "MLF")) {
    # check if the method is WL
    # if TRUE, ji = TRUE
    ji <- ifelse(method == "WL", TRUE, FALSE)

    # add two fence items and their responses when MLF is used
    if (method == "MLF") {
      # when fence.b = NULL, use the range argument as the fence.b argument
      if (is.null(fence.b)) {
        fence.b <- range
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

    # pre-populate elm_item with pars, model, cats, id from item metadata
    elm_item <- breakdown(x)

    # pre-compute quadrature points once for EAP method (not per examinee)
    popdist <- if (method == "EAP") {
      if (is.null(weights)) {
        gen.weight(n = nquad, dist = "norm", mu = norm.prior[1], sigma = norm.prior[2])
      } else {
        data.frame(weights)
      }
    } else NULL

    # check the number of CPU cores
    if (ncore < 1) {
      stop("The number of logical CPU cores must not be less than 1.", call. = FALSE)
    }

    # warn when parallel overhead likely exceeds computation gain
    if (ncore > 1 && nstd < 5000) {
      warning(
        "ncore > 1 is not recommended for N < 5,000 ",
        "as parallel overhead exceeds computation time. ",
        "Consider using ncore = 1.",
        call. = FALSE
      )
    }

    # estimation
    if (ncore == 1L) {
      # score all examinees in the current process
      rst <- est_score_1core(
        elm_item = elm_item, data = data, D = D,
        method = method, max.cats = max.cats,
        range = range, norm.prior = norm.prior,
        tol = tol, max.iter = max.iter,
        se = se, stval.opt = stval.opt, ji = ji,
        popdist = popdist        # pre-computed quadrature (NULL for non-EAP)
      )
    } else {
      # create a parallel processing cluster
      cl <- parallel::makeCluster(ncore, ...)

      # stop the cluster even when a worker fails
      on.exit(parallel::stopCluster(cl), add = TRUE)

      # split the row indices into at most ncore nearly equal chunks
      chunk_idx <- parallel::splitIndices(nstd, ncore)

      # keep every chunk as a two-dimensional object even when it holds one row
      data_list <- lapply(chunk_idx, function(ix) data[ix, , drop = FALSE])

      # delete 'data' object
      rm(data, envir = environment(), inherits = FALSE)

      # export pre-populated elm_item, pre-computed popdist, and required functions.
      parallel::clusterExport(cl, c(
        "elm_item", "popdist", "D", "method",
        "max.cats", "range", "norm.prior",
        "tol", "max.iter", "se", "stval.opt", "ji",
        "est_score_1core", "est_score_indiv", "idxfinder",
        "ll_score", "drm", "prm", "gpcm", "grm",
        "logprior_deriv", "esprior_norm",
        "info_score", "info_drm", "info_prm",
        "gen.weight"
      ), envir = environment())
      # pre-load Rfast on workers (used in ll_score, info_drm, etc.)
      parallel::clusterEvalQ(cl, library(Rfast))

      # set a function for scoring
      fsm <- function(subdat) {
        est_score_1core(
          elm_item = elm_item, data = subdat, D = D,
          method = method, max.cats = max.cats,
          range = range, norm.prior = norm.prior,
          tol = tol, max.iter = max.iter,
          se = se, stval.opt = stval.opt, ji = ji,
          popdist = popdist        # pre-computed quadrature passed to workers
        )
      }

      # parallel scoring
      est <- parallel::parLapply(cl = cl, X = data_list, fun = fsm)

      # combine the results
      rst <- do.call(what = "rbind", args = est)
    }

    # return a warning message when some examinees have all missing responses
    loc.na <- which(is.na(rst$est.theta))
    if (length(loc.na) > 0) {
      memo <- "NA values are returned for examinees with all missing responses."
      warning(memo, call. = FALSE)
    }
  }

  if (method == "EAP.SUM") {
    rst <- eap_sum(
      x = x, data = data, norm.prior = norm.prior,
      nquad = nquad, weights = weights, D = D
    )
  }

  if (method == "INV.TCC") {
    rst <- inv_tcc(x, data,
      D = D, intpol = intpol, range.tcc = range.tcc,
      tol = tol, max.it = 500
    )
  }

  # return results
  rst
}


#' @describeIn est_score Method for an object of class `est_irt`. The item
#'  parameter estimates, response data, and scaling constant `D` are taken from
#'  `x`, so supplying `data` or `D` stops with an error.
#' @export
est_score.est_irt <- function(x,
                              method = "ML",
                              range = c(-5, 5),
                              norm.prior = c(0, 1),
                              nquad = 41,
                              weights = NULL,
                              fence.a = 3.0,
                              fence.b = NULL,
                              tol = 1e-4,
                              max.iter = 100,
                              se = TRUE,
                              stval.opt = 1,
                              intpol = TRUE,
                              range.tcc = c(-7, 7),
                              missing = NA,
                              ncore = 1,
                              ...) {
  # stop when 'data' or 'D' is supplied, since both are taken from the fitted object
  if (any(c("data", "D") %in% names(list(...)))) {
    stop(
      "'data' and 'D' are taken from the est_irt object and cannot be supplied.",
      call. = FALSE
    )
  }

  # score the stored response data with the stored item estimates and scaling constant
  est_score.default(
    x = x$par.est, data = x$data, D = x$scale.D, method = method,
    range = range, norm.prior = norm.prior, nquad = nquad,
    weights = weights, fence.a = fence.a, fence.b = fence.b,
    tol = tol, max.iter = max.iter, se = se, stval.opt = stval.opt,
    intpol = intpol, range.tcc = range.tcc, missing = missing,
    ncore = ncore, ...
  )
}


#' @describeIn est_score Method for an object of class `est_item`. The item
#'  parameter estimates, response data, and scaling constant `D` are taken from
#'  `x`, so supplying `data` or `D` stops with an error.
#' @export
est_score.est_item <- function(x,
                               method = "ML",
                               range = c(-5, 5),
                               norm.prior = c(0, 1),
                               nquad = 41,
                               weights = NULL,
                               fence.a = 3.0,
                               fence.b = NULL,
                               tol = 1e-4,
                               max.iter = 100,
                               se = TRUE,
                               stval.opt = 1,
                               intpol = TRUE,
                               range.tcc = c(-7, 7),
                               missing = NA,
                               ncore = 1,
                               ...) {
  # stop when 'data' or 'D' is supplied, since both are taken from the fitted object
  if (any(c("data", "D") %in% names(list(...)))) {
    stop(
      "'data' and 'D' are taken from the est_item object and cannot be supplied.",
      call. = FALSE
    )
  }

  # drop the items that were excluded from the calibration for having no responses
  par_est <- x$par.est
  if (length(x$deleted.item) > 0L) {
    par_est <- par_est[-x$deleted.item, , drop = FALSE]
  }

  # score the stored response data with the stored item estimates and scaling constant
  est_score.default(
    x = par_est, data = x$data, D = x$scale.D, method = method,
    range = range, norm.prior = norm.prior, nquad = nquad,
    weights = weights, fence.a = fence.a, fence.b = fence.b,
    tol = tol, max.iter = max.iter, se = se, stval.opt = stval.opt,
    intpol = intpol, range.tcc = range.tcc, missing = missing,
    ncore = ncore, ...
  )
}


# Score a block of examinees; used for all examinees when ncore = 1 and for
# each chunk on a worker when ncore > 1
est_score_1core <- function(elm_item,
                            data,
                            D = 1,
                            method = "ML",
                            max.cats,
                            range = c(-4, 4),
                            norm.prior = c(0, 1),
                            tol = 1e-4,
                            max.iter = 30,
                            se = TRUE,
                            stval.opt = 1,
                            ji = FALSE,
                            popdist = NULL) {
  # check the number of examinees in this chunk
  nstd <- nrow(data)

  # classify DRM/PRM items once for this chunk (called once per worker, not per examinee)
  idx_full     <- idxfinder(elm_item)
  idx_drm_full <- idx_full$idx.drm
  idx_prm_full <- idx_full$idx.prm

  # score each examinee directly from the raw response row
  est <- lapply(seq_len(nstd), function(i) {
    resp_vec_i <- data[i, ]

    # subset to non-NA items only
    na_mask <- !is.na(resp_vec_i)

    # return NA for examinees with all missing responses
    if (!any(na_mask)) {
      return(data.frame(est.theta = NA_real_, se.theta = NA_real_))
    }

    # subset elm_item to observed items
    elm_sub       <- elm_item
    elm_sub$pars  <- elm_item$pars[na_mask, , drop = FALSE]
    elm_sub$model <- elm_item$model[na_mask]
    elm_sub$cats  <- elm_item$cats[na_mask]

    resp_sub  <- as.numeric(resp_vec_i[na_mask])
    obs_sum_i <- if (stval.opt == 2L) sum(resp_sub) else NULL

    # map pre-computed full-item indices to the non-NA subset
    if (all(na_mask)) {
      idx_drm_i <- idx_drm_full
      idx_prm_i <- idx_prm_full
    } else {
      na_pos    <- which(na_mask)
      idx_drm_i <- if (!is.null(idx_drm_full)) {
        loc <- which(na_pos %in% idx_drm_full); if (length(loc) == 0L) NULL else loc
      } else NULL
      idx_prm_i <- if (!is.null(idx_prm_full)) {
        loc <- which(na_pos %in% idx_prm_full); if (length(loc) == 0L) NULL else loc
      } else NULL
    }

    est_score_indiv(
      resp_vec = resp_sub,
      elm_item = elm_sub,
      max.cats = max.cats,
      idx.drm  = idx_drm_i,
      idx.prm  = idx_prm_i,
      D = D, method = method,
      range = range, norm.prior = norm.prior,
      tol = tol, max.iter = max.iter, se = se,
      stval.opt = stval.opt, ji = ji,
      obs.sum = obs_sum_i,
      popdist = popdist        # pre-computed quadrature (NULL for non-EAP)
    )
  })

  # combine per-examinee results and return
  do.call(rbind, est)
}


# This function computes an ability estimate for a single examinee (ML, WL, MLF, MAP, EAP)
est_score_indiv <- function(resp_vec, elm_item, max.cats, idx.drm, idx.prm,
                            D = 1, method = "ML",
                            range = c(-4, 4), norm.prior = c(0, 1),
                            tol = 1e-4, max.iter = 30, se = TRUE,
                            stval.opt = 1, ji = FALSE, obs.sum = NULL,
                            popdist = NULL) {
  # elm_item is pre-populated (pars, model, cats) for the observed (non-NA) items only;
  # idx.drm/idx.prm are pre-computed by the caller, outside the loop
  n.resp <- nrow(elm_item$pars)

  # build the n.resp x max.cats one-hot freq.cat via direct matrix indexing
  freq.cat <- matrix(0L, nrow = n.resp, ncol = max.cats)
  resp_int  <- as.integer(resp_vec)          # 0-based integer responses (no NAs: caller subsets)
  freq.cat[cbind(seq_len(n.resp), resp_int + 1L)] <- 1L

  ## ----------------------------------------------------
  ## ML, WL, MLF, and MAP
  if (method %in% c("ML", "WL", "MLF", "MAP")) {
    # set a starting value
    if (stval.opt == 1) {
      # use a grid search over range to find the starting value
      # prepare the discrete theta values
      theta.nodes <- seq(from = range[1], to = range[2], by = 0.1)

      # compute the negative log-likelihood values for all the discrete theta values
      ll_tmp <- ll_score(
        theta = theta.nodes, elm_item = elm_item, freq.cat = freq.cat,
        method = method, idx.drm = idx.drm, idx.prm = idx.prm, D = D,
        norm.prior = norm.prior, logL = TRUE
      )

      # find the local minima of the negative log-likelihood on the grid
      loc_change <- which(diff(sign(diff(ll_tmp))) > 0L) + 1

      # select the local minimum with the smallest negative log-likelihood
      stval_tmp1 <- theta.nodes[loc_change][which.min(ll_tmp[loc_change])]

      # use 0 when there is no local minimum (the function is monotone on the grid)
      theta <- ifelse(length(stval_tmp1) > 0L, stval_tmp1, 0)
    } else if (stval.opt == 2) {
      # compute the maximum possible sum score of the observed items
      total.nc <- sum(elm_item$cats - 1)

      # obs.sum is passed in by the caller (sum of non-NA responses)
      if (obs.sum == 0) {
        theta <- log(1 / total.nc)
      } else if (obs.sum == total.nc) {
        theta <- log(total.nc / 1)
      } else {
        theta <- log(obs.sum / (total.nc - obs.sum))
      }
    } else if (stval.opt == 3) {
      # use 0 as a starting value
      theta <- 0
    }

    # estimate an ability using Fisher scoring
    # set the iteration number to 0
    i <- 0
    abs_delta <- 1
    # preserve the last protected finfo for SE reuse on clean convergence;
    # initialized to 1e-5 (the floor) in case the loop body never executes
    finfo_last <- 1e-5

    # record the visited theta values and their gradients for the fallback root search
    th_hist <- numeric(0)
    gr_hist <- numeric(0)
    while (abs_delta >= tol) {
      # update the iteration number
      i <- i + 1

      # compute the gradient of the negative log-likelihood and the Fisher
      # information (negative expected second derivative of the log-likelihood)
      gr_fi <-
        info_score(
          theta = theta, elm_item = elm_item, freq.cat = freq.cat,
          idx.drm = idx.drm, idx.prm = idx.prm, method = method, D = D,
          norm.prior = norm.prior, grad = TRUE, ji = ji
        )
      grad <- gr_fi$grad
      finfo <- gr_fi$finfo

      # store the current theta value and its gradient
      th_hist <- c(th_hist, theta)
      gr_hist <- c(gr_hist, grad)

      # floor the Fisher information at 1e-5 to avoid division by values near 0
      finfo[finfo < 1e-5 | is.nan(finfo)] <- 1e-5

      # save protected finfo at current theta before the theta update
      finfo_last <- finfo

      # compute the theta correction factor (delta)
      delta <- grad / finfo

      # cap the step at 1 in absolute value
      abs_delta <- abs(delta)
      delta[abs_delta > 1] <- sign(delta)

      # update the theta value
      theta <- theta - delta

      if (i == max.iter) break
    }
    # flag whether the loop exited via convergence (abs_delta < tol) or
    # hit the iteration ceiling; finfo_last is only safe to reuse when converged
    nr_converged <- (abs_delta < tol)

    # solve for the gradient root inside a visited bracket when the iterations do not converge
    if (!nr_converged) {
      # sort the distinct visited points with a negative gradient
      lo_v <- sort(unique(th_hist[which(gr_hist < 0)]))

      # take the visited points with a positive gradient in the order of the iterations
      hi_v <- th_hist[which(gr_hist > 0)]

      # locate the largest negative-gradient point strictly below each positive-gradient point
      k <- findInterval(hi_v, lo_v, left.open = TRUE)

      # search only when such a bracket exists
      if (any(k > 0L)) {
        # take the narrowest bracket
        ok <- which(k > 0L)
        j <- ok[which.min(hi_v[ok] - lo_v[k[ok]])]
        br <- c(lo_v[k[j]], hi_v[j])

        # define the gradient of the objective function as a function of theta
        f_grad <- function(t) {
          info_score(
            theta = t, elm_item = elm_item, freq.cat = freq.cat,
            idx.drm = idx.drm, idx.prm = idx.prm, method = method, D = D,
            norm.prior = norm.prior, grad = TRUE, ji = ji
          )$grad
        }

        # find the root of the gradient within the bracket
        theta <- stats::uniroot(f_grad, lower = br[1], upper = br[2], tol = tol)$root
      }
    }

    # truncate the estimate to the bounds of range
    theta[theta <= range[1]] <- range[1]
    theta[theta >= range[2]] <- range[2]
    est.theta <- theta

    # compute the standard error
    if (se) {
      if (est.theta %in% range) {
        se.theta <- 99.9999
      } else if (nr_converged) {
        # reuse finfo_last (at theta_prev = est.theta + delta where |delta| < tol);
        # avoids a second info_score() call; approximation error is O(tol)
        se.theta <- 1 / sqrt(finfo_last)
      } else {
        # max.iter reached without convergence: the last delta may be large, so
        # finfo_last could be far from finfo(est.theta); recompute it at est.theta
        finfo_se <-
          info_score(
            theta = est.theta, elm_item = elm_item, freq.cat = freq.cat,
            idx.drm = idx.drm, idx.prm = idx.prm, method = method, D = D,
            norm.prior = norm.prior, grad = FALSE, ji = FALSE
          )$finfo
        se.theta <- 1 / sqrt(finfo_se)
      }
    } else {
      se.theta <- NA
    }
  }

  ## ----------------------------------------------------
  ## EAP scoring
  if (method == "EAP") {
    # popdist is pre-computed by the caller (gen.weight() is not called per examinee)

    # compute the log-likelihood at the quadrature points
    loglik <-
      -ll_score(
        theta = popdist[, 1], elm_item = elm_item, freq.cat = freq.cat,
        idx.drm = idx.drm, idx.prm = idx.prm, D = D, logL = TRUE
      )

    # add the log prior weights to obtain the unnormalized log posterior
    lpost <- loglik + log(popdist[, 2])

    # shift by the maximum so that the exponential does not underflow
    posterior <- exp(lpost - max(lpost))

    # compute the posterior mean (EAP)
    posterior <- posterior / sum(posterior)
    est.theta <- sum(popdist[, 1] * posterior)

    if (se) {
      # compute the posterior standard deviation
      ex2 <- sum(popdist[, 1]^2 * posterior)
      var <- ex2 - (est.theta)^2
      se.theta <- sqrt(var)
    } else {
      se.theta <- NA
    }
  }

  # combine theta and se into a data frame
  rst <- data.frame(est.theta = est.theta, se.theta = se.theta)

  # return results
  rst
}
