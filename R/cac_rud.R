#' Classification Accuracy and Consistency Based on Rudner's (2001, 2005)
#' Approach
#'
#' This function computes classification accuracy and consistency indices with
#' the approach of Rudner (2001, 2005), which assumes that an ability estimate
#' is normally distributed around the true ability with a standard deviation
#' equal to its standard error. The indices can be computed over an ability
#' distribution given by quadrature points and weights, or over individual
#' ability estimates.
#'
#' @param x A data frame containing item metadata (e.g., item parameters, number
#'   of categories, IRT model types, etc.). See [irtQ::est_irt()] or
#'   [irtQ::simdat()] for more details about the item metadata. This data frame
#'   can be easily created using the [irtQ::shape_df()] function.
#'   If `x = NULL`, the `se` argument must be explicitly provided.
#'   Defaults to `NULL`.
#' @param cutscore A numeric vector of finite cut scores on the ability (theta)
#'   metric, in strictly ascending order. The K - 1 cut scores divide examinees
#'   into K performance levels. An ability value equal to a cut score is
#'   assigned to the higher level.
#' @param theta A numeric vector of ability estimates, one per examinee. It is
#'   ignored when `weights` is supplied. An examinee with a missing `theta` (or
#'   a missing `se`) is excluded with a warning. Default is `NULL`.
#' @param se A numeric vector of the same length as `theta` representing the
#'   standard errors associated with each ability estimate. All standard errors
#'   must be positive. If `NULL` and `x` is supplied, standard errors are
#'   computed using the test information function. See the **Details** section
#'   for more information. When `weights` is supplied, `se` must have one value
#'   per quadrature point. Standard errors from [irtQ::est_score()] that are
#'   set to 99.9999 (ability estimates at a limit of `range`) should be handled
#'   before they are supplied.
#' @param weights A two-column data frame or matrix. The first column holds the
#'   quadrature points (nodes), and the second column holds the corresponding
#'   weights, which should sum to 1; weights that do not sum to 1 are rescaled
#'   to sum to 1 with a warning (for example, frequencies can be given as
#'   weights). The weights must be nonnegative and finite and have a positive
#'   sum. [irtQ::gen.weight()] creates such a data frame. When `weights`
#'   is supplied, `theta` is ignored. Default is `NULL`.
#' @param D A scaling constant used in IRT models to make the logistic function
#'   closely approximate the normal ogive function. A value of 1.702 is commonly
#'   used for this purpose. Default is 1. It is used only when `se` is computed
#'   from `x`.
#'
#' @details This function first validates the input arguments. If both `theta`
#' and `weights` are `NULL`, the function will stop and return an error message.
#' Either `theta` or `weights` must be specified.
#'
#' Either `x` or `se` must be specified. If `se` is not provided (i.e., `se = NULL`),
#' it will be computed using the test information derived from the item metadata `x`.
#' The length of `se` must match the length of `theta`, or the number of quadrature
#' points in `weights`.
#'
#' It then computes the probability that an examinee with a given ability is
#' classified into each performance level using the normal distribution function
#' centered at each `theta` (or quadrature point) with standard deviation `se`.
#' The conditional classification accuracy is the probability of being
#' classified into the level that contains the ability value itself (Rudner,
#' 2001, 2005). The conditional classification consistency is the sum of the
#' squared level probabilities, as in Lee (2010); Rudner (2001, 2005) defines
#' only accuracy indices. When individual ability estimates are used, each
#' estimate is treated as the true ability, and each examinee receives the
#' weight 1/N, where N is the number of ability estimates after the examinees
#' with a missing `theta` or `se` are excluded.
#'
#' Finally, the function computes marginal classification accuracy and
#' consistency across all examinees by aggregating the conditional indices with
#' the associated weights.
#'
#' @return A list with the following elements:
#'  - `confusion`: A K x K matrix of expected proportions, with the true levels
#'    in rows (dimension name `True`) and the expected (observed) levels in
#'    columns (dimension name `Expected`). The entries are rounded to 7 decimal
#'    places and sum to 1 up to rounding.
#'  - `marginal`: A data frame with the columns `level`, `accuracy`, and
#'    `consistency`. The rows for levels 1 to K give the contribution of the
#'    ability values whose true level is that level, that is, the weighted sums
#'    of their conditional indices. These rows add up to the last row,
#'    `marginal`, which holds the marginal classification accuracy and
#'    consistency.
#'  - `conditional`: A data frame with one row per ability value and the columns
#'    `theta`, `weights`, `level` (true level, an integer), `accuracy`, and
#'    `consistency`.
#'  - `prob.level`: A data frame with the columns `theta`, `weights`, and
#'    `level`, followed by `p.level.1` to `p.level.K`, the probabilities of being
#'    classified into each level.
#'  - `cutscore`: The cut scores used, on the ability (theta) metric.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [irtQ::gen.weight()], [irtQ::est_score()], [irtQ::cac_lee()]
#'
#' @references Lee, W.-C. (2010). Classification consistency and accuracy for
#'   complex assessments using item response theory. *Journal of Educational
#'   Measurement, 47*(1), 1-17. \doi{10.1111/j.1745-3984.2009.00096.x}.
#'
#'   Rudner, L. M. (2001). Computing the expected proportions of
#'   misclassified examinees. *Practical Assessment, Research & Evaluation,
#'   7*(14). \doi{10.7275/an9m-2035}.
#'
#'   Rudner, L. M. (2005). Expected classification accuracy. *Practical
#'   Assessment, Research & Evaluation, 10*(13). \doi{10.7275/56a5-6b14}.
#'
#' @examples
#' \donttest{
#' ## --------------------------------------------------------------------------
#' ## 1. Using a population ability distribution (D method)
#' ## --------------------------------------------------------------------------
#'
#' # Import the "-prm.txt" output file from flexMIRT
#' flex_prm <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
#'
#' # Read item parameter estimates and convert them into item metadata
#' x <- bring.flexmirt(file = flex_prm, "par")$Group1$full_df
#'
#' # Define cut scores on the theta scale
#' cutscore <- c(-2, -0.5, 0.8)
#'
#' # Create quadrature points and corresponding weights
#' node <- seq(-4, 4, 0.25)
#' weights <- gen.weight(dist = "norm", mu = 0, sigma = 1, theta = node)
#'
#' # Compute classification accuracy and consistency
#' cac_1 <- cac_rud(
#'   x = x,
#'   cutscore = cutscore,
#'   weights = weights,
#'   se = NULL,
#'   D = 1)
#' print(cac_1)
#'
#' ## --------------------------------------------------------------------------
#' ## 2. Using individual ability estimates (P method)
#' ## --------------------------------------------------------------------------
#'
#' # Generate true abilities from N(0, 1)
#' set.seed(12)
#' theta <- rnorm(n = 1000, mean = 0, sd = 1)
#'
#' # Simulate item response data
#' data <- simdat(x = x, theta = theta, D = 1)
#'
#' # Estimate ability and standard errors using ML estimation
#' est_theta <- est_score(
#'   x = x, data = data, D = 1, method = "ML",
#'   range = c(-4, 4), se = TRUE
#' )
#' theta_hat <- est_theta$est.theta
#' se <- est_theta$se.theta
#'
#' # Replace the SEs of the estimates at a bound of range (set to 99.9999)
#' # with the SEs from the test information function
#' at_bound <- se == 99.9999
#' se[at_bound] <- 1 / sqrt(info(x = x, theta = theta_hat[at_bound], D = 1, tif = TRUE)$tif)
#'
#' # Compute classification accuracy and consistency using provided SEs
#' cac_2 <- cac_rud(
#'   cutscore = cutscore,
#'   theta = theta_hat,
#'   se = se)
#' print(cac_2)
#'
#' # Or compute classification accuracy and consistency using the item metadata
#' # instead of providing the SEs directly
#' cac_2 <- cac_rud(
#'   x = x,
#'   cutscore = cutscore,
#'   theta = theta_hat)
#' print(cac_2)
#' }
#'
#' @import dplyr
#' @export
cac_rud <- function(x = NULL,
                    cutscore,
                    theta = NULL,
                    se = NULL,
                    weights = NULL,
                    D = 1) {

  # check if the provided inputs are correct
  if (is.null(theta) && is.null(weights)) {
    stop("Either 'theta' or 'weights' argument must be provided; both cannot be NULL",
         call. = FALSE
    )
  }

  # stop when the cut scores are not finite values in strictly ascending order
  if (!is.numeric(cutscore) || length(cutscore) == 0L ||
      !all(is.finite(cutscore)) || is.unsorted(cutscore, strictly = TRUE)) {
    stop("'cutscore' must be a numeric vector of finite values in strictly ascending order.",
         call. = FALSE
    )
  }

  # exclude the examinees whose ability estimate or standard error is missing
  if (is.null(weights)) {
    # check if the provided inputs are correct
    if (!is.null(se) && length(se) != length(theta)) {
      stop("The numbers of thetas and the standard errors must be equal.",
           call. = FALSE
      )
    }

    # find the missing ability estimates and standard errors
    na.lg <- is.na(theta)
    if (!is.null(se)) na.lg <- na.lg | is.na(se)
    if (all(na.lg)) {
      stop("All values in 'theta' (or 'se') are missing.", call. = FALSE)
    }
    if (any(na.lg)) {
      warning(sprintf(
        "%d examinee(s) with a missing value in 'theta' or 'se' were excluded.",
        sum(na.lg)
      ), call. = FALSE)
      theta <- theta[!na.lg]
      if (!is.null(se)) se <- se[!na.lg]
    }
  }

  # compute standard errors if not provided
  if (is.null(se)) {
    if (is.null(x)) {
      stop("Either `se` or `x` argument must be supplied.", call. = FALSE)
    }
    if (!is.null(weights)) {
      se <- 1 / sqrt(info(x = x, theta = weights[, 1], D = D, tif = TRUE)$tif)
    } else {
      se <- 1 / sqrt(info(x = x, theta = theta, D = D, tif = TRUE)$tif)
    }
  }

  # stop when a standard error is missing or not positive
  if (anyNA(se) || any(se <= 0)) {
    stop("All standard errors in 'se' must be positive.", call. = FALSE)
  }

  # count the number of levels
  n.lev <- length(cutscore) + 1

  if (!is.null(weights)) {
    # (1) when the quadrature points and the corresponding ses are provided
    # extract nodes and weights
    nodes <- weights[, 1]
    wts <- weights[, 2]

    # check if the provided inputs are correct
    if (length(wts) != length(se)) {
      stop("The numbers of weights and the standard errors must be equal.",
           call. = FALSE
      )
    }

    # stop when the weights are not nonnegative finite values with a positive sum
    if (!all(is.finite(wts)) || any(wts < 0) || sum(wts) <= 0) {
      stop("The weights must be nonnegative finite values with a positive sum.",
           call. = FALSE)
    }

    # rescale the weights when they do not sum to one
    if (abs(sum(wts) - 1) > 1e-8) {
      warning("The weights do not sum to 1 and were rescaled to sum to 1.",
              call. = FALSE
      )
      wts <- wts / sum(wts)
    }
  } else {
    # (2) when individual ability estimates and ses are provided
    # use the ability estimates as nodes with uniform weights
    nodes <- theta
    wts <- rep(1 / length(theta), length(theta))
  }

  # count the number of thetas
  n.theta <- length(wts)

  # assign the levels to each theta; a theta equal to a cut score belongs to
  # the higher level
  level <- findInterval(nodes, cutscore) + 1L

  # the cumulative probability across all cut scores, from the normal
  # distribution centered at each theta with sd equal to its se
  cum_ps <- matrix(
    vapply(
      X = c(-Inf, cutscore, Inf),
      FUN = function(q) stats::pnorm(q = q, mean = nodes, sd = se),
      FUN.VALUE = numeric(n.theta)
    ),
    nrow = n.theta
  )

  # the probability that each examinee with each ability value is assigned to
  # each level category
  ps_tb <- cum_ps[, -1, drop = FALSE] - cum_ps[, -(n.lev + 1), drop = FALSE]
  colnames(ps_tb) <- paste0("p.level.", 1:n.lev)

  # conditional classification accuracy and consistency for each theta value
  cond_tb <- data.frame(
    theta = nodes,
    weights = wts,
    level = level,
    accuracy = ps_tb[cbind(1:n.theta, level)],
    consistency = rowSums(ps_tb^2)
  )

  # compute the marginal accuracy and consistency
  margin_tb <-
    cond_tb %>%
    # keep empty performance levels by fixing the factor levels
    dplyr::mutate(level = factor(.data$level, levels = seq_len(n.lev))) %>%
    dplyr::group_by(.data$level, .drop = FALSE) %>%
    dplyr::summarise(
      accuracy = sum(.data$accuracy * .data$weights),
      consistency = sum(.data$consistency * .data$weights),
      .groups = "drop"
    ) %>%
    janitor::adorn_totals(where = "row", name = "marginal")

  # add more variables to ps_tb
  ps_tb2 <-
    data.frame(
      theta = nodes, weights = wts,
      level = level, ps_tb
    )

  # create a cross table between true and expected levels
  cross_tb <-
    ps_tb2 %>%
    # keep empty performance levels by fixing the factor levels
    dplyr::mutate(level = factor(.data$level, levels = seq_len(n.lev))) %>%
    dplyr::group_by(.data$level, .drop = FALSE) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::starts_with("p.level."),
        ~ {
          sum(.x * .data$weights)
        }
      ),
      .groups = "drop"
    ) %>%
    dplyr::arrange(.data$level) %>%
    tibble::column_to_rownames("level") %>%
    data.matrix()
  dimnames(cross_tb) <- list(True = 1:n.lev, Expected = 1:n.lev)

  # return the results
  rst <- list(
    confusion = round(cross_tb, 7),
    marginal = margin_tb,
    conditional = cond_tb,
    prob.level = ps_tb2,
    cutscore = cutscore
  )
  return(rst)
}
