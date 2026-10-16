#' Classification Accuracy and Consistency Using Lee's (2010) Approach
#'
#' This function computes classification accuracy and consistency indices for
#' complex assessments using the item response theory (IRT) approach of Lee
#' (2010). Tests may contain dichotomous items, polytomous items, or both.
#'
#' @param x A data frame containing item metadata (e.g., item parameters, number
#'   of categories, IRT model types, etc.). See [irtQ::est_irt()] or
#'   [irtQ::simdat()] for more details about the item metadata. This data frame
#'   can be easily created using the [irtQ::shape_df()] function.
#' @param cutscore A numeric vector of finite cut scores in strictly ascending
#'   order. The K - 1 cut scores divide examinees into K performance levels.
#'   With `cut.obs = TRUE`, the cut scores are on the observed summed score
#'   metric. With `cut.obs = FALSE`, they are on the ability (theta) metric. A
#'   score equal to a cut score is assigned to the higher level.
#' @param theta A numeric vector of ability estimates, one per examinee, used by
#'   the P method. It is ignored when `weights` is supplied. Missing values
#'   (`NA`) are excluded with a warning. Default is `NULL`.
#' @param weights A two-column data frame or matrix used by the D method. The
#'   first column holds the quadrature points (nodes), and the second column
#'   holds the corresponding weights, which should sum to 1; weights that do
#'   not sum to 1 are rescaled to sum to 1 with a warning (for example,
#'   frequencies can be given as weights). The weights must be nonnegative and
#'   finite and have a positive sum. [irtQ::gen.weight()] creates such a
#'   data frame. When `weights` is supplied, `theta` is ignored. Default is
#'   `NULL`.
#' @param D A scaling constant used in IRT models to make the logistic function
#'   closely approximate the normal ogive function. A value of 1.702 is commonly
#'   used for this purpose. Default is 1.
#' @param cut.obs Logical. If `TRUE`, `cutscore` is on the observed summed score
#'   metric. If `FALSE`, `cutscore` is on the ability (theta) metric and is
#'   converted to the summed score metric with the test characteristic curve
#'   (TCC). Default is `TRUE`.
#'
#' @details
#' This function first validates the input arguments. If both `theta` and `weights`
#' are `NULL`, the function will stop and return an error message. Either `theta`
#' (a vector of ability estimates) or `weights` (a quadrature-based weight matrix)
#' must be specified.
#'
#' If `cut.obs = FALSE`, the provided cut scores are assumed to be on the theta (ability)
#' scale, and they are internally converted to the observed summed score scale using the
#' test characteristic curve (TCC). This transformation allows classification to be carried
#' out on the summed score metric, even if theta-based cut points are provided.
#'
#' For each ability value (a quadrature point in the D method, or an ability
#' estimate in the P method), the conditional distribution of the observed
#' summed score is computed with the Lord-Wingersky recursion (Lord &
#' Wingersky, 1984) and its extension to polytomous items (see
#' [irtQ::lwrc()]). Summing this distribution over the score ranges defined by
#' the cut scores gives the probability of being classified into each level. The
#' true level of an ability value is found by comparing its expected summed
#' score (the TCC value) with the same cut scores, so the observed cut scores
#' also serve as the true cut scores. The conditional classification accuracy is
#' the probability of being classified into the true level, and the conditional
#' classification consistency is the sum of the squared level probabilities,
#' which is the probability of the same classification on two independent
#' administrations (Lee, 2010).
#'
#' In the D method (`weights` supplied), the marginal indices are the weighted
#' sums of the conditional indices over the quadrature points. In the P method
#' (`theta` supplied), each examinee receives the weight 1/N, where N is the
#' number of ability estimates (after the missing values are excluded), so the
#' marginal indices are the averages of the conditional indices.
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
#'    `theta`, `weights`, `true.score` (expected summed score), `level` (true
#'    level, a factor), `accuracy`, and `consistency`.
#'  - `prob.level`: A data frame with the columns `theta`, `weights`,
#'    `true.score`, and `level`, followed by `p.level.1` to `p.level.K`, the
#'    probabilities of being classified into each level.
#'  - `cutscore`: The cut scores used, on the observed summed score metric. When
#'    `cut.obs = FALSE`, these are the converted cut scores.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [irtQ::gen.weight()], [irtQ::est_score()], [irtQ::cac_rud()]
#'
#' @references Lee, W.-C. (2010). Classification consistency and accuracy for
#'   complex assessments using item response theory. *Journal of Educational
#'   Measurement, 47*(1), 1-17. \doi{10.1111/j.1745-3984.2009.00096.x}.
#'
#'   Lord, F. M., & Wingersky, M. S. (1984). Comparison of IRT true-score and
#'   equipercentile observed-score "equatings." *Applied Psychological
#'   Measurement, 8*(4), 453-461. \doi{10.1177/014662168400800409}.
#'
#' @examples
#' \donttest{
#' ## --------------------------------------------------------------------------
#' ## 1. When the ability distribution is given by quadrature points and
#' ##    weights (D method)
#' ## --------------------------------------------------------------------------
#'
#' # Import the "-prm.txt" output file from flexMIRT
#' flex_prm <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
#'
#' # Read item parameter data and convert it to item metadata
#' x <- bring.flexmirt(file = flex_prm, "par")$Group1$full_df
#'
#' # Set the cut scores on the observed summed score scale
#' cutscore <- c(10, 20, 30, 50)
#'
#' # Create a data frame containing the quadrature points and corresponding weights
#' node <- seq(-4, 4, 0.25)
#' weights <- gen.weight(dist = "norm", mu = 0, sigma = 1, theta = node)
#'
#' # Calculate classification accuracy and consistency
#' cac_1 <- cac_lee(x = x, cutscore = cutscore, weights = weights, D = 1)
#' print(cac_1)
#'
#' ## -------------------------------------------------------------
#' ## 2. When individual ability estimates are available (P method)
#' ## -------------------------------------------------------------
#'
#' # Randomly draw true ability values from N(0, 1)
#' set.seed(12)
#' theta <- rnorm(n = 1000, mean = 0, sd = 1)
#'
#' # Simulate item response data
#' data <- simdat(x = x, theta = theta, D = 1)
#'
#' # Estimate ability parameters using maximum likelihood (ML)
#' est_th <- est_score(
#'   x = x, data = data, D = 1, method = "ML",
#'   range = c(-4, 4), se = FALSE
#' )$est.theta
#'
#' # Calculate classification accuracy and consistency
#' cac_2 <- cac_lee(x = x, cutscore = cutscore, theta = est_th, D = 1)
#'
#' # Marginal indices and confusion matrix
#' cac_2$marginal
#' cac_2$confusion
#'
#' ## ---------------------------------------------------------
#' ## 3. When individual ability estimates are available,
#' ##    but cut scores are specified on the IRT theta scale
#' ## ---------------------------------------------------------
#' # Set the cut scores on the theta scale
#' cutscore <- c(-2, -0.4, 0.2, 1.0)
#'
#' # Calculate classification accuracy and consistency
#' cac_3 <- cac_lee(
#'   x = x, cutscore = cutscore, theta = est_th, D = 1,
#'   cut.obs = FALSE
#' )
#'
#' # Marginal indices and confusion matrix
#' cac_3$marginal
#' cac_3$confusion
#' }
#'
#' @import dplyr
#' @export
cac_lee <- function(x,
                    cutscore,
                    theta = NULL,
                    weights = NULL,
                    D = 1,
                    cut.obs = TRUE) {

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

  # if the cutscores are on the theta metric, compute the expected cutscores
  # on the observed score metric
  if (!cut.obs) {
    cutscore <- irtQ::traceline(x = x, theta = cutscore, D = D)$tcc
  }

  # count the number of levels
  n.lev <- length(cutscore) + 1

  if (!is.null(weights)) {
    # (1) D method: quadrature points and weights are provided
    # extract nodes and weights
    nodes <- weights[, 1]
    wts <- weights[, 2]

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
    # (2) P method: individual ability estimates are provided
    # find the missing ability estimates
    na.lg <- is.na(theta)
    if (all(na.lg)) {
      stop("All values in 'theta' are missing.", call. = FALSE)
    }

    # exclude the missing ability estimates
    if (any(na.lg)) {
      warning(sprintf("%d missing value(s) in 'theta' were excluded.", sum(na.lg)),
              call. = FALSE
      )
      theta <- theta[!na.lg]
    }

    # use the ability estimates as nodes with uniform weights
    nodes <- theta
    wts <- rep(1 / length(theta), length(theta))
  }

  # count the number of thetas
  n.theta <- length(wts)

  # estimate likelihood functions using the Lord-Wingersky recursion
  # for each theta, which is the conditional observed score distribution
  # given each theta; each column is the conditional score distribution for
  # a given theta
  lkhd <- lwrc(x = x, theta = nodes, prob = NULL, D = D)

  # count the total number of observed scores for the test
  n.score <- nrow(lkhd)

  # check the maximum possible observed score
  max.score <- n.score - 1

  # compute the expected true score using the IRT models for each theta
  tscore <- traceline(x = x, theta = nodes, D = D)$tcc

  # assign the levels to each expected true score; a score equal to a cut
  # score belongs to the higher level
  level <- factor(findInterval(tscore, cutscore) + 1L, levels = 1:n.lev)

  # assign the level to each possible observed score
  loc.lev <- findInterval(0:max.score, cutscore) + 1L

  # the probability that each examinee with each ability value is assigned to
  # each level category; a level that contains no observed score has
  # probability 0
  ps_tb <- crossprod(lkhd, outer(loc.lev, 1:n.lev, "==") * 1)

  # name the level columns and drop the row names taken from the theta columns
  dimnames(ps_tb) <- list(NULL, paste0("p.level.", 1:n.lev))

  # conditional classification accuracy and consistency for each theta value
  cond_tb <- data.frame(
    theta = nodes,
    weights = wts,
    true.score = tscore,
    level = level,
    accuracy = ps_tb[cbind(1:n.theta, as.integer(level))],
    consistency = rowSums(ps_tb^2)
  )

  # compute the marginal accuracy and consistency
  margin_tb <-
    cond_tb %>%
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
      true.score = tscore, level = level,
      ps_tb
    )

  # create a cross table between true and expected levels
  cross_tb <-
    ps_tb2 %>%
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
