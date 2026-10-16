#' Dichotomous Response Model (DRM) Probabilities
#'
#' This function computes the probability of a correct response for multiple items
#' given a set of theta values using the 1PL, 2PL, or 3PL item response models.
#'
#' @param theta A numeric vector of ability values (latent traits).
#' @param a A numeric vector of item discrimination (slope) parameters.
#' @param b A numeric vector of item difficulty parameters.
#' @param g A numeric vector of item guessing parameters. Not required for 1PL
#'   or 2PL models. `NA` values are treated as 0.
#' @param D A scaling constant used in IRT models to make the logistic function
#'   closely approximate the normal ogive function. A value of 1.702 is commonly
#'   used for this purpose. Default is 1.
#'
#' @details
#' If `g` is not specified, the function assumes a guessing parameter of 0 for all items,
#' corresponding to the 1PL or 2PL model. The function automatically adjusts the model
#' form based on the presence of `g`.
#'
#' @return
#' A matrix of response probabilities, where rows represent ability values (theta)
#' and columns represent items.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [irtQ::prm()]
#'
#' @examples
#' ## Example 1: theta and item parameters for 3PL model
#' drm(c(-0.1, 0.0, 1.5), a = c(1, 2), b = c(0, 1), g = c(0.2, 0.1), D = 1)
#'
#' ## Example 2: single theta value with 2PL item parameters
#' drm(0.0, a = c(1, 2), b = c(0, 1), D = 1)
#'
#' ## Example 3: multiple theta values with a single item (3PL model)
#' drm(c(-0.1, 0.0, 1.5), a = 1, b = 1, g = 0.2, D = 1)
#'
#' @importFrom Rfast Outer
#'
#' @export
drm <- function(theta, a, b, g = NULL, D = 1) {
  # count the numbers of items
  nitem <- length(a)

  # check the item guessing parameters
  if (is.null(g)) g <- rep(0, nitem)

  # treat NA guessing parameters as zeros
  g[is.na(g)] <- 0

  # calculate probability of correct answer (rows are items, columns are thetas)
  z <- (D * a) * Rfast::Outer(x = theta, y = b, oper = "-")
  P <- g + (1 - g) / (1 + exp(-z))

  # keep the probabilities away from 1 and from g (0 for the 1PLM and 2PLM)
  P[P > 9999999999e-10] <- 9999999999e-10
  # compare each item with its own guessing parameter before transposing
  lg.lessg <- P < (g + 1e-10)
  P[lg.lessg] <- P[lg.lessg] + 1e-10

  # return the probability matrix with thetas in rows and items in columns
  t(P)
}


#' Polytomous Response Model (PRM) Probabilities (GRM and GPCM)
#'
#' This function computes the probability of selecting each response category
#' for an item, given a set of theta values, using the graded response model
#' (GRM) or the (generalized) partial credit model (GPCM).
#'
#' @inheritParams drm
#' @param a A numeric value of the item discrimination (slope) parameter.
#' @param d A numeric vector of the threshold parameters of the item. For the
#'   GRM, the thresholds must be in increasing order.
#' @param pr.model A character string specifying the polytomous IRT model:
#'   `"GRM"` for the graded response model or `"GPCM"` for the (generalized)
#'   partial credit model. Default is `"GRM"`.
#'
#' @details When computing category probabilities using the partial credit model
#' (PCM), set `a = 1`.
#'
#' For `pr.model = "GPCM"`, `d` contains the threshold parameters
#' \eqn{b_v = \beta - \tau_v}, that is, the overall item location minus the
#' threshold of each score category. An item with *K* score categories requires
#' *K - 1* threshold parameters. The term for the lowest score category is
#' fixed at 0 and is not supplied. For example, a GPCM item with five
#' categories requires four threshold parameters.
#'
#' The GRM probabilities are bounded to the range from 1e-10 to 1 - 1e-10,
#' whereas the GPCM probabilities are bounded below only by the smallest
#' positive normalized double value (about 2.2e-308). The category
#' probabilities in [irtQ::traceline()] and [irtQ::info()] are bounded below by
#' 1e-10 for both models.
#'
#' For more details on the parameterization of the (generalized) partial credit
#' model, refer to the *IRT Models* section in the [irtQ-package] documentation.
#'
#' @return A matrix of category response probabilities, where each row
#' corresponds to a theta value and each column represents a score category of
#' the item.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [irtQ::drm()]
#'
#' @examples
#' ## Category probabilities for an item with four categories
#' ## using the generalized partial credit model (GPCM)
#' prm(theta = c(-0.2, 0, 0.5), a = 1.4, d = c(-0.2, 0, 0.5), D = 1, pr.model = "GPCM")
#'
#' ## Category probabilities for an item with five categories
#' ## using the graded response model (GRM)
#' prm(theta = c(-0.2, 0, 0.5), a = 1.2, d = c(-0.4, -0.2, 0.4, 1.5), D = 1, pr.model = "GRM")
#'
#' @export
prm <- function(theta, a, d, D = 1, pr.model = c("GRM", "GPCM")) {
  # use the first listed model when the default vector is left unchanged
  pr.model <- match.arg(toupper(pr.model), choices = c("GRM", "GPCM"))
  if (pr.model == "GRM") {
    P <- grm(theta = theta, a = a, d = d, D = D)
  } else {
    P <- gpcm(theta = theta, a = a, d = d, D = D)
  }

  # return
  P
}

# IRT GPC model
#' @importFrom Rfast Outer colCumSums rowsums
gpcm <- function(theta, a, d, D = 1) {
  # prepend 0 as the term for the lowest score category
  d <- c(0, d)

  # calculate category probabilities
  z <- (D * a) * (Rfast::Outer(x = theta, y = d, oper = "-"))
  cumsum_z <- t(Rfast::colCumSums(z))
  # subtract the row maximum so that exp() does not overflow
  cumsum_z <- cumsum_z - Rfast::rowMaxs(cumsum_z, value = TRUE)
  numer <- exp(cumsum_z) # numerator
  denom <- Rfast::rowsums(numer) # denominator
  P <- numer / denom

  # keep the probabilities above zero so that their logarithms are finite
  P[P < .Machine$double.xmin] <- .Machine$double.xmin

  # return
  P
}

# IRT GRM model
grm <- function(theta, a, d, D = 1) {
  # calculate the probabilities of scoring at or above each threshold
  allP <- drm(theta = theta, a = a, b = d, g = 0, D = D)

  # calculate category probabilities
  P <- cbind(1, allP) - cbind(allP, 0)
  P[P > 9999999999e-10] <- 9999999999e-10
  P[P < 1e-10] <- 1e-10

  # return
  P
}
