#' The Bisection Method to Find a Root
#'
#' This function is a modified version of the `bisection` function
#' in the \pkg{cmna} R package (Howard, 2017), designed to find a root of the function
#' `.fun` with respect to its first argument. Unlike the original `bisection()` in
#' \pkg{cmna}, this version allows additional arguments to be passed to `.fun`.
#'
#' @param .fun A function for which the root is to be found.
#' @param ... Additional arguments to be passed to `.fun`.
#' @param lb A numeric value specifying the lower bound of the search interval.
#' @param ub A numeric value specifying the upper bound of the search interval.
#' @param tol A numeric value specifying the tolerance for convergence. Default is 1e-4.
#' @param max.it An integer specifying the maximum number of iterations. Default is 100.
#'
#' @details The bisection method is a root-finding algorithm for a continuous
#' function whose values at the lower (`lb`) and upper (`ub`) bounds have
#' opposite signs. The method repeatedly halves the interval until its width
#' is no greater than `tol` or the maximum number of iterations (`max.it`) is
#' reached, and returns the midpoint of the final interval. A warning is issued
#' when the function values at the bounds do not have opposite signs or when
#' `max.it` is reached before the interval is no wider than `tol`.
#'
#' @return A list with the following components:
#' - `root`: The estimated root of the function.
#' - `iter`: The number of iterations performed.
#' - `delta`: The width of the final search interval.
#'
#' @seealso [irtQ::est_score()]
#'
#' @references Howard, J. P. (2017). *Computational methods for numerical
#' analysis with R*. New York: Chapman and Hall/CRC.
#'
#' @examples
#' ## Example: Find the theta value corresponding to a given probability
#' ## of a correct response using the item response function of a 2PLM
#' ## (a = 1, b = 0.2)
#'
#' # Define a function of theta
#' find.th <- function(theta, p) {
#'   p - drm(theta = theta, a = 1, b = 0.2, D = 1)
#' }
#'
#' # Find the theta value corresponding to p = 0.2
#' bisection(.fun = find.th, p = 0.2, lb = -10, ub = 10)$root
#'
#' # Find the theta value corresponding to p = 0.8
#' bisection(.fun = find.th, p = 0.8, lb = -10, ub = 10)$root
#'
#' @export
bisection <- function(.fun, ..., lb, ub, tol = 1e-4, max.it = 100) {
  iter <- 0
  f.ub <- .fun(ub, ...)

  # warn when the function values at the two bounds do not have opposite signs
  if (isTRUE(.fun(lb, ...) * f.ub > 0)) {
    warning("The function values at 'lb' and 'ub' must have opposite signs.", call. = FALSE)
  }

  while (abs(lb - ub) > tol) {
    # stop at the maximum number of iterations when the interval is still wider than tol
    if (iter >= max.it) {
      warning("The maximum number of iterations is reached.", call. = FALSE)
      break
    }
    mb <- (lb + ub) / 2
    f.mb <- .fun(mb, ...)
    if (f.mb == 0) {
      ub <- mb + (tol / 2)
      lb <- mb - (tol / 2)
    } else if (f.ub * f.mb < 0) {
      lb <- mb
    } else {
      ub <- mb
      f.ub <- f.mb
    }
    iter <- iter + 1
  }
  root <- (lb + ub) / 2
  list(root = root, iter = iter, delta = abs(lb - ub))
}
