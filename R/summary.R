#' Summary of Item Calibration Results
#'
#' Summarizes the item calibration results in an object of class `est_irt`,
#' `est_mg`, or `est_item`, returned by [irtQ::est_irt()], [irtQ::est_mg()],
#' and [irtQ::est_item()], respectively.
#'
#' @param object An object of class `est_irt`, `est_mg`, or `est_item`.
#' @param ... Additional arguments passed to or from other methods (currently
#'   not used).
#'
#' @details \pkg{irtQ} exports its own generic `summary()`, which calls the
#'   method of the class of `object` and therefore works as [base::summary()]
#'   for other objects. Printing a summary object rounds the tables of item and
#'   group parameter estimates to the number of decimal places given in the
#'   `digits` argument of the print method (default
#'   `max(2L, getOption("digits") - 5L)`).
#'
#' @return An object of class `summary.est_irt`, `summary.est_mg`, or
#'   `summary.est_item`: a list of selected components of `object` (the call,
#'   the numbers of items and examinees, the estimation settings and
#'   convergence messages, the log-likelihood, AIC and BIC (not for
#'   `est_item`), and the item and group parameter estimates). Printing the
#'   object displays an overview of the calibration results.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [irtQ::est_irt()], [irtQ::est_mg()], [irtQ::est_item()]
#' @examples
#' \donttest{
#' # Fit the 1PL model to the LSAT6 data, constraining the slopes to be equal
#' fit.1pl <- est_irt(data = LSAT6, D = 1, model = "1PLM", cats = 2,
#'                    fix.a.1pl = FALSE, verbose = FALSE)
#'
#' # Display the calibration summary
#' summary(fit.1pl)
#' }
#'
#' @order 1
#' @export
summary <- function(object, ...) UseMethod("summary")

#' @describeIn summary An object created by the function [irtQ::est_irt()].
#' @order 2
#' @export
summary.est_irt <- function(object, ...) {
  call.expr <- deparse(object$call)
  nitem <- object$nitem
  ncase <- object$ncase
  MaxE <- object$MaxE
  Etol <- object$Etol
  weights <- object$weights
  npar.est <- object$npar.est
  npar.detail <- object$npar.detail
  fix.loc <- object$fix.loc
  niter <- object$niter
  maxpar.diff <- object$maxpar.diff
  EMtime <- object$EMtime
  SEtime <- object$SEtime
  TotalTime <- object$TotalTime
  test.1 <- object$test.1
  test.2 <- object$test.2
  var.note <- object$var.note
  loglikelihood <- object$loglikelihood
  aic <- object$aic
  bic <- object$bic
  estimates <- object$estimates
  group.par <- object$group.par

  out <- list(
    call.expr = call.expr, nitem = nitem, ncase = ncase, MaxE = MaxE, Etol = Etol,
    weights = weights, npar.est = npar.est, npar.detail = npar.detail, fix.loc = fix.loc, niter = niter,
    maxpar.diff = maxpar.diff, EMtime = EMtime, SEtime = SEtime, TotalTime = TotalTime,
    test.1 = test.1, test.2 = test.2, var.note = var.note, loglikelihood = loglikelihood,
    aic = aic, bic = bic, estimates = estimates, group.par = group.par
  )
  class(out) <- "summary.est_irt"
  out
}

#' @describeIn summary An object created by the function [irtQ::est_mg()].
#' @order 3
#' @export
summary.est_mg <- function(object, ...) {
  call.expr <- deparse(object$call)
  nitem <- object$nitem
  ncase <- object$ncase
  MaxE <- object$MaxE
  Etol <- object$Etol
  weights <- object$weights
  npar.est <- object$npar.est
  npar.detail <- object$npar.detail
  niter <- object$niter
  maxpar.diff <- object$maxpar.diff
  EMtime <- object$EMtime
  SEtime <- object$SEtime
  TotalTime <- object$TotalTime
  test.1 <- object$test.1
  test.2 <- object$test.2
  var.note <- object$var.note
  loglikelihood <- object$loglikelihood
  aic <- object$aic
  bic <- object$bic
  estimates <- object$estimates
  group.par <- object$group.par
  group.name <- names(object$loglikelihood$group)
  ngroup <- length(group.name)
  fix.loc <- object$fix.loc
  if (is.null(fix.loc$overall)) {
    fix.loc$group <- purrr::map(.x = 1:ngroup, ~ {
      NULL
    })
  }
  names(fix.loc$group) <- group.name

  out <- list(
    call.expr = call.expr, nitem = nitem, ncase = ncase, MaxE = MaxE, Etol = Etol,
    weights = weights, npar.est = npar.est, npar.detail = npar.detail, fix.loc = fix.loc, niter = niter,
    maxpar.diff = maxpar.diff, EMtime = EMtime, SEtime = SEtime, TotalTime = TotalTime,
    test.1 = test.1, test.2 = test.2, var.note = var.note, loglikelihood = loglikelihood,
    aic = aic, bic = bic, estimates = estimates, group.par = group.par, group.name = group.name,
    ngroup = ngroup
  )
  class(out) <- "summary.est_mg"
  out
}


#' @describeIn summary An object created by the function [irtQ::est_item()].
#' @order 4
#' @export
summary.est_item <- function(object, ...) {
  call.expr <- deparse(object$call)
  nitem <- object$nitem
  deleted.item <- object$deleted.item
  nitem.del <- length(deleted.item)
  npar.est <- object$npar.est
  n.response <- object$n.response
  TotalTime <- object$TotalTime
  convergence <- object$convergence
  loglikelihood <- object$loglikelihood
  estimates <- object$estimates
  group.par <- object$group.par

  out <- list(
    call.expr = call.expr, nitem = nitem, deleted.item = deleted.item, nitem.del = nitem.del, npar.est = npar.est,
    n.response = n.response, TotalTime = TotalTime, convergence = convergence, loglikelihood = loglikelihood,
    estimates = estimates, group.par = group.par
  )
  class(out) <- "summary.est_item"
  out
}

#' Summarize a Combined CTT Analysis
#'
#' Prepares the full, detailed report for an object of class `"ctt"`
#' returned by [ctt()]. Following the convention used by
#' [irtQ::est_irt()]/`summary.est_irt()`, this method returns an object of
#' class `"summary.ctt"` whose own `print.summary.ctt()` method displays the
#' complete item-level table, test-level reliability summary, and total-score
#' frequency distribution.
#'
#' @param object An object of class `"ctt"`, as returned by [ctt()].
#' @param ... Additional arguments passed to or from other methods (currently
#'   not used).
#'
#' @return An object of class `"summary.ctt"`: a list with the same `item`,
#'   `crit`, `alpha`, `freq`, and `call` elements as `object` (see
#'   [ctt()]'s **Value**), to be displayed by `print.summary.ctt()`.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [ctt()], [print.ctt()]
#'
#' @examples
#' # simulate the responses of 300 examinees to 15 dichotomous 3PLM items
#' set.seed(1)
#' x <- shape_df(
#'   par.drm = list(a = rep(1.5, 15), b = seq(-1.5, 1.5, length.out = 15),
#'                  g = rep(0.2, 15)),
#'   cats = 2, model = "3PLM"
#' )
#' dat <- simdat(x = x, theta = rnorm(300), D = 1)
#'
#' # run the CTT analysis
#' out <- ctt(data = dat)
#'
#' # create the full report object and print it
#' summary(out)
#'
#' @export
summary.ctt <- function(object, ...) {

  # ctt() already computes every element needed for the full report, so
  # summary.ctt() simply re-classes the same content for dispatch to
  # print.summary.ctt(), mirroring how summary.est_irt() repackages est_irt()
  # output rather than recomputing anything
  out <- list(item = object$item, crit = object$crit, alpha = object$alpha,
              freq = object$freq, call = object$call)
  class(out) <- "summary.ctt"
  out
}
