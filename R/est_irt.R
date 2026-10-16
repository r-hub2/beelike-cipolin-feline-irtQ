#' Item parameter estimation using MMLE-EM algorithm
#'
#' This function fits unidimensional item response theory (IRT) models to
#' mixed-format data comprising both dichotomous and polytomous items, using
#' marginal maximum likelihood estimation via the expectation-maximization
#' (MMLE-EM) algorithm (Bock & Aitkin, 1981). When prior distributions are
#' specified for item parameters, the estimates are the modes of the marginal
#' posterior distribution (Bayes modal estimates). The function also supports
#' fixed item parameter calibration (FIPC; Kim, 2006), a practical method for
#' pretest (or newly developed) item calibration in computerized adaptive
#' testing (CAT). FIPC enables the parameter estimates of pretest items to be
#' placed on the same scale as those of operational items (Ban et al., 2001).
#' For dichotomous items, the function supports the one-, two-, and
#' three-parameter logistic models. For polytomous items, it supports the
#' graded response model (GRM) and the (generalized) partial credit model
#' (GPCM).
#'
#' @param x A data frame containing item metadata (item IDs, numbers of score
#'   categories, IRT models, and item parameters); see Details for its format.
#'   You can create an empty item metadata frame using the function
#'   [irtQ::shape_df()]. The item parameters in `x` are used as starting values
#'   when `use.startval = TRUE` and are otherwise ignored, except for the fixed
#'   items in FIPC, whose parameters are held at the values in `x`. If
#'   `x = NULL`, the models and numbers of score categories are taken from
#'   `model` and `cats`. When `fipc = TRUE`, `x` must be supplied. Default is
#'   `NULL`.
#' @param data A matrix of examinees' item responses, with rows for examinees
#'   and columns for items in the same order as the items in `x`. Responses
#'   to an item with K score categories must be integers from 0 to K - 1.
#'   Missing responses are coded as `NA` or, in functions that have a
#'   `missing` argument, as the value given in `missing`.
#' @param D A scaling constant used in IRT models to make the logistic function
#'   closely approximate the normal ogive function. A value of 1.702 is commonly
#'   used for this purpose. Default is 1.
#' @param model A character vector specifying the IRT model to fit each item.
#'   Available values are:
#'   - `"1PLM"`, `"2PLM"`, `"3PLM"`, `"DRM"` for dichotomous items
#'   - `"GRM"`, `"GPCM"` for polytomous items
#'
#'   Here, `"GRM"` denotes the graded response model and `"GPCM"` the
#'   (generalized) partial credit model. Note that `"DRM"` is a general label
#'   for dichotomous items; in this function, items labeled `"DRM"` are fitted
#'   with the 3PLM, and a warning is issued. If a single model name is provided,
#'   it is recycled for all items. This argument is only used when `x = NULL`
#'   and `fipc = FALSE`. Default is `NULL`.
#' @param cats Numeric vector specifying the number of score categories per
#'   item. For dichotomous items, this should be 2. If a single value is
#'   supplied, it will be recycled across all items. When `cats = NULL` and all
#'   models specified in the `model` argument are dichotomous (`"1PLM"`,
#'   `"2PLM"`, `"3PLM"`, or `"DRM"`), the function defaults to 2 categories per
#'   item. This argument is used only when `x = NULL` and `fipc = FALSE`.
#'   Default is `NULL`.
#' @param item.id Character vector of item identifiers. If `NULL`, the IDs in
#'   `x` are used, or IDs are generated automatically when `x = NULL`. A
#'   provided `item.id` replaces any IDs in `x`, with or without FIPC, and must
#'   have one value per item. With `fipc = TRUE`, `fix.id` then refers to the
#'   new IDs. Default is `NULL`.
#' @param fix.a.1pl Logical. If `TRUE`, the slope parameters of all 1PLM items
#'   are fixed to `a.val.1pl`; otherwise, they are constrained to be equal and
#'   estimated. Default is `FALSE`.
#' @param fix.a.gpcm Logical. If `TRUE`, GPCM items are calibrated as PCM with
#'   slopes fixed to `a.val.gpcm`; otherwise, each item's slope is estimated.
#'   Default is `FALSE`.
#' @param fix.g Logical. If `TRUE`, all 3PLM guessing parameters are fixed to
#'   `g.val`; otherwise, each guessing parameter is estimated. The prior in
#'   `gprior` is then not used. Default is `FALSE`.
#' @param a.val.1pl Numeric. Value to which the slope parameters of 1PLM items
#'   are fixed when `fix.a.1pl = TRUE`. Default is 1.
#' @param a.val.gpcm Numeric. Value to which the slope parameters of GPCM items
#'   are fixed when `fix.a.gpcm = TRUE`. Default is 1.
#' @param g.val Numeric. Value to which the guessing parameters of 3PLM items
#'   are fixed when `fix.g = TRUE`. Default is 0.2.
#' @param use.aprior Logical. If `TRUE`, applies a prior distribution to all
#'   item discrimination (slope) parameters during calibration. Default is
#'   `FALSE`.
#' @param use.bprior Logical. If `TRUE`, applies a prior distribution to all
#'   item difficulty (or threshold) parameters during calibration. Default is
#'   `FALSE`.
#' @param use.gprior Logical. If `TRUE`, applies the prior distribution in
#'   `gprior` to the guessing parameters of 3PLM items when they are estimated;
#'   the prior is not used when `fix.g = TRUE`. Default is `TRUE`.
#' @param aprior,bprior,gprior A list specifying the prior distribution for all
#'   item discrimination (slope), difficulty (or threshold), and guessing
#'   parameters. The slope prior is placed on `D` times the slope, that is,
#'   the density in `aprior` is evaluated at `D * a`; with `D = 1` (the
#'   default), this is the slope itself. The priors also enter the information
#'   matrix from which the standard errors are computed. Three distributions
#'   are supported: Beta, Log-normal, and Normal. The list must have two
#'   elements:
#'   - `dist`: A character string, one of `"beta"`, `"lnorm"`, or `"norm"`.
#'   - `params`: A numeric vector of length two giving the distribution's
#'   parameters. For details on each parameterization, see [stats::dbeta()],
#'   [stats::dlnorm()], and [stats::dnorm()].
#'
#'   Defaults are:
#'   - `aprior = list(dist = "lnorm", params = c(0.0, 0.5))`
#'   - `bprior = list(dist = "norm", params = c(0.0, 1.0))`
#'   - `gprior = list(dist = "beta", params = c(5, 16))`
#'
#'   for discrimination, difficulty, and guessing parameters, respectively.
#' @param missing A value indicating missing responses in the data set, in
#'   addition to `NA`, which is always treated as missing. Missing responses are
#'   omitted from the likelihood, which assumes that they are missing at random.
#'   Default is `NA`.
#' @param Quadrature A numeric vector of length two:
#'   - first element: number of quadrature points
#'   - second element: symmetric bound (absolute value) for those points
#'
#'   For example, `c(49, 6)` specifies 49 evenly spaced points from -6 to 6.
#'   The points are used to approximate the integral over the latent ability
#'   distribution in the E-step of the EM algorithm. `Quadrature` is ignored
#'   when `weights` is supplied. Default is `c(49, 6)`.
#' @param weights A two-column matrix or data frame of quadrature points (first
#'   column) and weights (second column) that defines the latent ability
#'   distribution; it can be created with [irtQ::gen.weight()]. The weights
#'   should sum to 1; weights that do not sum to 1 are rescaled to sum to 1
#'   with a warning. The number of points may differ from `Quadrature[1]`. If
#'   supplied, `Quadrature`, `group.mean`, and `group.var` are not used.
#'   Without FIPC, the mean and variance of this distribution fix the scale:
#'   the distribution is used as it is when `EmpHist = FALSE`, and its shape is
#'   updated but rescaled to its mean and variance when `EmpHist = TRUE`, in
#'   which case the quadrature points must be equally spaced, because the
#'   histogram is rescaled by interpolation on an equally spaced grid. With
#'   FIPC, it is only the initial distribution of the first EM cycle, which is
#'   then updated (see Details).
#'
#'   If `NULL`, a normal distribution with mean `group.mean` and variance
#'   `group.var` at the points given by `Quadrature` is used. Default is
#'   `NULL`.
#' @param group.mean A numeric value specifying the mean of the normal latent
#'   ability distribution when `weights = NULL`. Without FIPC, the mean is fixed
#'   at this value to identify the scale of the item parameters, also when
#'   `EmpHist = TRUE`. With FIPC, it is the mean of the initial distribution,
#'   and the mean is then estimated. Default is 0.
#' @param group.var A positive numeric value specifying the variance of the
#'   normal latent ability distribution when `weights = NULL`. Without FIPC,
#'   the variance is fixed at this value to identify the scale of the item
#'   parameters, also when `EmpHist = TRUE`. With FIPC, it is the variance of
#'   the initial distribution, and the variance is then estimated. Default is
#'   1.
#' @param EmpHist Logical. If `TRUE`, the latent ability distribution is
#'   estimated as an empirical histogram on the quadrature points together with
#'   the item parameters (Woods, 2007). Item calibration is conducted relative
#'   to the estimated empirical histogram. See Details. Default is `FALSE`.
#' @param use.startval Logical. If `TRUE`, the item parameters provided in the
#'   item metadata (i.e., the `x` argument) are used as starting values for the
#'   items whose parameters are estimated; `x` must then be supplied.
#'   Otherwise, internally generated starting values are used. Default is
#'   `FALSE`.
#' @param Etol A positive numeric value specifying the convergence criterion of
#'   the EM algorithm. The EM cycles stop when the largest absolute change in
#'   the estimated item parameters between two consecutive cycles is less than
#'   or equal to `Etol`. With FIPC, only the parameters of the items that are
#'   not fixed are compared. The mean and variance of the latent ability
#'   distribution estimated with FIPC, and the weights of an empirical
#'   histogram (`EmpHist = TRUE`), are not part of the criterion, so they can
#'   be less converged than the item parameters; a smaller `Etol` gives more
#'   precise estimates of them. In particular, an estimated histogram can keep
#'   changing very slowly over many cycles (toward a spiky distribution) while
#'   the log-likelihood barely changes, so `group.par` and `weights` then
#'   depend on `Etol`; with small samples, the item parameters can then also
#'   need many cycles, or more than `MaxE`, to meet `Etol`. When all items are
#'   fixed, the criterion is applied to the largest absolute change in the mean
#'   and variance of the latent ability distribution. `Etol` is not used when
#'   `fipc.method = "OEM"` with items to calibrate. The printed output of the
#'   fitted object refers to this criterion as the convergence criterion of the
#'   E-step. Default is `1e-4`.
#' @param MaxE A positive integer specifying the maximum number of EM cycles. It
#'   is set to 1 when `fipc.method = "OEM"`. The printed output of the fitted
#'   object refers to the cycles as E-step cycles. Default is `500`.
#' @param control A named list of options passed directly to [stats::nlminb()],
#'   which maximizes the M-step objective of each item. The default is
#'   `list(eval.max = 500, iter.max = 200, x.tol = 1e-4)`, where
#'   - `eval.max` is the maximum number of objective function evaluations,
#'   - `iter.max` is the maximum number of optimizer iterations, and
#'   - `x.tol` is the X-convergence tolerance of [stats::nlminb()].
#'
#'   A partial list may be supplied; options that are not supplied keep their
#'   default values. Other options of [stats::nlminb()], such as `rel.tol` or
#'   `trace`, can also be supplied.
#' @param fipc Logical. If `TRUE`, fixed item parameter calibration (FIPC) is
#'   applied during item parameter estimation. When `fipc = TRUE`, the
#'   information on which items are fixed must be provided via either `fix.loc`
#'   or `fix.id`. See Details.
#' @param fipc.method A character string specifying the FIPC method, used only
#'   when `fipc = TRUE`. Available options are:
#'   - `"OEM"`: No Prior Weights Updating and One EM Cycle (NWU-OEM; Wainer &
#'   Mislevy, 1990).
#'   - `"MEM"`: Multiple Prior Weights Updating and Multiple EM Cycles
#'   (MWU-MEM; Kim, 2006).
#'
#'   When `fipc.method = "OEM"`, a single EM cycle is run regardless of `MaxE`.
#'   See Details. Default is `"MEM"`.
#' @param fix.loc A vector of positive integers giving the row positions in `x`
#'   of the items to be fixed when `fipc = TRUE`. For example, suppose that five
#'   items located in the 1st, 2nd, 4th, 7th, and 9th rows of `x` should be
#'   fixed. Then use `fix.loc = c(1, 2, 4, 7, 9)`. It must not contain
#'   duplicated positions. If `fix.id` is also supplied, `fix.loc` is ignored
#'   with a warning. See Details.
#' @param fix.id A character vector of the IDs of the items to be fixed when
#'   `fipc = TRUE`. The IDs are matched to the first column (`id`) of `x`, or to
#'   `item.id` when it is supplied. For example, use
#'   `fix.id = c("CMC1", "CMC2", "CMC3", "CMC4", "CMC5")` to fix the items with
#'   these IDs. If `fix.id` is not `NULL`, `fix.loc` is ignored. See Details.
#' @param se Logical. If `FALSE`, the standard errors and the
#'   variance-covariance matrix of the item parameter estimates are not
#'   computed; `se.est` is then filled with `NA` and `covariance` is `NULL`.
#'   Default is `TRUE`.
#' @param verbose Logical. If `FALSE`, all progress messages, including
#'   information about the EM algorithm process, are suppressed. Default is
#'   `TRUE`.
#' @details A specific format of data frame should be used for the argument `x`.
#' The first column should contain item IDs, the second column should contain
#' the number of score categories of each item, and the third column
#' should specify the IRT model to be fitted to each item. Available IRT models
#' are:
#'   - `"1PLM"`, `"2PLM"`, `"3PLM"`, and `"DRM"` for dichotomous item data
#'   - `"GRM"` and `"GPCM"` for polytomous item data
#'
#' Note that `"DRM"` serves as a general label covering all dichotomous IRT
#' models (i.e., `"1PLM"`, `"2PLM"`, and `"3PLM"`), while `"GRM"` and `"GPCM"`
#' represent the graded response model and (generalized) partial credit model,
#' respectively. In `est_irt()`, items labeled `"DRM"` are calibrated with the
#' 3PLM.
#'
#' The subsequent columns should contain the item parameters for the specified
#' models. For dichotomous items, the fourth, fifth, and sixth columns represent
#' item discrimination (slope), item difficulty, and item guessing parameters,
#' respectively. When `"1PLM"` or `"2PLM"` is specified in the third column,
#' `NA`s must be entered in the sixth column for the guessing parameters.
#'
#' For polytomous items, the item discrimination (slope) parameter should appear
#' in the fourth column, and the item difficulty (or threshold) parameters for
#' category boundaries should occupy the fifth through the last columns. When
#' the number of score categories differs across items, unused parameter
#' cells should be filled with `NA`s.
#'
#' In the \pkg{irtQ} package, the threshold parameters for GPCM items are
#' expressed as the item location (or overall difficulty) minus the threshold
#' values for each score category. Note that a GPCM item with *K* score
#' categories has *K - 1* threshold parameters, because the term for the lowest
#' score category (score 0) is fixed at 0 and is not included in the metadata.
#' For example, if a GPCM item has five score categories, four threshold
#' parameters must be provided.
#'
#' An example of a data frame for a single-format test is shown below:
#' \tabular{lrlrrr}{
#'   id \tab cats \tab model \tab par.1 \tab par.2 \tab par.3 \cr
#'   ITEM1 \tab 2 \tab 1PLM \tab 1.000 \tab 1.461 \tab NA \cr
#'   ITEM2 \tab 2 \tab 2PLM \tab 1.921 \tab -1.049 \tab NA \cr
#'   ITEM3 \tab 2 \tab 3PLM \tab 1.736 \tab 1.501 \tab 0.203 \cr
#'   ITEM4 \tab 2 \tab 3PLM \tab 0.835 \tab -1.049 \tab 0.182 \cr
#'   ITEM5 \tab 2 \tab DRM \tab 0.926 \tab 0.394 \tab 0.099
#' }
#'
#' An example of a data frame for a mixed-format test is shown below:
#' \tabular{lrlrrrrr}{
#'   id \tab cats \tab model \tab par.1 \tab par.2 \tab par.3 \tab par.4
#'   \tab par.5 \cr
#'   ITEM1 \tab 2 \tab 1PLM \tab 1.000 \tab 1.461 \tab NA \tab NA \tab NA \cr
#'   ITEM2 \tab 2 \tab 2PLM \tab 1.921 \tab -1.049 \tab NA \tab NA \tab NA \cr
#'   ITEM3 \tab 2 \tab 3PLM \tab 0.926 \tab 0.394 \tab 0.099 \tab NA \tab NA \cr
#'   ITEM4 \tab 2 \tab DRM \tab 1.052 \tab -0.407 \tab 0.201 \tab NA \tab NA \cr
#'   ITEM5 \tab 4 \tab GRM \tab 1.913 \tab -1.869 \tab -1.238 \tab -0.714 \tab
#'   NA \cr
#'   ITEM6 \tab 5 \tab GRM \tab 1.278 \tab -0.724 \tab -0.068 \tab 0.568 \tab
#'   1.072 \cr
#'   ITEM7 \tab 4 \tab GPCM \tab 1.137 \tab -0.374 \tab 0.215 \tab 0.848 \tab
#'   NA \cr
#'   ITEM8 \tab 5 \tab GPCM \tab 1.233 \tab -2.078 \tab -1.347 \tab -0.705 \tab
#'   -0.116
#' }
#'
#' See the *IRT Models* section in the [irtQ-package] documentation for more
#' details about the IRT models used in the \pkg{irtQ} package. A convenient way
#' to create a data frame for the argument `x` is by using the function
#' [irtQ::shape_df()].
#'
#' To fit IRT models to data, the item response data must be accompanied by
#' information on the IRT model and the number of score categories for each
#' item. There are two ways to provide this information:
#'    1. Supply item metadata to the argument `x`. As explained above, such
#'    metadata can be easily created using [irtQ::shape_df()].
#'    2. Specify the IRT models and score category information directly
#'    through the arguments `model` and `cats`.
#'
#' If `x` is supplied, `model` and `cats` are ignored. If `x = NULL`, the
#' function uses `model` and `cats`; `cats` may be omitted when all items are
#' dichotomous.
#'
#' Item parameters are estimated with the EM algorithm of Bock and Aitkin
#' (1981). In the E-step, the posterior distribution of each examinee over the
#' quadrature points is computed from the current item parameters and the
#' current latent distribution, and the expected numbers of examinees and of
#' responses in each score category are obtained at each quadrature point. In
#' the M-step, the expected complete-data log-likelihood of each item, plus the
#' log prior densities when priors are used, is maximized with
#' [stats::nlminb()] using analytic gradients and Hessians; the 1PLM items with
#' a common slope are estimated jointly. The cycles are repeated until the
#' convergence criterion given by `Etol` is met or `MaxE` cycles are completed.
#' Missing responses are omitted from the likelihood. The slopes of dichotomous
#' items are not bounded during estimation, whereas those of GRM and GPCM items
#' are kept positive.
#'
#' To implement FIPC, the item metadata must be provided via the `x` argument.
#' This is because the item parameters of the fixed items in the metadata are
#' used to estimate the characteristics of the underlying latent variable prior
#' distribution when calibrating the remaining (freely estimated) items. More
#' specifically, the latent prior distribution is estimated based on the
#' fixed items, and then used to calibrate the new (pretest) items so that their
#' parameters are placed on the same scale as those of the fixed items (Kim,
#' 2006). The full item metadata, including both fixed and non-fixed items, can
#' be conveniently created using the [irtQ::shape_df_fipc()] function.
#'
#' In terms of approaches for FIPC, Kim (2006) described five different methods.
#' Among them, two methods are available in the [irtQ::est_irt()] function. The
#' first method is `"NWU-OEM"`, which uses a single E-step in the EM
#' algorithm (involving only the fixed items) followed by a single M-step
#' (involving only the non-fixed items). This method was proposed by Wainer and
#' Mislevy (1990) in the context of online calibration and can be implemented by
#' setting `fipc.method = "OEM"`. The item parameters are estimated with the
#' initial latent distribution (given by `group.mean` and `group.var`, or by
#' `weights`). The latent distribution returned in `group.par` and `weights` is
#' computed once afterward from the posterior distributions based on the fixed
#' items, and `EmpHist` affects only this returned distribution. When all items
#' are fixed, the function always warns that the convergence criteria are not
#' satisfied.
#'
#' The second method is `"MWU-MEM"`, which runs multiple EM cycles and updates
#' the weights of the latent ability distribution in every cycle (Kim, 2006;
#' Kim & Kolen, 2016). The first cycle is the same as in NWU-OEM: the posterior
#' distribution of each examinee is computed from the responses to the fixed
#' items only. From the second cycle onward, the posterior distribution is
#' computed from the responses to all items, using the fixed parameters and the
#' current estimates of the other items, and the latent distribution is updated
#' from the expected numbers of examinees at the quadrature points. With
#' `EmpHist = TRUE`, the updated distribution is this empirical histogram, as in
#' Kim (2006); with `EmpHist = FALSE`, it is a normal distribution with the
#' mean and variance of the histogram. This method can be implemented by
#' setting `fipc.method = "MEM"`.
#'
#' When `fipc = TRUE`, information about which items are to be fixed must be
#' provided via either the `fix.loc` or `fix.id` argument. For example, suppose
#' that five items with IDs "CMC1", "CMC2", "CMC3", "CMC4", and "CMC5" should
#' be fixed, and all item IDs are provided via the `x` or `item.id` argument.
#' Also, assume these five items are located in the 1st through 5th rows of the
#' item metadata (i.e., `x`). In this case, the fixed items can be specified
#' using either `fix.loc = c(1, 2, 3, 4, 5)` or
#' `fix.id = c("CMC1", "CMC2", "CMC3", "CMC4", "CMC5")`. Note that if both
#' `fix.loc` and `fix.id` are not `NULL`, the information in `fix.loc` is
#' ignored.
#'
#' When `EmpHist = TRUE`, the latent distribution is an empirical histogram: in
#' each EM cycle, the weight at each quadrature point is updated to the average
#' of the posterior probabilities of the examinees at that point (Woods, 2007).
#' Without FIPC, the updated histogram is then rescaled so that its mean and
#' variance equal `group.mean` and `group.var` (or the mean and variance of
#' `weights`, if supplied): the quadrature points are standardized and
#' transformed to the target scale, the weights are interpolated back to the
#' original points (and extrapolated at the ends), and the weights are
#' renormalized. With FIPC, the histogram is not rescaled, so its mean and
#' variance are estimated. When `EmpHist = FALSE`, a normal distribution is
#' used. With FIPC, its mean and variance are updated in each cycle to those of
#' the histogram; otherwise, they are fixed at `group.mean` and `group.var`.
#'
#' A maximum likelihood estimate does not exist for an item with a score
#' category that has no responses, which includes an item whose responses are
#' all the same. The estimates of such an item can be extreme or remain near
#' their starting values, and their standard errors are unreliable (a standard
#' error of 99999 means that it could not be obtained). Before the estimation,
#' one warning lists the IDs of these items; the items fixed in FIPC are not
#' checked. A parameter with a prior distribution, such as the guessing
#' parameter of a 3PLM item with `use.gprior = TRUE`, can still be stable.
#'
#' When prior distributions are used, the item parameters are Bayes modal
#' estimates, and `loglikelihood` is the log-likelihood at these estimates
#' without the log prior densities.
#'
#' The information criteria are `aic = -2 * loglikelihood + 2 * npar.est` and
#' `bic = -2 * loglikelihood + npar.est * log(ncase)`. Here `loglikelihood` is
#' the marginal log-likelihood without the prior densities of the item
#' parameters, and `ncase` is the number of examinees. `npar.est` is the number
#' of freely estimated parameters, which consists of (a) the item parameters
#' that are estimated (a slope shared by 1PLM items is counted once; fixed
#' parameters and the items fixed in FIPC are not counted), (b) the mean and
#' variance of the latent distribution when they are estimated, which is the
#' case with FIPC and `EmpHist = FALSE` (2 values), and (c) the weights of the
#' empirical histogram when `EmpHist = TRUE`. Without FIPC, the histogram is
#' rescaled to a fixed mean and variance (given by `group.mean` and `group.var`,
#' or by `weights`), so its Q weights, where Q is the number of quadrature
#' points, have Q - 3 free values (the sum of the weights and the two moments
#' are fixed). With FIPC, the histogram is not rescaled, so the weights have
#' Q - 1 free values. The counts by type are returned in `npar.detail`.
#'
#'   The print method of this class has an argument `digits` that is not used.
#'   It prints the log-likelihood with the number of significant digits given
#'   by `getOption("digits")` (7 by default).
#'
#' @return This function returns an object of class `est_irt`. The returned
#'   object contains the following components:
#'
#'   \item{estimates}{A data frame with the item ID, number of score
#'   categories, and model of each item, followed by the item parameter
#'   estimates and their standard errors in alternating columns (`par.1`,
#'   `se.1`, `par.2`, `se.2`, ...).}
#'
#'   \item{par.est}{A data frame of item parameter estimates, structured
#'   according to the item metadata format.}
#'
#'   \item{se.est}{A data frame of standard errors for the item parameter
#'   estimates, computed from the inverse of the information matrix. The
#'   information matrix is the cross-product (empirical) information matrix of
#'   the marginal log-likelihood evaluated at the final estimates (Meilijson,
#'   1989), plus the prior information (the negative second derivatives of the
#'   log prior densities) for the parameters that have prior distributions.
#'   With FIPC, the latent ability distribution is treated as known in this
#'   matrix, so the uncertainty of the estimated distribution is not reflected
#'   in these standard errors. Standard errors are `NA` for parameters that are
#'   not estimated (for example, fixed slopes or guessing parameters, and the
#'   items fixed in FIPC) and for all parameters when `se = FALSE`. The
#'   standard error of a common slope of 1PLM items appears only in the row of
#'   the first 1PLM item. A value of 99999 indicates that the standard error
#'   could not be obtained.}
#'
#'   \item{pos.par}{A data frame indicating, for each estimated item parameter,
#'   its row (and column) position in `covariance`. Parameters that are not
#'   estimated are `NA`, and the position of a common slope of 1PLM items is
#'   given only in the row of the first 1PLM item. `NULL` when FIPC is used
#'   with all items fixed.}
#'
#'   \item{covariance}{The variance-covariance matrix of the item parameter
#'   estimates (the inverse of the information matrix described in `se.est`),
#'   without row or column names; the order of the parameters is given by
#'   `pos.par`. `NULL` when `se = FALSE` or when no item parameter is
#'   estimated.}
#'
#'   \item{loglikelihood}{The marginal log-likelihood of the response data at
#'   the final estimates: the sum over examinees of the logarithm of the
#'   probability of each examinee's response pattern, integrated over the
#'   latent ability distribution with the quadrature points and weights. With
#'   FIPC, the responses to all items, including the fixed items, are used. The
#'   log prior densities of the item parameters are not included.}
#'
#'   \item{aic}{Akaike information criterion,
#'   `-2 * loglikelihood + 2 * npar.est`.}
#'
#'   \item{bic}{Bayesian information criterion,
#'   `-2 * loglikelihood + npar.est * log(ncase)`.}
#'
#'   \item{group.par}{A data frame with columns `mu`, `sigma2`, and `sigma`
#'   (mean, variance, and standard deviation of the latent ability
#'   distribution) and rows `estimates` and `se`. Without FIPC, the values are
#'   fixed at `group.mean` and `group.var` (or at the mean and variance of
#'   `weights`, if supplied) and `se` is `NA`. When `EmpHist = TRUE`, the
#'   empirical histogram is rescaled in each cycle to have this mean and
#'   variance. Because the rescaling interpolates on the quadrature points, the
#'   mean and variance computed from the returned `weights` can differ slightly
#'   from the reported values. With FIPC, the values are the mean and variance
#'   of the final latent distribution (these are the mean and variance of the
#'   returned `weights`), and `se` gives
#'   approximate standard errors, `sigma / sqrt(N)` for the mean and
#'   `sigma2 * sqrt(2 / (N - 1))` for the variance, with the delta method for
#'   `sigma`, where `N` is the number of examinees; these treat the abilities as
#'   observed and are therefore usually too small.}
#'
#'   \item{weights}{A data frame with columns `theta` and `weight` giving the
#'   quadrature points and the weights of the final latent ability
#'   distribution. Without FIPC and with `EmpHist = FALSE`, it is the
#'   distribution given by `group.mean` and `group.var`, or `weights` as
#'   supplied (rescaled to sum to 1 if needed). With `fipc.method = "OEM"`, it
#'   is the distribution updated once after the item parameters were
#'   estimated.}
#'
#'   \item{posterior.dist}{A matrix of the posterior probabilities of the
#'   quadrature points for each examinee, computed at the final estimates, with
#'   one row per examinee and one column per quadrature point; each row sums to
#'   1.}
#'
#'   \item{data}{A matrix of the examinees' response data, with the missing
#'   values coded as `NA`.}
#'
#'   \item{scale.D}{The scaling factor used in the IRT model.}
#'
#'   \item{ncase}{The number of examinees (rows of `data`).}
#'
#'   \item{nitem}{The total number of items in the response data.}
#'
#'   \item{Etol}{The convergence criterion of the EM algorithm.}
#'
#'   \item{MaxE}{The maximum number of EM cycles (1 when
#'   `fipc.method = "OEM"`).}
#'
#'   \item{aprior}{The prior distribution specified for the slope parameters
#'   when `use.aprior = TRUE`; `NULL` otherwise.}
#'
#'   \item{bprior}{The prior distribution specified for the difficulty (or
#'   threshold) parameters when `use.bprior = TRUE`; `NULL` otherwise.}
#'
#'   \item{gprior}{The prior distribution specified for the guessing parameters
#'   when `use.gprior = TRUE`; `NULL` otherwise.}
#'
#'   \item{npar.est}{The number of freely estimated parameters used for `aic`
#'   and `bic`: the estimated item parameters, the mean and variance of the
#'   latent distribution when they are estimated, and the weights of the
#'   empirical histogram when `EmpHist = TRUE`. See Details.}
#'
#'   \item{niter}{The number of completed EM cycles.}
#'
#'   \item{maxpar.diff}{The largest absolute change in the estimated item
#'   parameters in the last EM cycle (for FIPC with all items fixed, in the
#'   mean and variance of the latent ability distribution). With
#'   `fipc.method = "OEM"`, it is the change from the starting values.}
#'
#'   \item{EMtime}{Time (in seconds) spent on EM cycles.}
#'
#'   \item{SEtime}{Time (in seconds) spent computing standard errors; `NULL`
#'   when they are not computed.}
#'
#'   \item{TotalTime}{Total computation time (in seconds).}
#'
#'   \item{test.1}{A message indicating whether the convergence criteria were
#'   met: the M-step optimization converged for every item and the largest
#'   absolute change in the parameter estimates between two consecutive EM
#'   cycles was less than or equal to `Etol`. When `fipc.method = "OEM"` with
#'   new items to calibrate, the single EM cycle is not judged by `Etol` and
#'   only the M-step check applies. For FIPC with all items fixed, see `Etol`.}
#'
#'   \item{test.2}{A message indicating whether the information matrix used
#'   for the standard errors (see `se.est`) is positive definite. Printed
#'   results label it "Second-order test", but because the cross-product
#'   information matrix is positive semidefinite by construction, this check
#'   mainly detects a singular (rank-deficient) information matrix, for
#'   example due to unidentified parameters; it is not a second-order test of
#'   the log-likelihood surface. The message reports a possible local maximum
#'   only when the first-order test is also satisfied. When FIPC is used with
#'   all items fixed, no information matrix is computed, and the message says
#'   so.}
#'
#'   \item{var.note}{A note indicating whether the variance-covariance matrix
#'   was successfully obtained from the information matrix.}
#'
#'   \item{fipc}{Logical. Indicates whether FIPC was used.}
#'
#'   \item{fipc.method}{The FIPC method used; `NULL` when `fipc = FALSE`.}
#'
#'   \item{fix.loc}{An integer vector of the row positions of the fixed items
#'   in `x`, also when the items were specified with `fix.id`; `NULL` when
#'   `fipc = FALSE`.}
#'
#'   \item{call}{The matched call.}
#'
#'   \item{npar.detail}{An integer vector with the elements `item` (free item
#'   parameters), `dist` (mean and variance of the latent distribution), and
#'   `weight` (weights of the empirical histogram), whose sum is `npar.est`.
#'   It is the last element of the object; when it is absent, `summary()` shows
#'   the total only.}
#'
#'   Note that you can easily extract components from the output using the
#'   [irtQ::getirt()] function.
#'
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [irtQ::shape_df()], [irtQ::shape_df_fipc()], [irtQ::getirt()]
#'
#' @references Ban, J. C., Hanson, B. A., Wang, T., Yi, Q., & Harris, D. J.
#'   (2001). A comparative study of on-line pretest item calibration/scaling
#'   methods in computerized adaptive testing. *Journal of Educational
#'   Measurement, 38*(3), 191-212. \doi{10.1111/j.1745-3984.2001.tb01123.x}.
#'
#'   Bock, R. D., & Aitkin, M. (1981). Marginal maximum likelihood estimation of
#'   item parameters: Application of an EM algorithm. *Psychometrika, 46*,
#'   443-459.
#'
#'   Kim, S. (2006). A comparative study of IRT fixed parameter calibration
#'   methods. *Journal of Educational Measurement, 43*(4), 355-381.
#'
#'   Kim, S., & Kolen, M. J. (2016). *Multiple group IRT fixed-parameter
#'   estimation for maintaining an established ability scale* (CASMA Research
#'   Report No. 49). Center for Advanced Studies in Measurement and Assessment,
#'   University of Iowa.
#'
#'   Meilijson, I. (1989). A fast improvement to the EM algorithm on its own
#'   terms. *Journal of the Royal Statistical Society: Series B
#'   (Methodological), 51*, 127-138.
#'
#'   Wainer, H., & Mislevy, R. J. (1990). Item response theory, item
#'   calibration, and proficiency estimation. In H. Wainer (Ed.), *Computer
#'   adaptive testing: A primer* (Chap. 4, pp. 65-102). Hillsdale, NJ: Lawrence
#'   Erlbaum.
#'
#'   Woods, C. M. (2007). Empirical histograms in item response theory with
#'   ordinal data. *Educational and Psychological Measurement, 67*(1), 73-87.
#'
#' @examples
#' \donttest{
#'
#' ## --------------------------------------------------------------
#' ## 1. Item parameter estimation for dichotomous item data (LSAT6)
#' ## --------------------------------------------------------------
#' # Fit the 1PL model to LSAT6 data and estimate a common slope parameter
#' # (i.e., constrain slope parameters to be equal)
#' (mod.1pl.c <- est_irt(data = LSAT6, D = 1, model = "1PLM", cats = 2,
#'                       fix.a.1pl = FALSE))
#'
#' # Display a summary of the estimation results
#' summary(mod.1pl.c)
#'
#' # Extract the item parameter estimates
#' getirt(mod.1pl.c, what = "par.est")
#'
#' # Extract the standard error estimates
#' getirt(mod.1pl.c, what = "se.est")
#'
#' # Fit the 1PL model to LSAT6 data and fix slope parameters to 1.0
#' (mod.1pl.f <- est_irt(data = LSAT6, D = 1, model = "1PLM", cats = 2,
#'                       fix.a.1pl = TRUE, a.val.1pl = 1))
#'
#' # Display a summary of the estimation results
#' summary(mod.1pl.f)
#'
#' # Fit the 2PL model to LSAT6 data
#' (mod.2pl <- est_irt(data = LSAT6, D = 1, model = "2PLM", cats = 2))
#'
#' # Display a summary of the estimation results
#' summary(mod.2pl)
#'
#' # Assess the model fit for the 2PL model using the S-X2 fit statistic
#' (sx2fit.2pl <- sx2_fit(x = mod.2pl))
#'
#' # Compute item and test information functions at a range of theta values
#' theta <- seq(-4, 4, 0.1)
#' (info.2pl <- info(x = mod.2pl, theta = theta))
#'
#' # Plot the test characteristic curve (TCC)
#' (trace.2pl <- traceline(x = mod.2pl, theta = theta))
#' plot(trace.2pl)
#'
#' # Plot the item characteristic curve (ICC) for the first item
#' plot(trace.2pl, item.loc = 1)
#'
#' # Fit the 2PL model and simultaneously estimate an empirical histogram
#' # of the latent variable prior distribution
#' # Also apply a looser convergence threshold for the EM algorithm
#' (mod.2pl.hist <- est_irt(data = LSAT6, D = 1, model = "2PLM", cats = 2,
#'                          EmpHist = TRUE, Etol = 0.001))
#' (emphist <- getirt(mod.2pl.hist, what = "weights"))
#' plot(emphist$weight ~ emphist$theta, type = "h")
#'
#' # Fit the 3PL model and apply a Beta prior to the guessing parameters
#' (mod.3pl <- est_irt(
#'   data = LSAT6, D = 1, model = "3PLM", cats = 2, use.gprior = TRUE,
#'   gprior = list(dist = "beta", params = c(5, 16))
#' ))
#'
#' # Display a summary of the estimation results
#' summary(mod.3pl)
#'
#' # Fit the 3PL model and fix the guessing parameters at 0.2
#' (mod.3pl.f <- est_irt(data = LSAT6, D = 1, model = "3PLM", cats = 2,
#'                       fix.g = TRUE, g.val = 0.2))
#'
#' # Display a summary of the estimation results
#' summary(mod.3pl.f)
#'
#' # Fit different dichotomous models to each item in the LSAT6 data:
#' # Fit the constrained 1PL model to items 1-3, the 2PL model to item 4,
#' # and the 3PL model with a Beta prior on guessing to item 5
#' (mod.drm.mix <- est_irt(
#'   data = LSAT6, D = 1, model = c("1PLM", "1PLM", "1PLM", "2PLM", "3PLM"),
#'   cats = 2, fix.a.1pl = FALSE, use.gprior = TRUE,
#'   gprior = list(dist = "beta", params = c(5, 16))
#' ))
#'
#' # Display a summary of the estimation results
#' summary(mod.drm.mix)
#'
#' ## -------------------------------------------------------------------
#' ## 2. Item parameter estimation for mixed-format data (simulated data)
#' ## -------------------------------------------------------------------
#' ## Import the "-prm.txt" output file from flexMIRT
#' flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt",
#'                         package = "irtQ")
#'
#' # Extract item metadata
#' x <- bring.flexmirt(file = flex_sam, "par")$Group1$full_df
#'
#' # Modify the item metadata so that the 39th and 40th items use the GPCM
#' x[39:40, 3] <- "GPCM"
#'
#' # Generate 1,000 examinees' latent abilities from N(0, 1)
#' set.seed(37)
#' score1 <- rnorm(1000, mean = 0, sd = 1)
#'
#' # Simulate item response data
#' sim.dat1 <- simdat(x = x, theta = score1, D = 1)
#'
#' # Fit the 3PL model to all dichotomous items, the GPCM to items 39 and 40,
#' # and the GRM to items 53, 54, and 55.
#' # Use a Beta prior for guessing parameters, a log-normal prior for slope
#' # parameters, and a normal prior for difficulty (threshold) parameters.
#' # Also, specify the argument `x` to provide IRT model and score category
#' # information.
#' item.meta <- shape_df(item.id = x$id, cats = x$cats, model = x$model,
#'   default.par = TRUE)
#' (mod.mix1 <- est_irt(
#'   x = item.meta, data = sim.dat1, D = 1, use.aprior = TRUE,
#'   use.bprior = TRUE,
#'   use.gprior = TRUE,
#'   aprior = list(dist = "lnorm", params = c(0.0, 0.5)),
#'   bprior = list(dist = "norm", params = c(0.0, 2.0)),
#'   gprior = list(dist = "beta", params = c(5, 16))
#' ))
#'
#' # Display a summary of the estimation results
#' summary(mod.mix1)
#'
#' # Estimate examinees' latent scores using MLE and the estimated item
#' # parameters
#' (score.mle <- est_score(x = mod.mix1, method = "ML", range = c(-4, 4)))
#'
#' # Compute traditional model-fit statistics
#' (fit.mix1 <- irtfit(
#'   x = mod.mix1, score = score.mle$est.theta, group.method = "equal.width",
#'   n.width = 10, loc.theta = "middle"
#' ))
#'
#' # Residual plot for the first item (dichotomous)
#' plot(
#'   x = fit.mix1, item.loc = 1, type = "both", ci.method = "wald",
#'   show.table = TRUE, ylim.sr.adjust = TRUE
#' )
#'
#' # Residual plot for the last item (polytomous)
#' plot(
#'   x = fit.mix1, item.loc = 55, type = "both", ci.method = "wald",
#'   show.table = FALSE, ylim.sr.adjust = TRUE
#' )
#'
#' # Fit the 2PL model to all dichotomous items, the GPCM to items 39 and 40,
#' # and the GRM to items 53, 54, and 55.
#' # Provide IRT model and score category information via `model` and `cats`
#' # arguments.
#' (mod.mix2 <- est_irt(
#'   data = sim.dat1, D = 1,
#'   model = c(rep("2PLM", 38), rep("GPCM", 2), rep("2PLM", 12), rep("GRM", 3)),
#'   cats = c(rep(2, 38), rep(5, 2), rep(2, 12), rep(5, 3))
#' ))
#'
#' # Display a summary of the estimation results
#' summary(mod.mix2)
#'
#' # Fit the 2PL model to all dichotomous items, the GPCM to items 39 and 40,
#' # and the GRM to items 53, 54, and 55.
#' # Also estimate the empirical histogram of the latent prior distribution.
#' # Provide IRT model and score category information via `model` and `cats`
#' # arguments.
#' (mod.mix3 <- est_irt(
#'   data = sim.dat1, D = 1,
#'   model = c(rep("2PLM", 38), rep("GPCM", 2), rep("2PLM", 12), rep("GRM", 3)),
#'   cats = c(rep(2, 38), rep(5, 2), rep(2, 12), rep(5, 3)), EmpHist = TRUE
#' ))
#' (emphist <- getirt(mod.mix3, what = "weights"))
#' plot(emphist$weight ~ emphist$theta, type = "h")
#'
#' # Fit the 2PL model to all dichotomous items, the PCM to items 39 and 40 by
#' # fixing slope parameters to 1, and the GRM to items 53, 54, and 55.
#' # Provide IRT model and score category information via `model` and `cats`
#' # arguments.
#' (mod.mix4 <- est_irt(
#'   data = sim.dat1, D = 1,
#'   model = c(rep("2PLM", 38), rep("GPCM", 2), rep("2PLM", 12), rep("GRM", 3)),
#'   cats = c(rep(2, 38), rep(5, 2), rep(2, 12), rep(5, 3)),
#'   fix.a.gpcm = TRUE, a.val.gpcm = 1
#' ))
#'
#' # Display a summary of the estimation results
#' summary(mod.mix4)
#'
#' ## ----------------------------------------------------------------
#' ## 3. Fixed item parameter calibration (FIPC) for mixed-format data
#' ##    (simulated)
#' ## ----------------------------------------------------------------
#' ## Import the "-prm.txt" output file from flexMIRT
#' flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt",
#'                         package = "irtQ")
#'
#' # Select item metadata
#' x <- bring.flexmirt(file = flex_sam, "par")$Group1$full_df
#'
#' # Generate 1,000 examinees' latent abilities from N(0.4, 1.3)
#' set.seed(20)
#' score2 <- rnorm(1000, mean = 0.4, sd = 1.3)
#'
#' # Simulate response data
#' sim.dat2 <- simdat(x = x, theta = score2, D = 1)
#'
#' # Fit the 3PL model to all dichotomous items and the GRM to all polytomous
#' # items
#' # Fix five 3PL items (1st - 5th) and three GRM items (53rd - 55th)
#' # Also estimate the empirical histogram of the latent variable distribution
#' # Use the MEM method
#' fix.loc <- c(1:5, 53:55)
#' (mod.fix1 <- est_irt(
#'   x = x, data = sim.dat2, D = 1, use.gprior = TRUE,
#'   gprior = list(dist = "beta", params = c(5, 16)), EmpHist = TRUE,
#'   Etol = 1e-3, fipc = TRUE, fipc.method = "MEM", fix.loc = fix.loc
#' ))
#'
#' # Extract group-level parameter estimates
#' (prior.par <- mod.fix1$group.par)
#'
#' # Visualize the empirical prior distribution
#' (emphist <- getirt(mod.fix1, what = "weights"))
#' plot(emphist$weight ~ emphist$theta, type = "h")
#'
#' # Display a summary of the estimation results
#' summary(mod.fix1)
#'
#' # Alternatively, fix the same items by providing their item IDs
#' # using the `fix.id` argument. In this case, set `fix.loc = NULL`
#' fix.id <- c(x$id[1:5], x$id[53:55])
#' (mod.fix1 <- est_irt(
#'   x = x, data = sim.dat2, D = 1, use.gprior = TRUE,
#'   gprior = list(dist = "beta", params = c(5, 16)), EmpHist = TRUE,
#'   Etol = 1e-3, fipc = TRUE, fipc.method = "MEM", fix.loc = NULL,
#'   fix.id = fix.id
#' ))
#'
#' # Display a summary of the estimation results
#' summary(mod.fix1)
#'
#' # Fit the 3PL model to all dichotomous items and the GRM to all polytomous
#' # items
#' # Fix the same items as before (1st - 5th and 53rd - 55th)
#' # This time, do not estimate the empirical histogram of the latent prior
#' # Instead, estimate the scale of the normal prior distribution
#' # Use the MEM method
#' fix.loc <- c(1:5, 53:55)
#' (mod.fix2 <- est_irt(
#'   x = x, data = sim.dat2, D = 1, use.gprior = TRUE,
#'   gprior = list(dist = "beta", params = c(5, 16)), EmpHist = FALSE,
#'   Etol = 1e-3, fipc = TRUE, fipc.method = "MEM", fix.loc = fix.loc
#' ))
#'
#' # Extract group-level parameter estimates
#' (prior.par <- mod.fix2$group.par)
#'
#' # Visualize the prior distribution
#' (emphist <- getirt(mod.fix2, what = "weights"))
#' plot(emphist$weight ~ emphist$theta, type = "h")
#'
#' # Fit the 3PL model to all dichotomous items and the GRM to all polytomous
#' # items
#' # Fix only the five 3PL items (1st to 5th) and estimate the empirical
#' # histogram
#' # Use the OEM method (i.e., only one EM cycle is used)
#' fix.loc <- c(1:5)
#' (mod.fix3 <- est_irt(
#'   x = x, data = sim.dat2, D = 1, use.gprior = TRUE,
#'   gprior = list(dist = "beta", params = c(5, 16)), EmpHist = TRUE,
#'   Etol = 1e-3, fipc = TRUE, fipc.method = "OEM", fix.loc = fix.loc
#' ))
#'
#' # Extract group-level parameter estimates
#' (prior.par <- mod.fix3$group.par)
#'
#' # Visualize the prior distribution
#' (emphist <- getirt(mod.fix3, what = "weights"))
#' plot(emphist$weight ~ emphist$theta, type = "h")
#'
#' # Display a summary of the estimation results
#' summary(mod.fix3)
#'
#' # Fit the 3PL model to all dichotomous items and the GRM to all polytomous
#' # items
#' # Fix all 55 items and estimate only the latent ability distribution
#' # Use the MEM method
#' fix.loc <- c(1:55)
#' (mod.fix4 <- est_irt(
#'   x = x, data = sim.dat2, D = 1, EmpHist = TRUE,
#'   Etol = 1e-3, fipc = TRUE, fipc.method = "MEM", fix.loc = fix.loc
#' ))
#'
#' # Extract group-level parameter estimates
#' (prior.par <- mod.fix4$group.par)
#'
#' # Visualize the prior distribution
#' (emphist <- getirt(mod.fix4, what = "weights"))
#' plot(emphist$weight ~ emphist$theta, type = "h")
#'
#' # Display a summary of the estimation results
#' summary(mod.fix4)
#'
#' # Alternatively, fix all 55 items by providing their item IDs
#' # using the `fix.id` argument. In this case, set `fix.loc = NULL`
#' fix.id <- x$id
#' (mod.fix4 <- est_irt(
#'   x = x, data = sim.dat2, D = 1, EmpHist = TRUE,
#'   Etol = 1e-3, fipc = TRUE, fipc.method = "MEM", fix.loc = NULL,
#'   fix.id = fix.id
#' ))
#'
#' # Display a summary of the estimation results
#' summary(mod.fix4)
#'
#' }
#'
#' @importFrom utils modifyList
#'
#' @export
#'
est_irt <- function(x = NULL,
                    data,
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
                    aprior = list(dist = "lnorm", params = c(0.0, 0.5)),
                    bprior = list(dist = "norm", params = c(0.0, 1.0)),
                    gprior = list(dist = "beta", params = c(5, 16)),
                    missing = NA,
                    Quadrature = c(49, 6.0),
                    weights = NULL,
                    group.mean = 0.0,
                    group.var = 1.0,
                    EmpHist = FALSE,
                    use.startval = FALSE,
                    Etol = 1e-04,
                    MaxE = 500,
                    control = list(eval.max = 500, iter.max = 200, x.tol = 1e-4),
                    fipc = FALSE,
                    fipc.method = "MEM",
                    fix.loc = NULL,
                    fix.id = NULL,
                    se = TRUE,
                    verbose = TRUE) {
  
  # match.call
  cl <- match.call()

  # Merge user-supplied control list with defaults, enabling partial specification
  # (e.g., control = list(iter.max = 500) keeps eval.max and x.tol at defaults)
  default_control <- list(eval.max = 500, iter.max = 200, x.tol = 1e-4)
  control <- modifyList(default_control, control)

  # the FIPC method must be one of the two available methods
  if (fipc && !(length(fipc.method) == 1L && fipc.method %in% c("MEM", "OEM"))) {
    stop("'fipc.method' must be either \"MEM\" or \"OEM\".", call. = FALSE)
  }

  # item parameter estimation
  if (!fipc) {
    # item parameter estimation using MMLE-EM algorithm
    est_par <- est_irt_em(
      x = x, data = data, D = D, model = model, cats = cats, item.id = item.id,
      fix.a.1pl = fix.a.1pl, fix.a.gpcm = fix.a.gpcm, fix.g = fix.g,
      a.val.1pl = a.val.1pl, a.val.gpcm = a.val.gpcm, g.val = g.val,
      use.aprior = use.aprior, use.bprior = use.bprior, use.gprior = use.gprior,
      aprior = aprior, bprior = bprior, gprior = gprior, missing = missing,
      Quadrature = Quadrature, weights = weights, group.mean = group.mean,
      group.var = group.var, EmpHist = EmpHist, use.startval = use.startval,
      Etol = Etol, MaxE = MaxE, control = control, se = se, verbose = verbose
    )
  } else {
    # implement FIPC method
    est_par <- est_irt_fipc(
      x = x, data = data, D = D, item.id = item.id, fix.a.1pl = fix.a.1pl,
      fix.a.gpcm = fix.a.gpcm, fix.g = fix.g, a.val.1pl = a.val.1pl,
      a.val.gpcm = a.val.gpcm, g.val = g.val, use.aprior = use.aprior,
      use.bprior = use.bprior, use.gprior = use.gprior, aprior = aprior,
      bprior = bprior, gprior = gprior, missing = missing, Quadrature = Quadrature,
      weights = weights, group.mean = group.mean, group.var = group.var,
      EmpHist = EmpHist, use.startval = use.startval, Etol = Etol, MaxE = MaxE,
      control = control, fipc = TRUE, fipc.method = fipc.method, fix.loc = fix.loc,
      fix.id = fix.id, se = se, verbose = verbose
    )
  }
  
  # return the estimation results
  class(est_par) <- "est_irt"
  est_par$call <- cl

  # place the counts of the free parameters after the call
  npar.detail <- est_par$npar.detail
  est_par$npar.detail <- NULL
  est_par$npar.detail <- npar.detail
  est_par
}



# This function estimates item parameters using MMLE-EM algorithm
#' @importFrom Rfast colsums
#' @import dplyr
est_irt_em <- function(x = NULL,
                       data,
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
                       aprior = list(dist = "lnorm", params = c(0.0, 0.5)),
                       bprior = list(dist = "norm", params = c(0.0, 1.0)),
                       gprior = list(dist = "beta", params = c(5, 16)),
                       missing = NA,
                       Quadrature = c(49, 6.0),
                       weights = NULL,
                       group.mean = 0,
                       group.var = 1,
                       EmpHist = FALSE,
                       use.startval = FALSE,
                       Etol = 1e-04,
                       MaxE = 500,
                       control = list(eval.max = 500, iter.max = 200, x.tol = 1e-4),
                       se = TRUE,
                       verbose = TRUE) {
  
  # check start time
  start.time <- Sys.time()
  
  ## ---------------------------------------------------------------------
  # prepare the item parameter estimation
  ## ---------------------------------------------------------------------
  # check if the starting values are available
  if (use.startval & is.null(x)) {
    stop(paste0(
      "To use starting values for item parameter estimation, \n",
      "the item metadata must be specified in the argument 'x'."
    ), call. = FALSE)
  }
  
  # recode a non-numeric missing code before the responses are read as numbers
  if (!is.na(missing)) data[data == missing] <- NA

  # transform a data set to matrix (character and factor columns are read by their values)
  data <- resp_as_numeric(data)
  
  # extract information about the number of score categories and models
  if (verbose) {
    cat("Parsing input...", "\n")
  }
  if (!is.null(x)) {
    # confirm and correct all item metadata information
    x <- confirm_df(x)

    # override the item id in x if item.id is not NULL
    if (!is.null(item.id)) {
      # stop when item.id does not have one value per item of x
      if (length(item.id) != nrow(x)) {
        stop("The length of 'item.id' must be equal to the number of items in 'x'.", call. = FALSE)
      }
      x$id <- item.id
    }

    # extract item ids for all items
    id <- x$id
    
    # extract score categories for all items
    cats <- x$cats
    
    # extract model names for all items
    model <- x$model
    
    # if use.startval = FALSE
    # create a new item metadata object with starting values
    if (!use.startval) {
      x <- startval_df(cats = cats, model = model, item.id = id)
    }
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
    x <- startval_df(cats = cats, model = model, item.id = item.id)
    
    # confirm and correct all item metadata information
    x <- confirm_df(x)
    
    # extract item ids for all items
    id <- x$id
    
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
  
  # check the number of item responses across all items
  n.resp <- Rfast::colsums(!is.na(data))
  
  # check the items which have all missing responses
  loc_allmiss <- which(n.resp == 0L)
  if (length(loc_allmiss) > 0L) {
    memo2 <- paste0(paste0("item ", loc_allmiss, collapse = ", "), " has/have no item response data. \n")
    stop(memo2, call. = FALSE)
  }

  # stop when a response is not one of the score categories 0, ..., cats - 1
  check_resp_range(data = data, cats = cats, id = id)

  # warn about the items whose maximum likelihood estimates do not exist
  check_item_mle(data = data, cats = cats, id = id)

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
  
  ## ---------------------------------------------------------------------
  # conduct item parameter estimation using MMLE-EM algorithm
  ## ---------------------------------------------------------------------
  # create initial weights of prior ability distribution when it is not specified
  if (is.null(weights)) {
    # create a vector of quad-points
    quadpt <- seq(-Quadrature[2], Quadrature[2], length.out = Quadrature[1])
    
    # create a two column data frame to contain the quad-points and weights
    weights <- gen.weight(dist = "norm", mu = group.mean, sigma = sqrt(group.var), theta = quadpt)
    n.quad <- length(quadpt)
  } else {
    # check the weights and rescale them to sum to one, as the moments below, the
    # empirical histogram rescaling, and the marginal log-likelihood assume
    weights <- check_weights(weights)
    quadpt <- weights[, 1]
    n.quad <- length(quadpt)
    moments.tmp <- cal_moment(node = quadpt, weight = weights[, 2])
    group.mean <- moments.tmp[1]
    group.var <- moments.tmp[2]
  }
  
  # build the per-item one-hot frequency-category list used by
  # divide_data() below and by info_xpd() in the SE step.  See
  # build_freqcat() (R/util.R) for the output structure.
  freq.cat <- build_freqcat(data, cats)
  
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
  
  # divide the data set for the mixed-item format
  datlist <- divide_data(data = data, idx.item = idx.item, freq.cat = freq.cat)
  data_drm <- cbind(datlist$data_drm_q, datlist$data_drm_p)
  data_prm <- datlist$data_prm
  data_all <- datlist$data_all
  
  # delete 'datlist' object
  rm(datlist, envir = environment(), inherits = FALSE)
  
  # create the lower and upper bounds of the item parameters
  parbd <- lubound(model, cats, n.1PLM, idx4est, fix.a.1pl, fix.g, fix.a.gpcm)
  
  # find the columns of the frequency matrix corresponding to all items
  cols.item <- cols4item(nitem, cats, loc_1p_const)
  
  # estimation
  if (verbose) {
    cat("Estimating item parameters...", "\n")
  }
  
  # implement EM algorithm
  time1 <- Sys.time()
  for (r in 1:MaxE) {
    # implement E-step
    estep <- Estep(
      elm_item = elm_item, idx.drm = idx.drm, idx.prm = idx.prm,
      data_drm = data_drm, data_prm = data_prm, data_all = data_all,
      weights = weights, D = D
    )
    
    # implement M-step
    mstep <- Mstep(
      estep = estep, id = id, cats = cats, model = model, quadpt = quadpt,
      n.quad = n.quad, D = D, cols.item = cols.item, loc_1p_const = loc_1p_const,
      loc_else = loc_else, idx4est = idx4est, n.1PLM = n.1PLM, EmpHist = EmpHist,
      weights = weights, fix.a.1pl = fix.a.1pl, fix.a.gpcm = fix.a.gpcm, fix.g = fix.g,
      a.val.1pl = a.val.1pl, a.val.gpcm = a.val.gpcm, g.val = g.val,
      use.aprior = use.aprior, use.bprior = use.bprior, use.gprior = use.gprior,
      aprior = aprior, bprior = bprior, gprior = gprior, group.mean = group.mean,
      group.var = group.var, nstd = nstd, control = control,
      iter = r, fipc = FALSE, reloc.par = param_loc$reloc.par, parbd = parbd
    )
    
    # compute the difference between previous and updated item parameter estimates
    diff_par <- mstep$elm_item$pars - elm_item$pars
    # convergence statistic: largest absolute parameter change
    max.diff <- suppressWarnings(max(abs(diff_par), na.rm = TRUE))
    
    # loglikelihood value
    llike <- mstep$loglike
    
    # print
    if (verbose) {
      cat("\r", paste0(
        "EM iteration: ", r, ", Loglike: ", format(round(mstep$loglike, 4), nsmall = 4),
        ", Max-Change: ", format(round(max.diff, 6), nsmall = 5)
      ))
    }
    
    # check the convergence of EM algorithm; a non-finite change is not converged
    converge <- is.finite(max.diff) && max.diff <= Etol
    
    # extract the updated item parameter estimates
    elm_item$pars <- mstep$elm_item$pars

    # extract the updated quadrature points and
    # the corresponding weights of the prior population density
    weights <- mstep$weights
    
    # terminate the EM step if the convergence criterion is satisfied
    if (converge | r == MaxE) {
      break
    }
  }
  
  if (verbose) {
    cat("", "\n")
  }
  time2 <- Sys.time()
  
  # record the item parameter estimation time
  est_time1 <- round(as.numeric(difftime(time2, time1, units = "secs")), 2)
  
  # the first order test: check convergence-criteria test
  # require M-step convergence for every item and the EM criterion
  test_1st <- all(mstep$convergence == 0L) && converge
  if (test_1st) {
    memo3 <- "Convergence criteria are satisfied."
  } else {
    memo3 <- "Convergence criteria are not satisfied."
    warning(paste0(memo3, " \n"), call. = FALSE)
  }
  
  # conduct one more E-step to update the posterior distribution
  # using the final item parameter estimates
  estep <- Estep(
    elm_item = elm_item, idx.drm = idx.drm, idx.prm = idx.prm,
    data_drm = data_drm, data_prm = data_prm, data_all = data_all,
    weights = weights, D = D
  )
  
  # compute the final log of marginal likelihood
  # (adding back the log-scale offsets of the rows rescaled in likelihood())
  llike <- sum(log(estep$likehd %*% matrix(weights[, 2]))) + sum(estep$offset)
  
  ## ---------------------------------------------------------------------
  # estimates the information matrix and standard errors of item parameter estimates
  ## ---------------------------------------------------------------------
  # extract the finalized posterior density
  post_dist <- estep$post_dist
  
  # delete 'estep' and 'mstep' object
  rm(estep, mstep, data_drm, data_prm, data_all, envir = environment(), inherits = FALSE)
  
  # compute the information matrix of the item parameter estimates:
  # the cross-product information plus the prior information
  if (se) {
    if (verbose) {
      cat("Computing item parameter var-covariance matrix...", "\n")
    }
    time1 <- Sys.time()

    # compute the information matrix of item parameters on the quadrature grid
    info.data <- info_xpd(
      elm_item = elm_item, freq.cat = freq.cat, post_dist = post_dist,
      quadpt = quadpt, nstd = nstd,
      D = D, loc_1p_const = loc_1p_const, loc_else = loc_else, n.1PLM = n.1PLM,
      fix.a.1pl = fix.a.1pl, fix.a.gpcm = fix.a.gpcm, fix.g = fix.g, a.val.1pl = a.val.1pl,
      a.val.gpcm = a.val.gpcm, g.val = g.val, reloc.par = param_loc$reloc.par
    )

    # compute the information matrix of item parameter priors
    info.prior <- info_prior(
      elm_item = elm_item, D = D, loc_1p_const = loc_1p_const,
      loc_else = loc_else, n.1PLM = n.1PLM, fix.a.1pl = fix.a.1pl, fix.a.gpcm = fix.a.gpcm, fix.g = fix.g,
      aprior = aprior, bprior = bprior,
      gprior = gprior, use.aprior = use.aprior, use.bprior = use.bprior, use.gprior = use.gprior,
      reloc.par = param_loc$reloc.par
    )

    # sum of two information matrices
    info.mat <- info.data + info.prior

    # delete 'info.data' and 'info.prior' objects
    rm(info.data, info.prior, envir = environment(), inherits = FALSE)
    
    # second-order test + variance-covariance matrix in one Cholesky:
    # chol(info.mat) succeeds iff info.mat is positive-definite, and the
    # cached factor R lets chol2inv(R) compute the inverse.  When chol() fails
    # (rare for converged solutions), fall back to an eigenvalue check and solve().
    chol_R <- suppressWarnings(tryCatch(chol(info.mat), error = function(e) NULL))
    if (!is.null(chol_R)) {
      test_2nd <- TRUE
      cov_mat <- chol2inv(chol_R)
    } else {
      test_2nd <- all(eigen(info.mat, only.values = TRUE)$values > 1e-20)
      cov_mat <- suppressWarnings(tryCatch(
        {
          solve(info.mat, tol = 1e-200)
        },
        error = function(e) {
          NULL
        }
      ))
    }
    if (test_2nd) {
      if (test_1st) {
        memo4 <- "Solution is a possible local maximum."
      } else {
        memo4 <- "Information matrix of item parameter estimates is positive definite."
      }
    } else {
      memo4 <- "Information matrix of item parameter estimates is not positive definite; unstable solution."
      warning(paste0(memo4, " \n"), call. = FALSE)
    }

    # compute the standard errors of item parameter estimates
    if (is.null(cov_mat)) {
      se_par <- rep(99999, length(diag(info.mat)))
      memo5 <- "Variance-covariance matrix of item parameter estimates is not obtainable; unstable solution."
      warning(paste0(memo5, " \n"), call. = FALSE)
    } else {
      se_par <- suppressWarnings(sqrt(diag(cov_mat)))
      memo5 <- "Variance-covariance matrix of item parameter estimates is obtainable."
    }

    # prevent showing NaN values of standard errors
    if (any(is.nan(se_par))) {
      se_par[is.nan(se_par)] <- 99999
    }

    # set an upper bound of standard error
    se_par <- ifelse(se_par > 99999, 99999, se_par)
    time2 <- Sys.time()

    # record the standard error computation time
    est_time2 <- round(as.numeric(difftime(time2, time1, units = "secs")), 2)
  } else {
    memo4 <- "Information matrix of item parameter estimates is not computed."
    memo5 <- "Variance-covariance matrix of item parameter estimates is not computed."
    cov_mat <- NULL
    se_par <- NULL
    est_time2 <- NULL
  }
  
  # create an item metadata including the final item parameter estimates
  par_df <-
    data.frame(x[, 1:3], elm_item$pars) %>%
    confirm_df(g2na = TRUE)
  
  # deploy the standard errors into the same row/column layout that
  # holds the item parameter estimates.  loc.par is a nitem x ncols
  # matrix where each non-NA cell stores the position into the flat
  # se_par vector (i.e. the 1-based index of that estimated
  # parameter), and NA cells correspond to parameter slots that this
  # item's model does not use.
  se_df <- loc.par <- param_loc$loc.par
  mask  <- !is.na(loc.par)
  if (se) {
    se_df[mask] <- se_par[loc.par[mask]]
  } else {
    se_df[mask] <- NA_real_
  }
  
  # create a full data.frame for the standard error estimates
  se_df <-
    data.frame(x[, 1:3], se_df) %>%
    confirm_df(g2na = TRUE)
  
  # create a full data.frame containing the position of item parameter estimates
  # this data.frame is useful when interpreting the variance-covariance matrix of item parameter estimates
  loc_df <-
    data.frame(x[, 1:3], loc.par) %>%
    confirm_df(g2na = TRUE)
  
  ## ---------------------------------------------------------------------
  # summarize the estimation results
  ## ---------------------------------------------------------------------
  # create a full data.frame including both the item parameter estimates and standard error estimates
  all_df <- data.frame(matrix(NA, nrow = nrow(loc.par), ncol = 2 * ncol(loc.par)))
  all_df[, seq(1, 2 * ncol(loc.par), 2)] <- par_df[, -c(1:3)]
  all_df[, seq(2, 2 * ncol(loc.par), 2)] <- se_df[, -c(1:3)]
  col.names <- rep(NA, 2 * ncol(loc.par))
  col.names[seq(1, 2 * ncol(loc.par), 2)] <- paste0("par.", 1:ncol(loc.par))
  col.names[seq(2, 2 * ncol(loc.par), 2)] <- paste0("se.", 1:ncol(loc.par))
  colnames(all_df) <- col.names
  full_all_df <- data.frame(x[, 1:3], all_df)
  
  # population density parameters
  moments <- c(mu = group.mean, sigma2 = group.var, sigma = sqrt(group.var))
  moments.se <- rep(NA, 3)
  
  # data.frame for the population density parameter estimates
  group.par <- data.frame(rbind(moments, moments.se))
  colnames(group.par) <- c("mu", "sigma2", "sigma")
  rownames(group.par) <- c("estimates", "se")
  
  # prior information
  if (use.aprior) aprior.dist <- aprior else aprior.dist <- NULL
  if (use.bprior) bprior.dist <- bprior else bprior.dist <- NULL
  if (use.gprior) gprior.dist <- gprior else gprior.dist <- NULL
  
  # count the free parameters: the item parameters, and the weights of the empirical
  # histogram, whose number of free values is the number of quadrature points minus
  # three (sum to one, mean, and variance are fixed); the normal distribution is fixed
  npar.item <- length(param_loc$reloc.par)
  npar.dist <- 0L
  npar.weight <- if (EmpHist) n.quad - 3L else 0L
  npar.est <- npar.item + npar.dist + npar.weight
  npar.detail <- c(item = npar.item, dist = npar.dist, weight = npar.weight)

  # statistics based on the loglikelihood of the fitted model:
  neg2llke <- -2 * llike
  aic <- 2 * npar.est + neg2llke
  bic <- npar.est * log(nstd) + neg2llke
  
  ## ---------------------------------------------------------------
  # check end time
  end.time <- Sys.time()
  
  # record total computation time
  est_time3 <- round(as.numeric(difftime(end.time, start.time, units = "secs")), 2)
  
  # return results
  rst <- list(
    estimates = full_all_df, par.est = par_df, se.est = se_df, pos.par = loc_df, covariance = cov_mat, loglikelihood = llike, aic = aic, bic = bic,
    group.par = group.par, weights = weights, posterior.dist = post_dist, data = data, scale.D = D, ncase = nstd, nitem = nrow(par_df),
    Etol = Etol, MaxE = MaxE, aprior = aprior.dist, bprior = bprior.dist, gprior = gprior.dist, npar.est = npar.est, niter = r, maxpar.diff = max.diff,
    EMtime = est_time1, SEtime = est_time2, TotalTime = est_time3, test.1 = memo3, test.2 = memo4, var.note = memo5, fipc = FALSE,
    fipc.method = NULL, fix.loc = NULL, npar.detail = npar.detail
  )
  
  if (verbose) {
    cat("Estimation is finished in", est_time3, "seconds.", "\n")
  }
  return(rst)
}



# This function implement the fixed item parameter calibration (FIPC) using MMLE-EM algorithm
#' @import dplyr
est_irt_fipc <- function(x = NULL,
                         data,
                         D = 1,
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
                         aprior = list(dist = "lnorm", params = c(0.0, 0.5)),
                         bprior = list(dist = "norm", params = c(0.0, 1.0)),
                         gprior = list(dist = "beta", params = c(5, 16)),
                         missing = NA,
                         Quadrature = c(49, 6.0),
                         weights = NULL,
                         group.mean = 0.0,
                         group.var = 1.0,
                         EmpHist = FALSE,
                         use.startval = FALSE,
                         Etol = 1e-04,
                         MaxE = 500,
                         control = list(eval.max = 500, iter.max = 200, x.tol = 1e-4),
                         fipc = TRUE,
                         fipc.method = "MEM",
                         fix.loc = NULL,
                         fix.id = NULL,
                         se = TRUE,
                         verbose = TRUE) {
  
  # check start time
  start.time <- Sys.time()
  
  ## ---------------------------------------------------------------------
  # prepare the item parameter estimation
  ## ---------------------------------------------------------------------
  # check if item metadata argument of 'x' is provided
  if (is.null(x)) {
    stop(paste0(
      "To implement the fixed item parameter calibration, \n",
      "the item metadata must be specified in the argument 'x'."
    ), call. = FALSE)
  }
  
  # recode a non-numeric missing code before the responses are read as numbers
  if (!is.na(missing)) data[data == missing] <- NA

  # transform a data set to matrix (character and factor columns are read by their values)
  data <- resp_as_numeric(data)
  
  # extract information about the number of score categories and models
  if (verbose) {
    cat("Parsing input...", "\n")
  }
  
  # confirm and correct all item metadata information
  x <- confirm_df(x)
  
  # override the item id in x if item.id is not NULL
  if (!is.null(item.id)) {
    # stop when item.id does not have one value per item of x
    if (length(item.id) != nrow(x)) {
      stop("The length of 'item.id' must be equal to the number of items in 'x'.", call. = FALSE)
    }
    x$id <- item.id
  }

  # check the number of items
  nitem.all <- nrow(x)

  # check that 'x' and 'data' contain the same number of items
  if (nitem.all != ncol(data)) {
    stop("The number of items included in 'x' and 'data' must be the same.", call. = FALSE)
  }

  # find the items that should be fixed
  if (is.null(fix.loc) & is.null(fix.id)) {
    stop("When 'FIPC = TRUE', the information of which items are fixed \n",
         "must be provided via the 'fix.loc' or 'fix.id' argument.",
         call. = FALSE
    )
  }
  if (!is.null(fix.loc) & !is.null(fix.id)) {
    warning(paste0(
      "The information given to the 'fix.id' argument was used to fix the item parameters and \n",
      "the information given to the 'fix.loc' argument was ignored."
    ), call. = FALSE)
  }
  if (!is.null(fix.id)) {
    # stop when an ID in 'fix.id' is not found in the item IDs
    id.notfound <- setdiff(unique(fix.id), x$id)
    if (length(id.notfound) > 0L) {
      stop("The item ID(s) in 'fix.id' not found in the item IDs: ", paste(id.notfound, collapse = ", "), ".",
           call. = FALSE)
    }
    fix.loc <- which(x$id %in% unique(fix.id))
  }

  # stop when 'fix.loc' is empty or points outside the rows of 'x'
  if (length(fix.loc) == 0L || anyNA(fix.loc) || any(!fix.loc %in% seq_len(nitem.all))) {
    stop("'fix.loc' must contain row positions of 'x' between 1 and ", nitem.all, ".", call. = FALSE)
  }

  # stop when an item is listed more than once (it would be counted twice in the first E-step)
  if (anyDuplicated(fix.loc) > 0L) {
    stop("'fix.loc' must not contain duplicated row positions.", call. = FALSE)
  }

  # check the location of items whose item parameters are estimated
  nofix.loc <- c(1:nitem.all)[!c(1:nitem.all) %in% fix.loc]
  if (length(nofix.loc) == 0L) {
    nofix.loc <- NULL
  }
  
  # divide the item metadata into two groups: fixed (x_fix) and non-fixed (x_new)
  x_fix <- x[fix.loc, ]
  nitem.fix <- nrow(x_fix)
  if (!is.null(nofix.loc)) {
    x_new <- x[nofix.loc, ]
    nitem.new <- nrow(x_new)
  } else {
    x_new <- NULL
  }
  
  # clear the two item metadata sets
  x_fix <- confirm_df(x_fix)
  if (!is.null(nofix.loc)) {
    x_new <- confirm_df(x_new)
  }
  
  # record the score categories and model information of the new items to be estimated
  if (!is.null(x_new)) {
    id <- x_new$id
    cats <- x_new$cats
    model <- x_new$model
  } else {
    id <- x_fix$id
    cats <- x_fix$cats
    model <- x_fix$model
  }
  
  # generate the empty metadata with starting values
  if (!is.null(x_new) && !use.startval) {
    x_new <- startval_df(cats = cats, model = model, item.id = id)
  }
  
  # create the total item metadata to be used in the further estimation process
  x_all <- startval_df(cats = x$cats, model = x$model, item.id = x$id)
  x_all[fix.loc, 1:ncol(x_fix)] <- x_fix
  if (!is.null(nofix.loc)) {
    x_all[nofix.loc, 1:ncol(x_new)] <- x_new
  }
  
  # recode missing values
  if (!is.na(missing)) {
    data[data == missing] <- NA
  }
  
  # check the number of item responses across all items
  n.resp <- Rfast::colsums(!is.na(data))
  
  # check the items which have all missing responses
  loc_allmiss <- which(n.resp == 0L)
  if (length(loc_allmiss) > 0L) {
    memo2 <- paste0(paste0("item ", loc_allmiss, collapse = ", "), " has/have no item response data. \n")
    stop(memo2, call. = FALSE)
  }

  # stop when a response is not one of the score categories 0, ..., cats - 1
  check_resp_range(data = data, cats = x$cats, id = x$id)

  # stop when the responses to all fixed items show no variation
  check_fix_info(data = data[, fix.loc, drop = FALSE], id = x$id[fix.loc])

  # warn about the items to be estimated whose maximum likelihood estimates do not exist
  if (!is.null(nofix.loc)) {
    check_item_mle(data = data[, nofix.loc, drop = FALSE], cats = x$cats[nofix.loc], id = x$id[nofix.loc])
  }

  # save the item response data into the different object
  data_all <- data
  rm(data, envir = environment(), inherits = FALSE)
  
  # divide the item response data into two groups: fixed (data_fix) and non-fixed (data_new)
  data_fix <- data_all[, fix.loc, drop = FALSE]
  if (!is.null(x_new)) {
    data_new <- data_all[, nofix.loc, drop = FALSE]
  } else {
    data_new <- NULL
  }
  
  # find the location of 1PLM items in which slope parameters should be constrained to be equal
  # also, find the location of other items
  if ("1PLM" %in% model & !fix.a.1pl) {
    loc_1p_const <- which(model == "1PLM")
    loc_else <- which(model != "1PLM")
    
    # count the number of 1PLM items to be constrained
    n.1PLM <- length(loc_1p_const)
  } else {
    loc_1p_const <- NULL
    n.1PLM <- NULL
    if (!is.null(x_new)) {
      loc_else <- 1:nrow(x_new)
    } else {
      loc_else <- 1:nrow(x_all)
    }
  }
  
  # record the original location of item parameters to be estimated, and
  # the relocated position of item parameters when computing
  # the variance-covariance matrix of item parameter estimates
  if (!is.null(x_new)) {
    param_loc <- parloc(
      x = x_new, loc_1p_const = loc_1p_const, loc_else = loc_else,
      fix.a.1pl = fix.a.1pl, fix.a.gpcm = fix.a.gpcm, fix.g = fix.g
    )
  } else {
    param_loc <- parloc(
      x = x_all, loc_1p_const = loc_1p_const, loc_else = loc_else,
      fix.a.1pl = FALSE, fix.a.gpcm = FALSE, fix.g = FALSE
    )
    param_loc$loc.par[!is.na(param_loc$loc.par)] <- NA_real_
    param_loc$reloc.par <- NULL
  }
  
  
  ## ---------------------------------------------------------------------
  # conduct item parameter estimation using MMLE-EM algorithm
  ## ---------------------------------------------------------------------
  # check the total number of examinees
  nstd <- nrow(data_all)
  
  # create initial weights of prior ability distribution when it is not specified
  if (is.null(weights)) {
    # create quadrature points
    quadpt <- seq(-Quadrature[2], Quadrature[2], length.out = Quadrature[1])
    n.quad <- length(quadpt)
    
    # create the data.frame containing the quadrature points and weights
    weights <- gen.weight(dist = "norm", mu = group.mean, sigma = sqrt(group.var), theta = quadpt)
  } else {
    # check the weights and rescale them to sum to one, as the first EM cycle uses them
    weights <- check_weights(weights)
    quadpt <- weights[, 1]
    n.quad <- length(quadpt)

    # use the moments of the supplied weights as the initial moments of the distribution
    moments.tmp <- cal_moment(node = quadpt, weight = weights[, 2])
    group.mean <- moments.tmp[1]
    group.var <- moments.tmp[2]
  }

  # build the per-item one-hot frequency-category list once on the combined
  # response matrix; data_fix and data_new are the fix.loc and nofix.loc
  # columns of data_all, so list subsetting gives their freq.cat without
  # copying the underlying matrices
  freq_all.cat <- build_freqcat(data_all, x_all$cats)
  if (!is.null(x_new)) {
    freq_fix.cat <- freq_all.cat[fix.loc]
    freq_new.cat <- freq_all.cat[nofix.loc]
  } else {
    freq_fix.cat <- NULL
    freq_new.cat <- NULL
  }
  
  # break down the item metadata into several elements
  if (!is.null(x_new)) {
    elm_item_new <- breakdown(x_new)
    elm_item_fix <- breakdown(x_fix)
  }
  elm_item_all <- breakdown(x_all)
  
  # classify the items into DRM and PRM item groups
  if (!is.null(x_new)) {
    idx.item.new <- idxfinder(elm_item_new)
    idx.item.fix <- idxfinder(elm_item_fix)
    idx.drm.new <- idx.item.new$idx.drm
    idx.prm.new <- idx.item.new$idx.prm
    idx.drm.fix <- idx.item.fix$idx.drm
    idx.prm.fix <- idx.item.fix$idx.prm
  }
  idx.item.all <- idxfinder(elm_item_all)
  idx.drm.all <- idx.item.all$idx.drm
  idx.prm.all <- idx.item.all$idx.prm
  
  # then, divide items in the new test form into three groups
  # (1) DRM 1PL items with the constrained slope
  # (2) other DRM items
  # (3) PRM items
  # (no item other than the constrained 1PLM items when loc_else is empty)
  if (!is.null(x_new)) {
    if (length(loc_else) == 0) {
      drm.else <- NULL
    } else {
      drm.else <- loc_else
    }
    idx4est <- list(
      drm.slc = loc_1p_const,
      drm.else = drm.else,
      prm = idx.prm.new
    )
  } else {
    idx4est <- NULL
  }
  
  # divide the data set for the mixed-item format
  if (!is.null(x_new)) {
    datlist.new <- divide_data(data = data_new, idx.item = idx.item.new, freq.cat = freq_new.cat)
    datlist.fix <- divide_data(data = data_fix, idx.item = idx.item.fix, freq.cat = freq_fix.cat)
    data_drm.new <- cbind(datlist.new$data_drm_q, datlist.new$data_drm_p)
    data_drm.fix <- cbind(datlist.fix$data_drm_q, datlist.fix$data_drm_p)
    data_prm.new <- datlist.new$data_prm
    data_prm.fix <- datlist.fix$data_prm
    data.new <- datlist.new$data_all
    data.fix <- datlist.fix$data_all
  }
  datlist.all <- divide_data(data = data_all, idx.item = idx.item.all, freq.cat = freq_all.cat)
  data_drm.all <- cbind(datlist.all$data_drm_q, datlist.all$data_drm_p)
  data_prm.all <- datlist.all$data_prm
  data.all <- datlist.all$data_all
  
  # delete 'datlist' object
  if (!is.null(x_new)) {
    rm(datlist.new, datlist.fix, freq_fix.cat, envir = environment(), inherits = FALSE)
  }
  rm(datlist.all, freq_all.cat, envir = environment(), inherits = FALSE)
  
  # moments of population prior distribution
  if (is.null(x_new)) {
    mmt_dist_old <- c(group.mean, group.var)
  }

  # create the lower and upper bounds of the item parameters
  if (!is.null(x_new)) {
    parbd <- lubound(model, cats, n.1PLM, idx4est, fix.a.1pl, fix.g, fix.a.gpcm)
  } else {
    parbd <- NULL
  }
  
  # find the columns of the frequency matrix corresponding to all items
  if (!is.null(x_new)) {
    cols.item <- cols4item(nitem.new, cats, loc_1p_const)
  } else {
    cols.item <- NULL
  }
  
  # estimation
  if (verbose) {
    cat("Estimating item parameters...", "\n")
  }
  
  # set the number of EM iteration to one when OEM (Wainer & Mislevy, 1990) method is used
  if (fipc.method == "OEM") MaxE <- 1
  
  # implement EM algorithm
  time1 <- Sys.time()
  for (r in 1:MaxE) {
    # implement E-step
    if (!is.null(x_new)) {
      if (r == 1L) {
        estep <- Estep_fipc(
          elm_item1 = elm_item_new, elm_item2 = elm_item_fix, idx.drm2 = idx.drm.fix,
          idx.prm2 = idx.prm.fix, data_drm2 = data_drm.fix, data_prm2 = data_prm.fix,
          data_all1 = data.new, weights = weights, D = D
        )
      } else {
        estep <- Estep_fipc(
          elm_item1 = elm_item_new, elm_item2 = elm_item_all, idx.drm2 = idx.drm.all,
          idx.prm2 = idx.prm.all, data_drm2 = data_drm.all, data_prm2 = data_prm.all,
          data_all1 = data.new, weights = weights, D = D
        )
      }
    } else {
      estep <- Estep_fipc(
        elm_item1 = elm_item_all, elm_item2 = elm_item_all, idx.drm2 = idx.drm.all,
        idx.prm2 = idx.prm.all, data_drm2 = data_drm.all, data_prm2 = data_prm.all,
        data_all1 = data.all, weights = weights, D = D
      )
      estep$elm_item <- NULL
    }
    
    # implement M-step
    mstep <- Mstep(
      estep = estep, id = id, cats = cats, model = model, quadpt = quadpt, n.quad = n.quad,
      D = D, cols.item = cols.item, loc_1p_const = loc_1p_const, loc_else = loc_else, idx4est = idx4est,
      n.1PLM = n.1PLM, EmpHist = EmpHist, weights = weights, fix.a.1pl = fix.a.1pl, fix.a.gpcm = fix.a.gpcm, fix.g = fix.g,
      a.val.1pl = a.val.1pl, a.val.gpcm = a.val.gpcm, g.val = g.val, use.aprior = use.aprior, use.bprior = use.bprior,
      use.gprior = use.gprior, aprior = aprior, bprior = bprior, gprior = gprior, group.mean = group.mean,
      group.var = group.var, nstd = nstd, control = control,
      iter = r, fipc = TRUE, reloc.par = param_loc$reloc.par, parbd = parbd
    )
    
    if (!is.null(x_new)) {
      # compute the difference between previous and updated item parameter estimates
      diff_par <- mstep$elm_item$pars - elm_item_new$pars
      # convergence statistic: largest absolute parameter change
      max.diff <- suppressWarnings(max(abs(diff_par), na.rm = TRUE))
    } else {
      # compute the mean and variance of the updated prior distribution
      mmt_dist_new <- cal_moment(node = mstep$weights$theta, weight = mstep$weights$weight)
      diff_par <- mmt_dist_new - mmt_dist_old
      # convergence statistic: largest absolute change in prior mean and variance
      max.diff <- suppressWarnings(max(abs(diff_par), na.rm = TRUE))
    }

    # log-likelihood value
    llike <- mstep$loglike
    
    # print
    if (verbose) {
      cat("\r", paste0(
        "EM iteration: ", r, ", Loglike: ", format(round(mstep$loglike, 4), nsmall = 4),
        ", Max-Change: ", format(round(max.diff, 6), nsmall = 5)
      ))
    }
    
    # check the convergence of EM algorithm; a non-finite change is not converged
    converge <- is.finite(max.diff) && max.diff <= Etol
    
    # extract the updated item (or group) parameter estimates
    # and update the new and all item parameters
    if (!is.null(x_new)) {
      elm_item_new$pars <- mstep$elm_item$pars
      elm_item_all$pars[nofix.loc, 1:ncol(elm_item_new$pars)] <- elm_item_new$pars
    } else {
      mmt_dist_old <- mmt_dist_new
    }

    # extract the updated quadrature points and
    # the corresponding weights of the prior population density
    weights <- mstep$weights
    
    # terminate the EM step if the convergence criterion is satisfied
    if (converge | r == MaxE) {
      break
    }
  }
  
  if (verbose) {
    cat("", "\n")
  }
  time2 <- Sys.time()
  
  # record the item parameter estimation time
  est_time1 <- round(as.numeric(difftime(time2, time1, units = "secs")), 2)
  
  # the first order test: check convergence-criteria test
  # OEM runs a single EM cycle by design, so only the M-step check applies;
  # with all items fixed (no new items), the EM criterion must still be met
  test_1st <- all(mstep$convergence == 0L) &&
    (converge || (fipc.method == "OEM" && !is.null(x_new)))
  if (test_1st) {
    memo3 <- "Convergence criteria are satisfied."
  } else {
    memo3 <- "Convergence criteria are not satisfied."
    warning(paste0(memo3, " \n"), call. = FALSE)
  }
  
  # conduct one more E-step to update the posterior distribution
  # using the final item parameter estimates
  if (!is.null(x_new)) {
    estep <- Estep_fipc(
      elm_item1 = elm_item_new, elm_item2 = elm_item_all, idx.drm2 = idx.drm.all,
      idx.prm2 = idx.prm.all, data_drm2 = data_drm.all, data_prm2 = data_prm.all,
      data_all1 = data.new, weights = weights, D = D
    )
  } else {
    estep <- Estep_fipc(
      elm_item1 = elm_item_all, elm_item2 = elm_item_all, idx.drm2 = idx.drm.all,
      idx.prm2 = idx.prm.all, data_drm2 = data_drm.all, data_prm2 = data_prm.all,
      data_all1 = data.all, weights = weights, D = D
    )
  }
  
  # compute the final log of marginal likelihood
  # (adding back the log-scale offsets of the rows rescaled in likelihood())
  llike <- sum(log(estep$likehd %*% matrix(weights[, 2]))) + sum(estep$offset)
  
  # compute the mean and variance of the estimated density distribution
  pop_moments <- cal_moment(node = quadpt, weight = weights[, 2])
  
  ## ---------------------------------------------------------------------
  # estimates the standard errors of item parameter estimates
  ## ---------------------------------------------------------------------
  if (!is.null(x_new)) {
    # extract the finalized posterior density
    post_dist <- estep$post_dist
    
    # delete 'estep' and 'mstep' object
    rm(estep, mstep, data_drm.new, data_drm.fix,
       data_prm.new, data_prm.fix, data.new, data.fix,
       data_drm.all, data_prm.all, data.all,
       envir = environment(), inherits = FALSE
    )
    
    # compute the information matrix of the item parameter estimates:
    # the cross-product information plus the prior information
    if (se) {
      if (verbose) {
        cat("Computing item parameter var-covariance matrix...", "\n")
      }
      time1 <- Sys.time()

      # compute the information matrix of item parameters
      info.data <- info_xpd(
        elm_item = elm_item_new, freq.cat = freq_new.cat, post_dist = post_dist,
        quadpt = quadpt, nstd = nstd,
        D = D, loc_1p_const = loc_1p_const, loc_else = loc_else, n.1PLM = n.1PLM,
        fix.a.1pl = fix.a.1pl, fix.a.gpcm = fix.a.gpcm, fix.g = fix.g, a.val.1pl = a.val.1pl,
        a.val.gpcm = a.val.gpcm, g.val = g.val, reloc.par = param_loc$reloc.par
      )
      
      # compute the information matrix of item parameter priors
      info.prior <- info_prior(
        elm_item = elm_item_new, D = D, loc_1p_const = loc_1p_const,
        loc_else = loc_else, n.1PLM = n.1PLM, fix.a.1pl = fix.a.1pl, fix.a.gpcm = fix.a.gpcm, fix.g = fix.g,
        aprior = aprior, bprior = bprior,
        gprior = gprior, use.aprior = use.aprior, use.bprior = use.bprior, use.gprior = use.gprior,
        reloc.par = param_loc$reloc.par
      )
      
      # sum of two information matrices
      info.mat <- info.data + info.prior
      
      # second-order test + variance-covariance matrix in one Cholesky:
      # see est_irt_em()
      chol_R <- suppressWarnings(tryCatch(chol(info.mat), error = function(e) NULL))
      if (!is.null(chol_R)) {
        test_2nd <- TRUE
        cov_mat <- chol2inv(chol_R)
      } else {
        test_2nd <- all(eigen(info.mat, only.values = TRUE)$values > 1e-20)
        cov_mat <- suppressWarnings(tryCatch(
          {
            solve(info.mat, tol = 1e-200)
          },
          error = function(e) {
            NULL
          }
        ))
      }
      if (test_2nd) {
        if (test_1st) {
          memo4 <- "Solution is a possible local maximum."
        } else {
          memo4 <- "Information matrix of item parameter estimates is positive definite."
        }
      } else {
        memo4 <- "Information matrix of item parameter estimates is not positive definite; unstable solution."
        warning(paste0(memo4, " \n"), call. = FALSE)
      }

      # compute the standard errors of item parameter estimates
      if (is.null(cov_mat)) {
        se_par <- rep(99999, length(diag(info.mat)))
        memo5 <- "Variance-covariance matrix of item parameter estimates is not obtainable; unstable solution."
        warning(paste0(memo5, " \n"), call. = FALSE)
      } else {
        se_par <- suppressWarnings(sqrt(diag(cov_mat)))
        memo5 <- "Variance-covariance matrix of item parameter estimates is obtainable."
      }
      
      # prevent showing NaN values of standard errors
      if (any(is.nan(se_par))) {
        se_par[is.nan(se_par)] <- 99999
      }
      
      # set an upper bound of standard error
      se_par <- ifelse(se_par > 99999, 99999, se_par)
      time2 <- Sys.time()
      
      # record the standard error computation time
      est_time2 <- round(as.numeric(difftime(time2, time1, units = "secs")), 2)
    } else {
      memo4 <- "Information matrix of item parameter estimates is not computed."
      memo5 <- "Variance-covariance matrix of item parameter estimates is not computed."
      cov_mat <- NULL
      se_par <- NULL
      est_time2 <- NULL
    }
    
    # create an item metadata including the final item parameter estimates for all items
    x_all <-
      data.frame(x_all[, 1:3], elm_item_all$pars) %>%
      confirm_df(g2na = TRUE)
    
    # deploy the standard errors into the same row/column layout that
    # holds the item parameter estimates.
    # 1) for the new items only (see est_irt_em() for the loc.par layout).
    se_df <- loc.par <- param_loc$loc.par
    mask  <- !is.na(loc.par)
    if (se) {
      se_df[mask] <- se_par[loc.par[mask]]
    } else {
      se_df[mask] <- NA_real_
    }
    
    # 2) for the total test form
    all.col <- max(x_all$cats)
    all.col <- ifelse(all.col == 2, 3, all.col)
    se_all_df <- loc_all.par <- matrix(NA, nrow = nitem.all, ncol = all.col)
    if (ncol(se_df) < ncol(se_all_df)) {
      n2add <- ncol(se_all_df) - ncol(se_df)
      se_df <- cbind(se_df, matrix(NA, nrow = nrow(se_df), ncol = n2add))
    }
    se_all_df[nofix.loc, 1:ncol(se_df)] <- se_df
    loc_all.par[nofix.loc, 1:ncol(loc.par)] <- loc.par
    
    # create a full data.frame for the standard error estimates
    se_all_df <-
      data.frame(x_all[, 1:3], se_all_df) %>%
      confirm_df(g2na = TRUE)
    
    # create a full data.frame containing the position of item parameter estimates
    # this data.frame is useful when interpreting the variance-covariance matrix of item parameter estimates
    loc_all_df <-
      data.frame(x_all[, 1:3], loc_all.par) %>%
      confirm_df(g2na = TRUE)
  } else {
    # extract the finalized posterior density
    post_dist <- estep$post_dist
    
    # the second-order test is not available because the information matrix is not computed
    memo4 <- "No item parameters are estimated, so the information matrix is not computed."

    # compute the standard errors of item parameter estimates
    cov_mat <- NULL
    memo5 <- "Variance-covariance matrix of item parameter estimates was not estimated."
    
    # record the standard error computation time (NULL)
    est_time2 <- NULL
    
    # all standard errors are NA because no item parameter is estimated
    se_df <- loc_all.par <- param_loc$loc.par
    
    # create a full data.frame for the standard error estimates
    # however, SEs are all NA when only group parameters are estimated
    se_all_df <-
      data.frame(x_all[, 1:3], se_df) %>%
      confirm_df(g2na = TRUE)
    
    # create a full data.frame containing the position of item parameter estimates
    # however, this should be NULL when only group parameters are estimated
    loc_all_df <- NULL
  }
  
  ## ---------------------------------------------------------------------
  # summarize the estimation results
  ## ---------------------------------------------------------------------
  # create a full data.frame including both the item parameter estimates and standard error estimates
  all_df <- data.frame(matrix(NA, nrow = nrow(loc_all.par), ncol = 2 * ncol(loc_all.par)))
  all_df[, seq(1, 2 * ncol(loc_all.par), 2)] <- x_all[, -c(1:3)]
  all_df[, seq(2, 2 * ncol(loc_all.par), 2)] <- se_all_df[, -c(1:3)]
  col.names <- rep(NA, 2 * ncol(loc_all.par))
  col.names[seq(1, 2 * ncol(loc_all.par), 2)] <- paste0("par.", 1:ncol(loc_all.par))
  col.names[seq(2, 2 * ncol(loc_all.par), 2)] <- paste0("se.", 1:ncol(loc_all.par))
  colnames(all_df) <- col.names
  full_all_df <- data.frame(x_all[, 1:3], all_df)
  
  # population density parameters
  mu <- pop_moments[1]
  sigma2 <- pop_moments[2]
  sigma <- sqrt(pop_moments[2])
  moments.est <- c(mu, sigma2, sigma)
  
  # compute the standard errors of population density parameters
  se.mu <- sigma / sqrt(nstd)
  se.sigma2 <- sigma2 * sqrt(2 / (nstd - 1))
  se.sigma <- (1 / (2 * sigma)) * se.sigma2 # using Delta method
  moments.se <- c(se.mu, se.sigma2, se.sigma)
  
  # data.frame for the population density parameter estimates
  group.par <- data.frame(rbind(moments.est, moments.se))
  colnames(group.par) <- c("mu", "sigma2", "sigma")
  rownames(group.par) <- c("estimates", "se")
  
  # prior information
  if (use.aprior) aprior.dist <- aprior else aprior.dist <- NULL
  if (use.bprior) bprior.dist <- bprior else bprior.dist <- NULL
  if (use.gprior) gprior.dist <- gprior else gprior.dist <- NULL
  
  # count the free parameters: the item parameters to be estimated, and the latent
  # distribution, which is estimated without rescaling (the mean and variance of the
  # normal distribution, or the quadrature points minus one weights of the empirical
  # histogram, whose sum is one)
  npar.item <- if (!is.null(x_new)) length(param_loc$reloc.par) else 0L
  npar.dist <- if (EmpHist) 0L else 2L
  npar.weight <- if (EmpHist) n.quad - 1L else 0L
  npar.est <- npar.item + npar.dist + npar.weight
  npar.detail <- c(item = npar.item, dist = npar.dist, weight = npar.weight)

  # statistics based on the loglikelihood of the fitted model:
  neg2llke <- -2 * llike
  aic <- 2 * npar.est + neg2llke
  bic <- npar.est * log(nstd) + neg2llke
  
  ## ---------------------------------------------------------------
  # check end time
  end.time <- Sys.time()
  
  # record total computation time
  est_time3 <- round(as.numeric(difftime(end.time, start.time, units = "secs")), 2)
  
  # return results
  rst <- list(
    estimates = full_all_df, par.est = x_all, se.est = se_all_df, pos.par = loc_all_df, covariance = cov_mat, loglikelihood = llike,
    aic = aic, bic = bic, group.par = group.par, weights = weights, posterior.dist = post_dist, data = data_all, scale.D = D, ncase = nstd,
    nitem = nitem.all, Etol = Etol, MaxE = MaxE, aprior = aprior.dist, bprior = bprior.dist, gprior = gprior.dist, npar.est = npar.est,
    niter = r, maxpar.diff = max.diff, EMtime = est_time1, SEtime = est_time2, TotalTime = est_time3, test.1 = memo3, test.2 = memo4,
    var.note = memo5, fipc = TRUE, fipc.method = fipc.method, fix.loc = fix.loc, npar.detail = npar.detail
  )

  if (verbose) {
    cat("Estimation is finished in", est_time3, "seconds.", "\n")
  }
  return(rst)
}


# This function converts a response data set to a numeric matrix by the values of the
# responses: character and factor columns are read as numbers, and a value that is not
# a number becomes -1 so that the range check flags it
resp_as_numeric <- function(data) {
  # read the values of a vector as numbers
  read_num <- function(v) {
    # numeric and logical values are used as they are
    if (is.numeric(v) || is.logical(v)) {
      return(v)
    }
    z <- suppressWarnings(as.numeric(as.character(v)))
    # a value that is not a number is flagged as out of range
    z[is.na(z) & !is.na(v)] <- -1
    z
  }

  if (is.data.frame(data)) {
    # read each column by its values
    data[] <- lapply(data, read_num)
  } else if (!is.numeric(data) && !is.logical(data)) {
    # read a character matrix by its values
    data <- matrix(read_num(as.vector(data)), nrow = nrow(data), dimnames = dimnames(data))
  }
  data.matrix(data)
}


# This function stops when a response in 'data' is not one of the score
# categories 0, 1, ..., cats - 1 of its item (missing responses are allowed)
check_resp_range <- function(data, cats, id) {
  # largest admissible score of each item, laid out like 'data'
  max.score <- matrix(cats - 1, nrow = nrow(data), ncol = ncol(data), byrow = TRUE)

  # flag negative, too large, or non-integer responses
  bad <- (data < 0) | (data > max.score) | (data != round(data))
  loc_bad <- which(Rfast::colsums(bad & !is.na(bad)) > 0L)
  if (length(loc_bad) > 0L) {
    stop(paste0(
      "The responses must be scored as 0, 1, ..., (cats - 1). ",
      "Out-of-range or non-integer responses were found for item(s) ",
      paste0(id[loc_bad], " (allowed 0 to ", cats[loc_bad] - 1, ")", collapse = ", "), "."
    ), call. = FALSE)
  }
  invisible(NULL)
}


# This function stops when the responses to every fixed item of FIPC show no variation
# (one observed response value per item). Such fixed items carry no information about
# differences in ability: the first E-step gives every examinee the same posterior, so the
# new items have no interior maximum in the M-step, and the scale cannot be placed
check_fix_info <- function(data, id) {
  # number of distinct observed responses of each fixed item
  n.value <- vapply(seq_len(ncol(data)), function(j) {
    length(unique(data[!is.na(data[, j]), j]))
  }, integer(1))

  # stop when no fixed item has two or more observed response values
  if (all(n.value <= 1L)) {
    stop(paste0(
      "The responses to all fixed items (", paste(id, collapse = ", "), ") show no variation, ",
      "so the fixed items carry no information to place the ability scale and FIPC cannot be applied."
    ), call. = FALSE)
  }
  invisible(NULL)
}


# This function warns, once, about the items that have a score category with no
# responses (which includes the items with all responses the same): their
# maximum likelihood estimates do not exist
check_item_mle <- function(data, cats, id) {
  # flag the items with at least one score category that nobody chose
  has.empty <- vapply(seq_len(ncol(data)), function(j) {
    resp <- data[, j]
    any(tabulate(resp[!is.na(resp)] + 1L, nbins = cats[j]) == 0L)
  }, logical(1))

  if (any(has.empty)) {
    warning(paste0(
      "Item(s) ", paste(id[has.empty], collapse = ", "),
      " have a score category with no responses or no variation in responses, ",
      "so their parameter estimates may be unstable and their standard errors unreliable."
    ), call. = FALSE)
  }
  invisible(NULL)
}


# This function checks the quadrature weights supplied by the user (second
# column of 'weights') and rescales them to sum to one with a warning
check_weights <- function(weights) {
  # extract the weights as a numeric vector
  wts <- as.numeric(unlist(weights[, 2]))

  # stop when the weights are not nonnegative finite values with a positive sum
  if (!all(is.finite(wts)) || any(wts < 0) || sum(wts) <= 0) {
    stop("The weights must be nonnegative finite values with a positive sum.", call. = FALSE)
  }

  # rescale the weights when they do not sum to one
  if (abs(sum(wts) - 1) > 1e-8) {
    warning("The weights do not sum to 1 and were rescaled to sum to 1.", call. = FALSE)
    weights[, 2] <- wts / sum(wts)
  }

  weights
}
