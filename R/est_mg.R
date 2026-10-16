#' Multiple-group item calibration using MMLE-EM algorithm
#'
#' This function performs multiple-group item calibration (Bock & Zimowski,
#' 1997) using marginal maximum likelihood estimation via the
#' expectation-maximization (MMLE-EM) algorithm (Bock & Aitkin, 1981). When
#' prior distributions are specified for item parameters, the estimates are the
#' modes of the marginal posterior distribution (Bayes modal estimates). It also
#' supports multiple-group fixed item parameter calibration (MG-FIPC; e.g., Kim
#' & Kolen, 2016), which extends the single-group FIPC method (Kim, 2006) to
#' multiple-group settings.
#' For dichotomous items, the function supports one-, two-, and three-parameter
#' logistic IRT models. For polytomous items, the graded response model (GRM)
#' and the (generalized) partial credit model (GPCM) are available.
#'
#' @inheritParams est_irt
#' @param x A list containing item metadata for all groups to be analyzed.
#'   For example, if five groups are analyzed, the list should contain five
#'   elements, each representing the item metadata for one group. The order of
#'   the elements in the list must match the order of group names specified in
#'   the `group.name` argument.
#'
#'   Each group's item metadata includes essential information for each item
#'   (e.g., number of score categories, IRT model type, etc.) required for
#'   calibration. See [irtQ::est_irt()] or [irtQ::simdat()] for more
#'   details about the item metadata.
#'
#'   When `use.startval = TRUE`, the item parameters specified in the metadata
#'   will be used as starting values for parameter estimation; otherwise they
#'   are ignored, except for the fixed items in MG-FIPC, whose parameters are
#'   held at the values in `x`. A common item must have the same model and
#'   number of score categories in every group that includes it. If
#'   `x = NULL`, the `model`, `cats`, and `item.id` arguments are used. Note
#'   that when `fipc = TRUE` to implement MG-FIPC, the `x` argument must be
#'   specified and cannot be NULL. Default is `NULL`.
#' @param data A list containing the item response data of all groups to be
#'   analyzed. For example, if five groups are analyzed, the list should include
#'   five elements, each representing the response data of one group. The
#'   elements in the list must be ordered to match the group names specified
#'   in the `group.name` argument. Each element is a matrix or data frame whose
#'   rows represent the examinees of the group and whose columns represent the
#'   items of the group's form, in the same order as in the group's item
#'   metadata. Responses to an item with K score categories must be integers
#'   from 0 to K - 1, and missing responses are coded as `NA` or as the value
#'   given in `missing`. Items are matched across groups by their IDs, and the
#'   items that are not in a group's form are treated as missing for the
#'   examinees of that group.
#' @param group.name A character vector indicating the names of the groups, in
#'   the order of the elements of `data`. For example, if five groups are
#'   analyzed, use `group.name = c("G1", "G2", "G3", "G4", "G5")`. Group names
#'   can be any valid character strings and label the group-level elements of
#'   the returned object. If `NULL`, the names `"g1"`, `"g2"`, and so on are
#'   used. Default is `NULL`.
#' @param model A list containing character vectors specifying the IRT models
#'   used to calibrate items across all groups. For example, if five groups are
#'   analyzed, the list should contain five elements, each being a character
#'   vector of IRT model names for one group. The elements in the list must be
#'   ordered to match the group names specified in the `group.name` argument.
#'
#'   Available IRT models include:
#'   - `"1PLM"`, `"2PLM"`, `"3PLM"`, `"DRM"` for dichotomous items
#'   - `"GRM"`, `"GPCM"` for polytomous items
#'
#'   Here, `"GRM"` denotes the graded response model and `"GPCM"` the
#'   (generalized) partial credit model. Note that `"DRM"` is a general label
#'   for dichotomous items; items labeled `"DRM"` are fitted with the 3PLM, and
#'   a warning is issued. If a single model name is provided in any element of
#'   the list, it will be recycled across all items within that group.
#'
#'   This argument is used only when `x = NULL` and `fipc = FALSE`. Default
#'   is `NULL`.
#' @param cats A list containing numeric vectors specifying the number of score
#'   categories for items in each group. For example, if five groups are
#'   analyzed, the list should contain five numeric vectors corresponding to
#'   the five groups. The elements in the list must be ordered consistently with
#'   the group names specified in the `group.name` argument.
#'
#'   If a single numeric value is specified in any element of the list, it will
#'   be recycled across all items in the corresponding group. If `cats = NULL`
#'   and all models specified in the `model` argument are dichotomous
#'   (i.e., "1PLM", "2PLM", "3PLM", or "DRM"), the function assumes that all
#'   items have two score categories across all groups.
#'
#'   This argument is used only when `x = NULL` and `fipc = FALSE`. Default
#'   is `NULL`.
#' @param item.id A list containing character vectors of item IDs for each group
#'   to be analyzed. For example, if five groups are analyzed, the list should
#'   contain five character vectors corresponding to the five groups. The
#'   elements in the list must be ordered consistently with the group names
#'   specified in the `group.name` argument.
#'
#'   Items with the same ID in different groups are treated as common items. A
#'   provided `item.id` replaces any IDs in `x`, with or without FIPC, and each
#'   element must have one value per item of its group. With `fipc = TRUE`,
#'   `fix.id` then refers to the new IDs. It is required when `x = NULL`.
#'   Default is `NULL`.
#' @param free.group A numeric or character vector of the positions or names of
#'   the groups whose latent ability distributions are freely estimated: their
#'   means and variances, or their empirical histograms when `EmpHist = TRUE`.
#'   The distributions of the remaining groups (reference groups, those not
#'   specified in this argument) are fixed using the values provided in the
#'   `group.mean` and `group.var` arguments, or from the `weights` argument.
#'   Each group in `free.group` must share common items with another group or,
#'   with FIPC, include fixed items.
#'
#'   For example, suppose that five groups are analyzed with group names "G1",
#'   "G2", "G3", "G4", and "G5". To freely estimate the distributions of groups
#'   2 through 5, set `free.group = c(2, 3, 4, 5)` or
#'   `free.group = c("G2", "G3", "G4", "G5")`. In this case, the first group
#'   ("G1") will have a fixed distribution (e.g., a mean of 0 and variance of 1
#'   when `group.mean = 0`, `group.var = 1`, and `weights = NULL`).
#'
#'   Without FIPC, at least one group must be a reference group to identify the
#'   scale; if `free.group = NULL`, all groups are reference groups. With FIPC,
#'   `free.group` must be supplied, and it usually includes all groups so that
#'   every distribution is estimated on the scale of the fixed items (see
#'   Details). Default is `NULL`.
#' @param weights A two-column matrix or data frame containing the quadrature
#'   points (in the first column) and the corresponding weights (in the second
#'   column) that define the initial latent ability distribution of every group.
#'   The weights and points can be conveniently generated using the function
#'   [irtQ::gen.weight()]. The weights should sum to 1; weights that do not sum
#'   to 1 are rescaled to sum to 1 with a warning. The number of points may
#'   differ from `Quadrature[1]`. If supplied, `Quadrature`, `group.mean`, and
#'   `group.var` are not used. For the groups not specified in `free.group`, the
#'   mean and variance of this distribution fix the scale: the distribution is
#'   used as it is when `EmpHist = FALSE`, and its shape is updated but
#'   rescaled to this mean and variance when `EmpHist = TRUE`, in which case the
#'   quadrature points must be equally spaced, because the histogram is rescaled
#'   by interpolation on an equally spaced grid. For the groups in `free.group`,
#'   it is only the initial distribution.
#'
#'   If `NULL`, a normal prior density is used instead, based on the
#'   information provided in the `Quadrature`, `group.mean`, and `group.var`
#'   arguments. Default is `NULL`.
#' @param group.mean A numeric value specifying the mean of the normal latent
#'   ability distribution when `weights = NULL`. Default is 0. The means of the
#'   groups not specified in the `free.group` argument are fixed at this value
#'   to resolve the indeterminacy of the item parameter scale, also when
#'   `EmpHist = TRUE`. For the groups in `free.group`, it is the mean of the
#'   initial distribution.
#' @param group.var A positive numeric value specifying the variance of the
#'   normal latent ability distribution when `weights = NULL`. Default is 1.
#'   The variances of the groups not specified in the `free.group` argument are
#'   fixed at this value to resolve the indeterminacy of the item parameter
#'   scale, also when `EmpHist = TRUE`. For the groups in `free.group`, it is
#'   the variance of the initial distribution.
#' @param EmpHist Logical. If `TRUE`, the latent ability distribution of each
#'   group is estimated as an empirical histogram on the quadrature points
#'   together with the item parameters, using the approach proposed by Woods
#'   (2007). The histograms of the groups not in `free.group` are rescaled to
#'   the fixed mean and variance in every EM cycle, and those of the groups in
#'   `free.group` are not rescaled. If `FALSE`, the distribution of each group
#'   in `free.group` is a normal distribution whose mean and variance are
#'   updated in every EM cycle. Default is `FALSE`.
#' @param Etol A positive numeric value specifying the convergence criterion
#'   of the EM algorithm. The EM cycles stop when the largest absolute change
#'   in the estimated item parameters between two consecutive cycles is less
#'   than or equal to `Etol`. With MG-FIPC, only the parameters of the items
#'   that are not fixed are compared. The means and variances of the latent
#'   ability distributions of the groups in `free.group`, and the weights of
#'   the empirical histograms (`EmpHist = TRUE`), are not part of the
#'   criterion, so they can be less converged than the item parameters; a
#'   smaller `Etol` gives more precise estimates of them. In particular, an
#'   estimated histogram can keep changing very slowly over many cycles
#'   (toward a spiky distribution) while the log-likelihood barely changes, so
#'   `group.par` and `weights` then depend on `Etol`; with small samples, the
#'   item parameters can then also need many cycles, or more than `MaxE`, to
#'   meet `Etol`. When all items are fixed, the criterion is applied to the
#'   largest absolute change in the means and variances of the latent ability
#'   distributions of all groups; those of the groups not in `free.group` stay
#'   at their fixed values, apart from small changes due to the rescaling of
#'   their histograms when `EmpHist = TRUE`. `Etol` is not used when
#'   `fipc.method = "OEM"` with items to calibrate. The printed output of the
#'   fitted object refers to this criterion as the convergence criterion of the
#'   E-step. Default is `1e-3`, which is larger than the default of
#'   [irtQ::est_irt()].
#' @param fipc Logical. If `TRUE`, multiple-group fixed item parameter
#'   calibration (MG-FIPC) is applied during item parameter estimation.
#'   When `fipc = TRUE`, the information on which items are fixed
#'   must be provided via either `fix.loc` or `fix.id`. See Details.
#' @param fix.loc A list of positive integer vectors. Each internal vector
#'   specifies the positions of the items to be fixed in the item metadata
#'   (i.e., `x`) for each group when MG-FIPC is implemented (i.e.,
#'   `fipc = TRUE`). The internal objects in the list must follow the same order
#'   as the group names provided in the `group.name` argument. Use `NULL` for a
#'   group without fixed items. An item fixed in one group is fixed in every
#'   group that includes it.
#'
#'   For example, suppose three groups are analyzed. In the first group,
#'   the 1st, 3rd, and 5th items are fixed; in the second group, the 2nd, 3rd,
#'   4th, and 7th items are fixed; and in the third group, the 1st, 2nd, and
#'   6th items are fixed. Then
#'   `fix.loc = list(c(1, 3, 5), c(2, 3, 4, 7), c(1, 2, 6))`.
#'   Note that if the `fix.id` argument is not NULL, the information in
#'   `fix.loc` will be ignored. See Details.
#' @param fix.id A vector of character strings specifying the IDs of items to be
#'   fixed when MG-FIPC is implemented (i.e., `fipc = TRUE`). It is a single
#'   character vector for all groups, not a list.
#'
#'   For example, suppose that three groups are analyzed. In the first group,
#'   three items with IDs G1I1, C1I1, and C1I2 are fixed. In the second group,
#'   four items with IDs C1I1, C1I2, C2I1, and C2I2 are fixed. In the third
#'   group, three items with IDs C2I1, C2I2, and G3I1 are fixed.
#'
#'   In this case, there are six unique items fixed across the groups
#'   (G1I1, C1I1, C1I2, C2I1, C2I2, and G3I1), because C1I1 and C1I2 appear in
#'   both the first and second groups, while C2I1 and C2I2 appear in both the
#'   second and third groups. Thus, you should specify
#'   `fix.id = c("G1I1", "C1I1", "C1I2", "C2I1", "C2I2", "G3I1")`.
#'   Note that if the `fix.id` argument is not NULL, the information provided in
#'   `fix.loc` is ignored. See Details.
#'
#' @details Multiple-group (MG) item calibration (Bock & Zimowski, 1997)
#' provides a unified framework for handling testing scenarios involving
#' multiple groups, such as nonequivalent groups equating, vertical scaling,
#' and the identification of differential item functioning (DIF). In such
#' applications, examinees from different groups typically respond to either
#' the same test form or to different forms that share common (anchor) items.
#'
#' The goal of MG item calibration is to estimate both item parameters and
#' latent ability distributions for all groups simultaneously (Bock & Zimowski,
#' 1997). The \pkg{irtQ} package implements MG calibration via the
#' [irtQ::est_mg()] function, which uses marginal maximum likelihood estimation
#' through the expectation-maximization (MMLE-EM) algorithm (Bock & Aitkin,
#' 1981).
#' In addition, the function supports multiple-group fixed item parameter
#' calibration (MG-FIPC; e.g., Kim & Kolen, 2016), which allows the parameters
#' of specific items to be fixed across groups.
#'
#' In MG IRT analyses, it is common for multiple groups' test forms to share
#' some common (anchor) items. The [irtQ::est_mg()] function treats items with
#' the same item ID in different groups as common items and constrains them to
#' have the same parameter estimates.
#'
#' Most of the features of the [irtQ::est_mg()] function are similar to those of
#' the [irtQ::est_irt()] function. The main difference is that the arguments
#' `x`, `data`, `model`, `cats`, `item.id`, and `fix.loc` accept list objects
#' containing elements for each group to be analyzed, in the order of the group
#' names in `group.name`. In contrast, `fix.id` is a single character vector of
#' the IDs of the fixed items across all groups.
#'
#' Additionally, [irtQ::est_mg()] introduces two new arguments: `group.name`
#' and `free.group`. The `group.name` argument assigns a unique identifier to
#' each group.
#'
#' The `free.group` argument indicates the groups whose latent ability
#' distributions are freely estimated; the distributions of the other groups
#' are fixed using the values provided in the `group.mean` and `group.var`
#' arguments or from the `weights` argument. When no item parameters are fixed
#' (i.e., `fipc = FALSE`), at least one group must have a fixed latent ability
#' distribution (e.g., mean = 0 and variance = 1) among the multiple groups
#' sharing common items, in order to resolve the scale indeterminacy inherent
#' in IRT estimation.
#'
#' Situations requiring the implementation of MG-FIPC typically arise when
#' new latent ability scales from multiple-group (MG) test data need to be
#' linked to an established scale (e.g., that of an existing item bank).
#' In a single run of the MG-FIPC procedure, the parameters of non-fixed
#' (freed) items across multiple test forms, as well as the latent ability
#' distributions for multiple groups, can be estimated on the same scale
#' as the fixed items (Kim & Kolen, 2016).
#'
#' For example, suppose that three different test forms (Form 1, Form 2, and
#' Form 3) are administered to three nonequivalent groups: Group1, Group2,
#' and Group3. Form 1 and Form 2 share 12 common items (C1I1 to C1I12),
#' while Form 2 and Form 3 share 10 common items (C2I1 to C2I10). There are
#' no common items between Form 1 and Form 3. Also, assume that all unique
#' items in Form 1 are from an existing item bank and have already been
#' calibrated on the item bank's scale.
#'
#' In this case, the goal of MG-FIPC is to estimate the parameters of all
#' items across the three test forms, except the unique items in Form 1,
#' and the latent ability distributions of the three groups, all on the
#' same scale as the item bank. To achieve this, the unique items in Form 1
#' must be fixed during MG-FIPC to link the current MG test data to the
#' item bank scale.
#'
#' The [irtQ::est_mg()] function can implement MG-FIPC by setting
#' `fipc = TRUE`. In this case, the information on which items to fix must be
#' provided through either the `fix.loc` or `fix.id` argument. When using
#' `fix.loc`, you must supply a list of item positions (locations)
#' to be fixed in each group's test form. For example, suppose that the test
#' data from the three groups above are analyzed. In the first group,
#' the 1st, 3rd, and 5th items are fixed; in the second group, the 2nd, 3rd,
#' 4th, and 7th items are fixed; and in the third group, the 1st, 2nd, and
#' 6th items are fixed. In this case, you would specify:
#' `fix.loc = list(c(1, 3, 5), c(2, 3, 4, 7), c(1, 2, 6))`.
#'
#' Alternatively, you can use `fix.id` to specify a character vector of item IDs
#' to be fixed across groups. In the first group, the items with IDs G1I1, C1I1,
#' and C1I2 are fixed; in the second group, the items with IDs C1I1, C1I2,
#' C2I1, and C2I2 are fixed; and in the third group, the items with IDs C2I1,
#' C2I2, and G3I1 are fixed. In this case, there are six unique items to be
#' fixed across all groups: G1I1, C1I1, C1I2, C2I1, C2I2, and G3I1. You would
#' then specify:
#' `fix.id = c("G1I1", "C1I1", "C1I2", "C2I1", "C2I2", "G3I1")`.
#'
#' Note that when both `fix.loc` and `fix.id` are provided, the information
#' in `fix.id` takes precedence and overrides `fix.loc`.
#'
#' Two methods are available for MG-FIPC through the `fipc.method` argument,
#' as in [irtQ::est_irt()]. The method `"MEM"` runs multiple EM cycles. The
#' posterior distributions of the first EM cycle are computed from the
#' responses to the fixed items only, and those of the later cycles from the
#' responses to all items, as in Kim and Kolen (2016, Equations 13 to 20). The
#' latent ability distribution of each group in `free.group` is updated in
#' every EM cycle without rescaling, so that the item parameters and the
#' distributions are placed on the scale of the fixed items. With
#' `EmpHist = TRUE`, the distribution is the histogram of the expected numbers
#' of examinees at the quadrature points, as in Kim and Kolen (2016); with
#' `EmpHist = FALSE`, it is a normal distribution with the mean and variance of
#' that histogram. The distributions of the groups not in `free.group` remain
#' fixed also in MG-FIPC, whereas Kim and Kolen (2016) estimate the
#' distributions of all groups; to follow their method, include all groups in
#' `free.group`. The method `"OEM"` (Wainer & Mislevy, 1990) runs a single EM
#' cycle: the posterior distributions are computed from the fixed items and the
#' initial distributions, the parameters of the other items are estimated once,
#' and the returned distributions are updated once afterward.
#'
#' A maximum likelihood estimate does not exist for an item with a score
#' category that has no responses, which includes an item whose responses are
#' all the same. The estimates of such an item can be extreme or remain near
#' their starting values, and their standard errors are unreliable (a standard
#' error of 99999 means that it could not be obtained). Before the estimation,
#' one warning lists the IDs of these items. Because the parameters of an item
#' are common to all the groups that include it, the responses of an item are
#' pooled over these groups; an item that is constant in one group but varies in
#' another is not flagged. The items fixed in MG-FIPC are not checked. A
#' parameter with a prior distribution, such as the guessing parameter of a
#' 3PLM item with `use.gprior = TRUE`, can still be stable.
#'
#' When prior distributions are used, the item parameters are Bayes modal
#' estimates, and `loglikelihood` is the log-likelihood at these estimates
#' without the log prior densities.
#'
#' The information criteria are computed as
#' `aic = -2 * loglikelihood$overall + 2 * npar.est` and
#' `bic = -2 * loglikelihood$overall + npar.est * log(ncase$overall)`,
#' where `loglikelihood$overall` is the overall marginal log-likelihood
#' without the prior densities of the item parameters and `ncase$overall`
#' is the overall number of examinees. `npar.est` is the number
#' of freely estimated parameters, which consists of (a) the item parameters
#' that are estimated across all unique items (a slope shared by 1PLM items is
#' counted once; fixed parameters and the items fixed in FIPC are not counted),
#' (b) the mean and variance of each group in `free.group` when
#' `EmpHist = FALSE` (2 values per group), and (c) the weights of the empirical
#' histograms when `EmpHist = TRUE`. The histogram of a group in `free.group`
#' is not rescaled, so its Q weights, where Q is the number of quadrature
#' points, have Q - 1 free values (the sum of the weights is fixed). The
#' histogram of any other group is rescaled to a fixed mean and variance, so
#' its weights have Q - 3 free values. The counts by type are returned in
#' `npar.detail`.
#'
#'   The print method of this class has an argument `digits` that is not used.
#'   It prints the log-likelihood with the number of significant digits given
#'   by `getOption("digits")` (7 by default).
#'
#' @return This function returns an object of class `est_mg`. The returned
#'   object contains the following components:
#'
#' \item{estimates}{A list containing two internal elements: `overall` and
#'   `group`. The `overall` element is a data frame with item parameter
#'   estimates and their standard errors, computed from the combined data across
#'   all groups. This data frame includes only the unique items across all
#'   groups. The `group` element is a list of group-specific data frames, each
#'   containing item parameter estimates and standard errors for that particular
#'   group.}
#'
#' \item{par.est}{A list with the same structure as `estimates`, containing only
#'   the item parameter estimates (excluding standard errors), formatted
#'   according to the item metadata structure.}
#'
#' \item{se.est}{A list with the same structure as `estimates`, but containing
#'   only the standard errors of the item parameter estimates. Note that the
#'   standard errors are computed from the inverse of the information matrix,
#'   which is the cross-product (empirical) information matrix of the marginal
#'   log-likelihood evaluated at the final estimates (Meilijson, 1989), plus
#'   the prior information (the negative second derivatives of the log prior
#'   densities) for the parameters that have prior distributions. The latent
#'   ability distributions of the groups are treated as known in this matrix,
#'   so the uncertainty of the estimated distributions is not reflected in
#'   these standard errors. Parameters that are not estimated, including those
#'   of the items fixed in FIPC, have `NA` standard errors, and all are `NA`
#'   when `se = FALSE`.}
#'
#' \item{pos.par}{A data frame of the positions of the estimated item
#'   parameters in `covariance`, for the unique items across all groups (the
#'   `overall` element of `estimates`). Parameters that are not estimated,
#'   including those of the items fixed in FIPC, are `NA`. It is `NULL` when
#'   all items are fixed in FIPC.}
#'
#' \item{covariance}{The variance-covariance matrix of the estimated item
#'   parameters of the unique items across all groups, in the order given by
#'   `pos.par`. It is `NULL` when `se = FALSE` or when all items are fixed in
#'   FIPC.}
#'
#' \item{loglikelihood}{A list with two elements. `overall` is the marginal
#'   log-likelihood of the response data at the final estimates: the sum over
#'   all examinees of the logarithm of the probability of each examinee's
#'   response pattern, integrated over the latent ability distribution of the
#'   examinee's group. `group` is a named list of the same sums within the
#'   groups, which add up to `overall`. The log prior densities of the item
#'   parameters are not included.}
#'
#' \item{aic}{Akaike information criterion,
#'   `-2 * loglikelihood$overall + 2 * npar.est`,
#'   calculated from the overall log-likelihood.}
#'
#' \item{bic}{Bayesian information criterion,
#'   `-2 * loglikelihood$overall + npar.est * log(ncase$overall)`, calculated
#'   from the overall log-likelihood and the overall number of examinees.}
#'
#' \item{group.par}{A named list with one data frame per group, with rows
#'   `estimates` and `se` and columns `mu`, `sigma2`, and `sigma` (the mean,
#'   variance, and standard deviation of the latent ability distribution). For
#'   a group not in `free.group`, the estimates are the fixed values given by
#'   `group.mean` and `group.var` (or by `weights`), and the standard errors
#'   are `NA`; with `EmpHist = TRUE`, the histogram of such a group is
#'   rescaled in each cycle to have this mean and variance, but because the
#'   rescaling interpolates on the quadrature points, the mean and variance
#'   computed from the returned `weights` can differ slightly from the
#'   reported values. For a group in `free.group`, the estimates are the mean
#'   and variance of the final distribution (the same as those computed from
#'   the returned `weights`), and the standard errors are approximations that
#'   treat the abilities as observed and are therefore usually too small:
#'   `sigma / sqrt(N)`, `sigma2 * sqrt(2 / (N - 1))`, and
#'   `se(sigma2) / (2 * sigma)`, where `N` is the number of examinees in the
#'   group.}
#'
#' \item{weights}{A list of two-column data frames, where the first column
#'   contains quadrature points and the second column contains the
#'   corresponding weights of the (updated) latent variable prior distributions
#'   for each group.}
#'
#' \item{posterior.dist}{A matrix of the posterior probabilities of the
#'   quadrature points for each examinee, computed with the final estimates and
#'   the latent distribution of the examinee's group. Rows represent the
#'   examinees of all groups, stacked in the order of `group.name`, and columns
#'   represent quadrature points; each row sums to 1.}
#'
#' \item{data}{A list with two elements. `overall` is a matrix of the
#'   responses of all examinees to all unique items, with `NA` for the items
#'   not in a group's form and for missing responses. `group` is the list of
#'   response data sets named by `group.name`, in which the missing responses
#'   are also recoded as `NA`.}
#'
#' \item{scale.D}{The scaling factor used in the IRT model.}
#'
#' \item{ncase}{A list with two elements: `overall`, the total number of
#'   examinees, and `group`, a named vector of the numbers of examinees in the
#'   groups.}
#'
#' \item{nitem}{A list with two elements: `overall`, the number of unique items
#'   across all groups, and `group`, a named vector of the numbers of items in
#'   the groups.}
#'
#' \item{Etol}{The convergence criterion of the EM algorithm; see `Etol` in
#'   Arguments.}
#'
#' \item{MaxE}{The maximum number of EM cycles (1 when
#'   `fipc.method = "OEM"`).}
#'
#' \item{aprior}{A list describing the prior distribution used for
#'   discrimination parameters. It is `NULL` when `use.aprior = FALSE`.}
#'
#' \item{bprior}{A list describing the prior distribution used for difficulty
#'   parameters. It is `NULL` when `use.bprior = FALSE`.}
#'
#' \item{gprior}{A list describing the prior distribution used for guessing
#'   parameters. It is `NULL` when `use.gprior = FALSE`.}
#'
#' \item{npar.est}{The number of freely estimated parameters used for `aic` and
#'   `bic`: the estimated item parameters of all unique items, the means and
#'   variances of the free groups when they have normal distributions, and the
#'   weights of the empirical histograms when `EmpHist = TRUE`. See Details.}
#'
#' \item{niter}{The number of completed EM cycles.}
#'
#' \item{maxpar.diff}{The largest absolute change in the estimated item
#'   parameters in the last EM cycle (for MG-FIPC with all items fixed, in the
#'   means and variances of the latent ability distributions of all groups; see
#'   `Etol`). With `fipc.method = "OEM"`, it is the change from the starting
#'   values.}
#'
#' \item{EMtime}{Time (in seconds) spent on EM cycles.}
#'
#' \item{SEtime}{Time (in seconds) spent computing standard errors; `NULL`
#'   when they are not computed.}
#'
#' \item{TotalTime}{Total computation time (in seconds).}
#'
#' \item{test.1}{A message indicating whether the convergence criteria were
#'   met: the M-step optimization converged for every item and the largest
#'   absolute change in the parameter estimates between two consecutive EM
#'   cycles was less than or equal to `Etol`. When `fipc.method = "OEM"` with
#'   new items to calibrate, the single EM cycle is not judged by `Etol` and
#'   only the M-step check applies. For FIPC with all items fixed, see `Etol`.}
#'
#' \item{test.2}{A message indicating whether the information matrix used
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
#' \item{var.note}{A note indicating whether the variance-covariance matrix
#'   was successfully obtained from the information matrix.}
#'
#' \item{fipc}{Logical. Indicates whether FIPC was used.}
#'
#' \item{fipc.method}{The method used for FIPC, or `NULL` when `fipc = FALSE`.}
#'
#' \item{fix.loc}{A list with two elements: `overall`, the positions of the
#'   fixed items among the unique items across all groups (the rows of
#'   `par.est$overall`), and `group`, a named list of the positions of the
#'   fixed items in each group's form. Both are `NULL` when `fipc = FALSE`.}
#'
#' \item{call}{The matched call.}
#'
#' \item{npar.detail}{An integer vector with the elements `item` (free item
#'   parameters), `dist` (means and variances of the latent distributions of
#'   all groups), and `weight` (weights of the empirical histograms of all
#'   groups), whose sum is `npar.est`. It is the last element of the object;
#'   when it is absent, `summary()` shows the total only.}
#'
#'   Note that you can easily extract components from the output using the
#'   [irtQ::getirt()] function.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [irtQ::est_irt()], [irtQ::shape_df()], [irtQ::shape_df_fipc()],
#' [irtQ::getirt()]
#'
#' @references Bock, R. D., & Aitkin, M. (1981). Marginal maximum likelihood
#' estimation of item parameters: Application of an EM algorithm.
#' *Psychometrika, 46*, 443-459.
#'
#'   Bock, R. D., & Zimowski, M. F. (1997). Multiple group IRT. In W. J. van der
#'   Linden & R. K. Hambleton (Eds.), *Handbook of modern item response theory*
#'   (pp. 433-448). Springer. \doi{10.1007/978-1-4757-2691-6_25}.
#'
#' Kim, S. (2006). A comparative study of IRT fixed parameter calibration
#' methods. *Journal of Educational Measurement, 43*(4), 355-381.
#'
#' Kim, S., & Kolen, M. J. (2016). *Multiple group IRT fixed-parameter
#' estimation for maintaining an established ability scale* (CASMA Research
#' Report No. 49). Center for Advanced Studies in Measurement and Assessment,
#' University of Iowa.
#'
#' Meilijson, I. (1989). A fast improvement to the EM algorithm on its own
#' terms. *Journal of the Royal Statistical Society: Series B
#' (Methodological), 51*, 127-138.
#'
#' Wainer, H., & Mislevy, R. J. (1990). Item response theory, item
#' calibration, and proficiency estimation. In H. Wainer (Ed.), *Computer
#' adaptive testing: A primer* (Chap. 4, pp. 65-102). Hillsdale, NJ: Lawrence
#' Erlbaum.
#'
#' Woods, C. M. (2007). Empirical histograms in item response theory with
#' ordinal data. *Educational and Psychological Measurement, 67*(1), 73-87.
#'
#' @examples
#' \donttest{
#' ## ------------------------------------------------------------------------------
#' # 1. MG calibration using the simMG data
#' #  - Details:
#' #    (a) Constrain common items between groups to have
#' #        identical item parameters (i.e., items C1I1 - C1I12 between
#' #        Groups 1 and 2, and items C2I1 - C2I10 between Groups 2 and 3).
#' #    (b) Freely estimate the means and variances of the ability
#' #        distributions for all groups except the reference group,
#' #        where the mean and variance are fixed to 0 and 1, respectively.
#' ## ------------------------------------------------------------------------------
#' # 1-(1). Freely estimate the means and variances of Groups 2 and 3
#' # Import the true item metadata for the three groups
#' x <- simMG$item.prm
#'
#' # Extract model, score category, and item ID information
#' # from the item metadata for the three groups
#' model <- list(x$Group1$model, x$Group2$model, x$Group3$model)
#' cats <- list(x$Group1$cats, x$Group2$cats, x$Group3$cats)
#' item.id <- list(x$Group1$id, x$Group2$id, x$Group3$id)
#'
#' # Import the simulated response data sets for the three groups
#' data <- simMG$res.dat
#'
#' # Import the group names for the three groups
#' group.name <- simMG$group.name
#'
#' # Specify Groups 2 and 3 as the free groups where the scale
#' # of the ability distributions will be freely estimated.
#' # Group 1 will serve as the reference group, where the scale
#' # of the ability distribution is fixed to the values specified
#' # via the 'group.mean' and 'group.var' arguments
#' free.group <- c(2, 3) # or use 'free.group <- group.name[2:3]'
#'
#' # Estimate IRT parameters:
#' # As long as common items across groups share the same item IDs,
#' # their item parameters will be constrained to be equal across groups
#' # unless FIPC is implemented.
#' fit.1 <-
#'   est_mg(
#'     data = data, group.name = group.name, model = model,
#'     cats = cats, item.id = item.id, D = 1, free.group = free.group,
#'     use.gprior = TRUE, gprior = list(dist = "beta", params = c(5, 16)),
#'     group.mean = 0, group.var = 1, EmpHist = TRUE, Etol = 0.001, MaxE = 500
#'   )
#'
#' # Summary of the estimation
#' summary(fit.1)
#'
#' # Extract the item parameter estimates
#' getirt(fit.1, what = "par.est")
#'
#' # Extract the standard error estimates
#' getirt(fit.1, what = "se.est")
#'
#' # Extract the group-level parameter estimates (i.e., scale parameters)
#' getirt(fit.1, what = "group.par")
#'
#' # Extract the posterior latent ability distributions for each group
#' getirt(fit.1, what = "weights")
#'
#' # 1-(2). Alternatively, the same parameter estimation can be performed by
#' # inserting a list of item metadata for the groups into the 'x' argument.
#' # If the item metadata contains item parameters to be used as starting values,
#' # set 'use.startval = TRUE'.
#' # Also, specify the groups in which the ability distribution scales
#' # will be freely estimated using their group names.
#' free.group <- group.name[2:3]
#' fit.2 <-
#'   est_mg(
#'     x = x, data = data, group.name = group.name, D = 1,
#'     free.group = free.group, use.gprior = TRUE,
#'     gprior = list(dist = "beta", params = c(5, 16)),
#'     group.mean = 0, group.var = 1, EmpHist = TRUE, use.startval = TRUE,
#'     Etol = 0.001, MaxE = 500
#'   )
#'
#' # Summary of the estimation
#' summary(fit.2)
#'
#' ## ------------------------------------------------------------------------------
#' # 2. MG calibration with FIPC using simMG data
#' #  - Details:
#' #    (a) Fix the parameters of the common items between the groups
#' #        (i.e., items C1I1 - C1I12 between Groups 1 and 2, and
#' #        items C2I1 - C2I10 between Groups 2 and 3)
#' #    (b) Freely estimate the means and variances of the ability
#' #        distributions for all three groups
#' ## ------------------------------------------------------------------------------
#' # 2-(1). Freely estimate the means and variances for all three groups
#' # Set all three groups as free groups in which the scales
#' # of the ability distributions are to be freely estimated
#' free.group <- 1:3 # or use 'free.group <- group.name'
#'
#' # Specify the locations of items to be fixed in each group's metadata
#' # For Group 1: C1I1 - C1I12 are located in rows 1-10 and 49-50
#' # For Group 2: C1I1 - C1I12 are in rows 1-12, and
#' #              C2I1 - C2I10 are in rows 41-50
#' # For Group 3: C2I1 - C2I10 are in rows 1-10
#' fix.loc <- list(
#'   c(1:10, 49:50),
#'   c(1:12, 41:50),
#'   c(1:10)
#' )
#'
#' # Estimate IRT parameters using MG-FIPC:
#' # When FIPC is implemented, item metadata for all groups
#' # must be provided via the 'x' argument.
#' # For fixed items, their item parameters must be specified
#' # in the metadata. For non-fixed items, any placeholder values
#' # can be used in the metadata.
#' # Also set fipc = TRUE and fipc.method = "MEM"
#' fit.3 <-
#'   est_mg(
#'     x = x, data = data, group.name = group.name, D = 1,
#'     free.group = free.group, use.gprior = TRUE,
#'     gprior = list(dist = "beta", params = c(5, 16)),
#'     EmpHist = TRUE, Etol = 0.001, MaxE = 500, fipc = TRUE,
#'     fipc.method = "MEM", fix.loc = fix.loc
#'   )
#'
#' # Summary of the estimation
#' summary(fit.3)
#'
#' # Extract the item parameter estimates
#' getirt(fit.3, what = "par.est")
#'
#' # Extract the standard error estimates
#' getirt(fit.3, what = "se.est")
#'
#' # Extract the group parameter estimates (i.e., scale parameters)
#' getirt(fit.3, what = "group.par")
#'
#' # Extract the posterior latent ability distributions of the groups
#' getirt(fit.3, what = "weights")
#'
#' # 2-(2). Alternatively, MG-FIPC can be implemented by specifying the
#' # IDs of the items to be fixed using the 'fix.id' argument.
#' # Provide a character vector of fixed item IDs to 'fix.id'
#' fix.id <- c(paste0("C1I", 1:12), paste0("C2I", 1:10))
#' fit.4 <-
#'   est_mg(
#'     x = x, data = data, group.name = group.name, D = 1,
#'     free.group = free.group, use.gprior = TRUE,
#'     gprior = list(dist = "beta", params = c(5, 16)),
#'     EmpHist = TRUE, Etol = 0.001, MaxE = 500, fipc = TRUE,
#'     fipc.method = "MEM", fix.id = fix.id
#'   )
#'
#' # Summary of the estimation
#' summary(fit.4)
#'
#' ## ------------------------------------------------------------------------------
#' # 3. MG calibration with FIPC using simMG data
#' #    (Estimate group parameters only)
#' #  - Details:
#' #    (a) Fix all item parameters across all three groups
#' #    (b) Freely estimate the means and variances of the ability
#' #        distributions for all three groups
#' ## ------------------------------------------------------------------------------
#' # 3-(1). Freely estimate the means and variances for all three groups
#' # Set all three groups as free groups in which the scales
#' # of the ability distributions will be freely estimated
#' free.group <- 1:3 # or use 'free.group <- group.name'
#'
#' # Specify the locations of all fixed items in each group's metadata
#' fix.loc <- list(1:50, 1:50, 1:38)
#'
#' # Estimate group parameters only using MG-FIPC
#' fit.5 <-
#'   est_mg(
#'     x = x, data = data, group.name = group.name, D = 1,
#'     free.group = free.group, use.gprior = TRUE,
#'     gprior = list(dist = "beta", params = c(5, 16)),
#'     EmpHist = TRUE, Etol = 0.001, MaxE = 500, fipc = TRUE,
#'     fipc.method = "MEM", fix.loc = fix.loc
#'   )
#'
#' # Summary of the estimation
#' summary(fit.5)
#'
#' # Extract the group parameter estimates (i.e., scale parameters)
#' getirt(fit.5, what = "group.par")
#'
#' ## ------------------------------------------------------------------------------
#' # 4. MG calibration with FIPC using simMG data
#' #    (Fix only the unique items of Group 1)
#' #  - Details:
#' #    (a) Fix item parameters of the unique items in Group 1 only
#' #    (b) Constrain the common items across groups to have
#' #        the same item parameters (i.e., C1I1 - C1I12 between
#' #        Groups 1 and 2, and C2I1 - C2I10 between Groups 2 and 3)
#' #    (c) Freely estimate the means and variances of the ability
#' #        distributions for all three groups
#' ## ------------------------------------------------------------------------------
#' # 4-(1). Freely estimate the means and variances for all three groups
#' # Set all three groups as free groups in which the scales
#' # of the ability distributions will be freely estimated
#' free.group <- group.name # or use 'free.group <- 1:3'
#'
#' # Specify the item IDs of the unique items in Group 1 to be fixed using
#' # the `fix.id` argument.
#' fix.id <- paste0("G1I", 1:38)
#'
#' # Alternatively, use the 'fix.loc' argument as
#' # 'fix.loc = list(11:48, NULL, NULL)'
#'
#' # Estimate IRT parameters using MG-FIPC
#' fit.6 <-
#'   est_mg(
#'     x = x, data = data, group.name = group.name, D = 1,
#'     free.group = free.group, use.gprior = TRUE,
#'     gprior = list(dist = "beta", params = c(5, 16)),
#'     EmpHist = TRUE, Etol = 0.001, MaxE = 500, fipc = TRUE,
#'     fipc.method = "MEM", fix.loc = NULL, fix.id = fix.id
#'   )
#'
#' # Summary of the estimation
#' summary(fit.6)
#'
#' # Extract the group parameter estimates (i.e., scale parameters)
#' getirt(fit.6, what = "group.par")
#'
#' }
#'
#'
#' @importFrom utils modifyList
#'
#' @export
est_mg <- function(x = NULL,
                   data,
                   group.name = NULL,
                   D = 1,
                   model = NULL,
                   cats = NULL,
                   item.id = NULL,
                   free.group = NULL,
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
                   Etol = 1e-03,
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

  # stop when the arguments given per group do not match the groups in 'data'
  check_mg_args(
    x = x, data = data, group.name = group.name, model = model, cats = cats,
    item.id = item.id, fipc = fipc, fipc.method = fipc.method
  )

  # item parameter estimation
  if (!fipc) {
    # item parameter estimation using MMLE-EM algorithm
    est_par <- est_mg_em(
      x = x, data = data, group.name = group.name, D = D, model = model, cats = cats, item.id = item.id,
      free.group = free.group, fix.a.1pl = fix.a.1pl, fix.a.gpcm = fix.a.gpcm, fix.g = fix.g, a.val.1pl = a.val.1pl,
      a.val.gpcm = a.val.gpcm, g.val = g.val, use.aprior = use.aprior, use.bprior = use.bprior, use.gprior = use.gprior,
      aprior = aprior, bprior = bprior, gprior = gprior, missing = missing, Quadrature = Quadrature, weights = weights,
      group.mean = group.mean, group.var = group.var, EmpHist = EmpHist, use.startval = use.startval, Etol = Etol, MaxE = MaxE,
      control = control, se = se, verbose = verbose
    )
  } else {
    # implement FIPC method
    est_par <- est_mg_fipc(
      x = x, data = data, group.name = group.name, D = D, model = model, cats = cats, item.id = item.id,
      free.group = free.group, fix.a.1pl = fix.a.1pl, fix.a.gpcm = fix.a.gpcm, fix.g = fix.g, a.val.1pl = a.val.1pl,
      a.val.gpcm = a.val.gpcm, g.val = g.val, use.aprior = use.aprior, use.bprior = use.bprior, use.gprior = use.gprior,
      aprior = aprior, bprior = bprior, gprior = gprior, missing = missing, Quadrature = Quadrature, weights = weights,
      group.mean = group.mean, group.var = group.var, EmpHist = EmpHist, use.startval = use.startval, Etol = Etol, MaxE = MaxE,
      control = control, fipc = fipc, fipc.method = fipc.method, fix.loc = fix.loc, fix.id = fix.id, se = se, verbose = verbose
    )
  }

  # return the estimation results
  class(est_par) <- "est_mg"
  est_par$call <- cl

  # place the counts of the free parameters after the call
  npar.detail <- est_par$npar.detail
  est_par$npar.detail <- NULL
  est_par$npar.detail <- npar.detail
  est_par
}


# multiple group calibration via EM
#' @import dplyr
est_mg_em <- function(x = NULL,
                      data,
                      group.name = NULL,
                      D = 1,
                      model = NULL,
                      cats = NULL,
                      item.id = NULL,
                      free.group = NULL,
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
                      Etol = 1e-03,
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

  # start parsing inputs
  if (verbose) {
    cat("Parsing input...", "\n")
  }

  # count the number of groups
  ngroup <- length(data)
  if (ngroup == 1) {
    stop("Use the est_irt() function when a single-group item calibration needs to be implemented.", call. = FALSE)
  }

  # create a group name vector when group.name = NULL
  if (is.null(group.name)) {
    group.name <- paste0("g", 1:ngroup)
  }

  # test length of each group
  nitem.gr <- purrr::map_dbl(.x = data, ncol)
  names(nitem.gr) <- group.name

  # sample size of each group
  nstd.gr <- purrr::map_dbl(.x = data, nrow)
  names(nstd.gr) <- group.name

  # extract information about the number of score categories and models
  if (!is.null(x)) {
    # confirm and correct all item metadata information
    x.gr <- warn_drm_once(purrr::map(.x = x, ~ {
      confirm_df(.x)
    }))
    names(x.gr) <- group.name

    # override the item IDs in x when item.id is supplied
    if (!is.null(item.id)) {
      x.gr <- purrr::map2(.x = x.gr, .y = item.id, ~ {
        .x$id <- .y
        .x
      })
    }

    # item, cats, model information for each group
    item.id <- purrr::map(.x = x.gr, ~ {
      .x$id
    })
    cats.gr <- purrr::map(.x = x.gr, ~ {
      .x$cats
    })
    model.gr <- purrr::map(.x = x.gr, ~ {
      .x$model
    })
    if (!use.startval) {
      x.gr <-
        purrr::map(
          .x = x.gr,
          ~ {
            startval_df(cats = .x$cats, model = .x$model, item.id = .x$id)
          }
        )
    }
  } else {
    # check if the item IDs are specified when x = NULL
    if (is.null(item.id)) {
      stop("A list containing the vectors of the item IDs across all groups must be specified in the argument 'item.id'.",
        call. = FALSE
      )
    }

    # assign group name to the item id list
    names(item.id) <- group.name

    # IRT models for each group test form
    model.gr <-
      purrr::map2(
        .x = model, .y = nitem.gr,
        .f = function(x, y) {
          if (length(x) == 1) {
            rep(x, y)
          } else {
            x
          }
        }
      ) %>%
      purrr::map(~ {
        toupper(.x)
      })
    names(model.gr) <- group.name

    # score categories for each group test form (cats = NULL is handled below)
    if (!is.null(cats)) {
      cats.gr <-
        purrr::map2(
          .x = cats, .y = nitem.gr,
          .f = function(x, y) {
            if (length(x) == 1) {
              rep(x, y)
            } else {
              x
            }
          }
        )
    }

    # create score category vectors when cats = NULL & all models are dichotomous item models
    if (is.null(cats)) {
      cats.gr <-
        purrr::map2(
          .x = nitem.gr, .y = model.gr,
          .f = function(x, y) {
            if (all(y %in% c("1PLM", "2PLM", "3PLM", "DRM"))) {
              rep(2, x)
            } else {
              stop("The number of score categories for the items should be specified in the argument 'cats'.", call. = FALSE)
            }
          }
        )
    }

    # create score category vectors when any of internal objects is Null in the cats argument &
    # all models are dichotomous item models
    cats.gr <-
      purrr::pmap(
        .l = list(x = cats.gr, y = nitem.gr, z = model.gr),
        .f = function(x, y, z) {
          if (is.null(x)) {
            if (all(z %in% c("1PLM", "2PLM", "3PLM", "DRM"))) {
              rep(2, y)
            } else {
              stop("The number of score categories for the items should be specified in the argument 'cats'.", call. = FALSE)
            }
          } else {
            x
          }
        }
      )
    names(cats.gr) <- group.name

    # create the item metadata for each group test form containing starting values
    x.gr <- warn_drm_once(
      purrr::pmap(
        .l = list(x = cats.gr, y = model.gr, z = item.id),
        ~ {
          # check the metadata as est_irt() does (e.g., "DRM" items become "3PLM" items)
          confirm_df(startval_df(cats = ..1, model = ..2, item.id = ..3))
        }
      )
    )
    names(x.gr) <- group.name
  }

  # stop when an item ID is repeated within a group or is defined differently across groups
  check_mg_items(x.gr)

  # copy data to data.gr object
  data.gr <- data
  names(data.gr) <- group.name

  # create combined response data sets across all groups
  # and assign item IDs
  data <-
    # read the responses of each group by their values before the binding
    purrr::map(.x = data.gr, ~ {
      # recode a non-numeric missing code before the responses are read as numbers
      if (!is.na(missing)) .x[.x == missing] <- NA
      data.frame(resp_as_numeric(.x))
    }) %>%
    purrr::map2(
      .y = item.id,
      .f = function(x, y) {
        dplyr::rename_all(.tbl = x, .funs = function(x) y)
      }
    ) %>%
    dplyr::bind_rows()

  # a vector of item ID for the combined data set
  id <- colnames(data)

  # create an item metadata for the combined data set
  x <-
    dplyr::bind_rows(x.gr) %>%
    dplyr::distinct(.data$id, .keep_all = TRUE)
  if (!all(id == x$id)) {
    x <- x[match(id, x$id), ]
  }
  rownames(x) <- 1:nrow(x)

  # create two vectors of model & score categories for the combined data set
  model <- x$model
  cats <- x$cats

  # recode missing values
  if (!is.na(missing)) {
    data[data == missing] <- NA
    # recode the missing values in the response data of each group as well
    data.gr <- purrr::map(.x = data.gr, ~ {
      .x[.x == missing] <- NA
      .x
    })
  }

  # transform a data set to matrix (character and factor columns are read by their values)
  data <- resp_as_numeric(data)

  # check the total number of examinees for the combined data
  nstd <- nrow(data)

  # check the total number of items in the combined data
  nitem <- ncol(data)

  # stop when a response is not one of the score categories 0, ..., cats - 1
  check_resp_range(data = data, cats = cats, id = id)

  # designate a reference group whose group parameters will be fixed
  # (stop when a free group does not exist or no reference group is left)
  grp <- set_mg_groups(free.group = free.group, group.name = group.name, fipc = FALSE)
  free.group <- grp$free.group
  ref.group <- grp$ref.group

  # find the groups that have no common items with other groups
  pair <-
    t(utils::combn(group.name, m = 2)) %>%
    data.frame()
  overlap <- c()
  for (i in 1:nrow(pair)) {
    freq.tmp <-
      unlist(item.id[unlist(pair[i, ])]) %>%
      table()
    overlap[i] <- any(freq.tmp > 1)
  }
  link.gr <-
    dplyr::mutate(pair, overlap = overlap) %>%
    dplyr::filter(.data$overlap == TRUE) %>%
    dplyr::select(1:2) %>%
    unlist() %>%
    unique()
  nolink.gr <- c(1:ngroup)[!group.name %in% link.gr]

  # check if there exist free groups which do not share common items with other groups
  nolink.free <- free.group[free.group %in% nolink.gr]

  # when there exist free groups that do not share any common items,
  # stop the further process
  if (length(nolink.free) > 0) {
    warning.memo <-
      paste0(
        paste(group.name[nolink.free], collapse = ", "),
        " group(s) whose ability distribution(s) is(are) freely estimated do(does) not share common items with other groups. \n",
        "Please specify the freely estimated group(s) correctly."
      )
    stop(paste0(warning.memo, " \n"), call. = FALSE)
  }

  # check the number of item responses across all items
  n.resp <- Rfast::colsums(!is.na(data))

  # check the items which have all missing responses
  loc_allmiss <- which(n.resp == 0L)
  if (length(loc_allmiss) > 0L) {
    allmiss_item <- id[loc_allmiss]
    memo2 <- paste0(paste0(allmiss_item, collapse = ", "), " has/have no item response data across all groups. \n")
    stop(memo2, call. = FALSE)
  }

  # warn about the items whose maximum likelihood estimates do not exist, using the
  # responses pooled over all groups because the parameters are common to the groups
  check_item_mle(data = data, cats = cats, id = id)

  # find the location of 1PLM items in which slope parameters should be constrained to be equal
  # also, find the location of other items
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
    # create quadrature points
    quadpt <- seq(-Quadrature[2], Quadrature[2], length.out = Quadrature[1])

    # create the data.frame containing the quadrature points and weights
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

  # a list containing the weights and densities for each group
  weights.gr <- replicate(n = ngroup, expr = weights, simplify = FALSE)
  names(weights.gr) <- group.name

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

  # create indices of samples (examinees) who belong to each group
  idx.tmp2 <- as.numeric(cumsum(nstd.gr))
  idx.tmp1 <- dplyr::lag(idx.tmp2, default = 0) + 1
  idx.std <- purrr::map2(.x = idx.tmp1, .y = idx.tmp2, ~ {
    .x:.y
  })
  names(idx.std) <- group.name

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
      weights = weights.gr, D = D, idx.std = idx.std
    )

    # implement M-step
    mstep <- Mstep(
      estep = estep, id = id, cats = cats, model = model, quadpt = quadpt, n.quad = n.quad,
      D = D, cols.item = cols.item, loc_1p_const = loc_1p_const, loc_else = loc_else, idx4est = idx4est,
      n.1PLM = n.1PLM, EmpHist = EmpHist, weights = weights.gr, fix.a.1pl = fix.a.1pl, fix.a.gpcm = fix.a.gpcm, fix.g = fix.g,
      a.val.1pl = a.val.1pl, a.val.gpcm = a.val.gpcm, g.val = g.val, use.aprior = use.aprior, use.bprior = use.bprior,
      use.gprior = use.gprior, aprior = aprior, bprior = bprior, gprior = gprior, group.mean = group.mean,
      group.var = group.var, nstd = nstd.gr, control = control,
      iter = r, fipc = FALSE, reloc.par = param_loc$reloc.par, ref.group = ref.group,
      free.group = free.group, parbd = parbd
    )

    # compute the difference between previous and updated item parameter estimates
    diff_par <- mstep$elm_item$pars - elm_item$pars
    # convergence statistic: largest absolute parameter change
    max.diff <- suppressWarnings(max(abs(diff_par), na.rm = TRUE))

    # loglikelihood value
    llike <- do.call(what = "sum", args = mstep$loglike)

    # print
    if (verbose) {
      cat("\r", paste0(
        "EM iteration: ", r, ", Loglike: ",
        format(round(llike, 4), nsmall = 4), ", Max-Change: ",
        format(round(max.diff, 6), nsmall = 5)
      ))
    }

    # check the convergence of EM algorithm; a non-finite change is not converged
    converge <- is.finite(max.diff) && max.diff <= Etol

    # extract the updated item parameter estimates
    elm_item$pars <- mstep$elm_item$pars

    # the quadrature points and the corresponding weights of the prior population density
    weights.gr <- mstep$weights

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
    weights = weights.gr, D = D, idx.std = idx.std
  )

  # compute the final log of marginal likelihood
  likehd.gr <- purrr::map(.x = idx.std, ~ {
    estep$likehd[.x, ]
  })

  # log-scale offsets of the rows rescaled in likelihood(), by group
  offset.gr <- purrr::map(.x = idx.std, ~ {
    estep$offset[.x]
  })
  llike.gr <- purrr::pmap(
    .l = list(likehd.gr, weights.gr, offset.gr),
    .f = function(lk, wt, os) {
      sum(log(lk %*% matrix(wt[, 2]))) + sum(os)
    }
  )
  llike <- do.call(what = "sum", args = llike.gr)

  # compute the mean and variance of the estimated density distributions
  pop_moments <-
    purrr::map(.x = weights.gr, ~ {
      cal_moment(node = quadpt, weight = .x[, 2])
    }) %>%
    purrr::map(~ {
      c(.x, sigma = sqrt(as.numeric(.x[2])))
    })
  if (length(ref.group) > 0) {
    pop_moments[ref.group] <-
      purrr::map(
        .x = 1:length(ref.group),
        ~ {
          c(mu = group.mean, sigma2 = group.var, sigma = sqrt(group.var))
        }
      )
  }

  ## ---------------------------------------------------------------------
  # estimates the information matrix and standard errors of item parameter estimates
  ## ---------------------------------------------------------------------
  # extract the finalized posterior density
  post_dist <- estep$post_dist

  # delete 'estep' and 'mstep' object
  rm(estep, mstep, data_drm, data_prm, data_all, envir = environment(), inherits = FALSE)

  # compute the information matrix of item parameter estimates using the cross-product method
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
  # holds the item parameter estimates (see est_irt_em())
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

  # divide the parameter estimation results to each group
  full_all_df_gr <-
    purrr::map(
      .x = item.id,
      .f = function(x) {
        df.tmp <-
          full_all_df %>%
          dplyr::filter(.data$id %in% x)
        df.tmp <- df.tmp[match(x, df.tmp$id), ]
        rownames(df.tmp) <- 1:nrow(df.tmp)
        df.tmp
      }
    )
  names(full_all_df_gr) <- group.name
  par_df_gr <-
    purrr::map(
      .x = item.id,
      .f = function(x) {
        df.tmp <-
          par_df %>%
          dplyr::filter(.data$id %in% x)
        df.tmp <- df.tmp[match(x, df.tmp$id), ]
        rownames(df.tmp) <- 1:nrow(df.tmp)
        df.tmp
      }
    )
  names(par_df_gr) <- group.name
  se_df_gr <-
    purrr::map(
      .x = item.id,
      .f = function(x) {
        df.tmp <-
          se_df %>%
          dplyr::filter(.data$id %in% x)
        df.tmp <- df.tmp[match(x, df.tmp$id), ]
        rownames(df.tmp) <- 1:nrow(df.tmp)
        df.tmp
      }
    )
  names(se_df_gr) <- group.name

  # population density parameters
  group.par <- vector("list", ngroup)
  for (i in 1:ngroup) {
    if (i %in% ref.group) {
      moments.se <- rep(NA_real_, 3)
      group.par.tmp <- data.frame(rbind(pop_moments[[i]], moments.se))
    } else {
      moments.est <- pop_moments[[i]]
      se.mu <- moments.est[3] / sqrt(nstd.gr[i])
      se.sigma2 <- moments.est[2] * sqrt(2 / (nstd.gr[i] - 1))
      se.sigma <- (1 / (2 * moments.est[3])) * se.sigma2 # using Delta method
      moments.se <- c(se.mu, se.sigma2, se.sigma)
      group.par.tmp <- data.frame(rbind(moments.est, moments.se))
    }
    colnames(group.par.tmp) <- c("mu", "sigma2", "sigma")
    rownames(group.par.tmp) <- c("estimates", "se")
    group.par[[i]] <- group.par.tmp
  }
  names(group.par) <- group.name

  # prior information
  if (use.aprior) aprior.dist <- aprior else aprior.dist <- NULL
  if (use.bprior) bprior.dist <- bprior else bprior.dist <- NULL
  if (use.gprior) gprior.dist <- gprior else gprior.dist <- NULL

  # count the free parameters: the item parameters, the mean and variance of each free
  # group with a normal distribution, and the weights of each empirical histogram, whose
  # number of free values is the number of quadrature points minus one for a free group
  # (sum to one) and minus three for a reference group (sum to one, mean, and variance fixed)
  npar.item <- length(param_loc$reloc.par)
  npar.dist <- if (EmpHist) 0L else length(free.group) * 2L
  npar.weight <- if (EmpHist) length(free.group) * (n.quad - 1L) + length(ref.group) * (n.quad - 3L) else 0L
  npar.est <- npar.item + npar.dist + npar.weight
  npar.detail <- c(item = npar.item, dist = npar.dist, weight = npar.weight)
  neg2llke <- -2 * llike
  aic <- 2 * npar.est + neg2llke
  bic <- npar.est * log(nstd) + neg2llke

  ## ---------------------------------------------------------------
  # check end time
  end.time <- Sys.time()

  # record total computation time
  est_time3 <- round(as.numeric(difftime(end.time, start.time, units = "secs")), 2)

  rst <- list(
    estimates = list(overall = full_all_df, group = full_all_df_gr),
    par.est = list(overall = par_df, group = par_df_gr),
    se.est = list(overall = se_df, group = se_df_gr),
    pos.par = loc_df, covariance = cov_mat, loglikelihood = list(overall = llike, group = llike.gr), aic = aic, bic = bic,
    group.par = group.par, weights = weights.gr, posterior.dist = post_dist, data = list(overall = data, group = data.gr),
    scale.D = D, ncase = list(overall = nstd, group = nstd.gr), nitem = list(overall = nitem, group = nitem.gr),
    Etol = Etol, MaxE = MaxE, aprior = aprior.dist, bprior = bprior.dist, gprior = gprior.dist, npar.est = npar.est,
    niter = r, maxpar.diff = max.diff, EMtime = est_time1, SEtime = est_time2, TotalTime = est_time3, test.1 = memo3,
    test.2 = memo4, var.note = memo5, fipc = FALSE, fipc.method = NULL, fix.loc = list(overall = NULL, group = NULL),
    npar.detail = npar.detail
  )

  if (verbose) {
    cat("Estimation is finished in", est_time3, "seconds.", "\n")
  }
  return(rst)
}


# multiple group FIPC via EM
#' @import dplyr
est_mg_fipc <- function(x = NULL,
                        data,
                        group.name = NULL,
                        D = 1,
                        model = NULL,
                        cats = NULL,
                        item.id = NULL,
                        free.group = NULL,
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
                        Etol = 1e-03,
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

  # start parsing inputs
  if (verbose) {
    cat("Parsing input...", "\n")
  }

  # count the number of groups
  ngroup <- length(data)
  if (ngroup == 1) {
    stop("Use the est_irt() function when a single-group item calibration is implemented.", call. = FALSE)
  }

  # create a group name vector when group.name = NULL
  if (is.null(group.name)) {
    group.name <- paste0("g", 1:ngroup)
  }

  # test length of each group
  nitem.gr <- purrr::map_dbl(.x = data, ncol)
  names(nitem.gr) <- group.name

  # sample size of each group
  nstd.gr <- purrr::map_dbl(.x = data, nrow)
  names(nstd.gr) <- group.name

  # confirm and correct all item metadata information
  x.gr <- warn_drm_once(purrr::map(.x = x, ~ {
    confirm_df(.x)
  }))
  names(x.gr) <- group.name

  # override the item IDs in x when item.id is supplied
  if (!is.null(item.id)) {
    x.gr <- purrr::map2(.x = x.gr, .y = item.id, ~ {
      .x$id <- .y
      .x
    })
  }

  # item, cats, model information for each group
  item.id <- purrr::map(.x = x.gr, ~ {
    .x$id
  })
  cats.gr <- purrr::map(.x = x.gr, ~ {
    .x$cats
  })
  model.gr <- purrr::map(.x = x.gr, ~ {
    .x$model
  })

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
    # stop when an ID in 'fix.id' is not found in the item IDs of any group
    id.notfound <- setdiff(unique(fix.id), unlist(item.id))
    if (length(id.notfound) > 0L) {
      stop("The item ID(s) in 'fix.id' not found in the item IDs: ", paste(id.notfound, collapse = ", "), ".",
           call. = FALSE)
    }
    fix.item <- unique(fix.id)
    fix.item.gr <- purrr::map(.x = item.id, ~ {
      .x[.x %in% fix.item]
    })
    fix.loc.gr <- purrr::map(.x = item.id, ~ {
      which(.x %in% fix.item)
    })
    names(fix.loc.gr) <- group.name
  } else {
    # stop when 'fix.loc' is not a list with one element (NULL or row positions) per group
    fix.loc.chk <- fix.loc
    if (!is.list(fix.loc.chk) || length(fix.loc.chk) != ngroup) {
      stop("'fix.loc' must be a list with one element for each of the ", ngroup, " groups.", call. = FALSE)
    }
    # groups whose positions are missing, outside the rows of the group's metadata, or repeated
    bad.loc <- which(vapply(seq_len(ngroup), function(g) {
      l <- fix.loc.chk[[g]]
      length(l) > 0L && (anyNA(l) || any(!l %in% seq_len(nitem.gr[g])) || anyDuplicated(l) > 0L)
    }, logical(1)))
    if (length(bad.loc) > 0L) {
      stop("Each element of 'fix.loc' must contain distinct row positions of the item metadata of its group. ",
           "Check group(s) ", paste(bad.loc, collapse = ", "), ".", call. = FALSE)
    }
    # stop when no item is fixed
    if (sum(lengths(fix.loc.chk)) == 0L) {
      stop("'fix.loc' must contain at least one fixed item.", call. = FALSE)
    }
    fix.item.gr <- purrr::map2(.x = item.id, .y = fix.loc, ~ {
      .x[.y]
    })
    fix.item <- unique(unlist(fix.item.gr))
    fix.loc.gr <- fix.loc
    names(fix.loc.gr) <- group.name
  }

  # stop when an item ID is repeated within a group or is defined differently across groups,
  # and warn when a fixed item has different parameter values across groups
  check_mg_items(x.gr, fix.item = fix.item)

  # copy data to data.gr object
  data.gr <- data
  names(data.gr) <- group.name

  # create combined response data sets across all groups
  # and assign item IDs
  data <-
    # read the responses of each group by their values before the binding
    purrr::map(.x = data.gr, ~ {
      # recode a non-numeric missing code before the responses are read as numbers
      if (!is.na(missing)) .x[.x == missing] <- NA
      data.frame(resp_as_numeric(.x))
    }) %>%
    purrr::map2(
      .y = item.id,
      .f = function(x, y) {
        dplyr::rename_all(.tbl = x, .funs = function(x) y)
      }
    ) %>%
    dplyr::bind_rows()

  # a vector of item ID for the combined data set
  id <- colnames(data)

  # find the locations of items that should be fixed in the combined data set
  fix.loc <- which(id %in% fix.item)

  # create an item metadata for the combined data set
  x <-
    dplyr::bind_rows(x.gr) %>%
    dplyr::distinct(.data$id, .keep_all = TRUE)
  if (!all(id == x$id)) {
    x <- x[match(id, x$id), ]
  }
  rownames(x) <- 1:nrow(x)

  # recode missing values
  if (!is.na(missing)) {
    data[data == missing] <- NA
    # recode the missing values in the response data of each group as well
    data.gr <- purrr::map(.x = data.gr, ~ {
      .x[.x == missing] <- NA
      .x
    })
  }

  # transform a data set to matrix (character and factor columns are read by their values)
  data <- resp_as_numeric(data)

  # check the total number of examinees for the combined data
  nstd <- nrow(data)

  # check the total number of items in the combined data
  nitem.all <- ncol(data)

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

  # clean the two item metadata sets
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

  # designate a reference group whose group parameters will be fixed
  if (is.null(free.group)) {
    stop("The free.group argument should not be NULL when FIPC = TRUE.", call. = FALSE)
  }
  # (stop when a free group does not exist)
  grp <- set_mg_groups(free.group = free.group, group.name = group.name, fipc = TRUE)
  free.group <- grp$free.group
  ref.group <- grp$ref.group

  # find the groups that have no common items with other groups
  pair <-
    t(utils::combn(group.name, m = 2)) %>%
    data.frame()
  overlap <- c()
  for (i in 1:nrow(pair)) {
    freq.tmp <-
      unlist(item.id[unlist(pair[i, ])]) %>%
      table()
    overlap[i] <- any(freq.tmp > 1)
  }
  link.gr <-
    dplyr::mutate(pair, overlap = overlap) %>%
    dplyr::filter(.data$overlap == TRUE) %>%
    dplyr::select(1:2) %>%
    unlist() %>%
    unique()
  nolink.gr <- c(1:ngroup)[!group.name %in% link.gr]

  # check if there exist free groups which do not share common items with other groups
  nolink.free <- free.group[free.group %in% nolink.gr]

  # count the number of fixed items in each group (an item fixed in one group is fixed in all groups)
  n.fix.gr <- purrr::map_dbl(.x = item.id, ~ {
    sum(.x %in% fix.item)
  })

  # when there exist groups that are specified as free groups but not have fixed items,
  # stop the further process
  if (!all(nolink.free %in% which(n.fix.gr != 0))) {
    warning.memo <-
      paste0(
        # names of the free groups that have neither common items nor fixed items
        paste(group.name[nolink.free[!nolink.free %in% which(n.fix.gr != 0)]], collapse = ", "),
        " group(s) in which ability distribution(s) is(are) freely estimated do(does) not share common items with other groups \n",
        "nor include fixed items. Please specify the freely estimated group(s) correctly."
      )
    stop(paste0(warning.memo, " \n"), call. = FALSE)
  }

  # check the number of item responses across all items
  n.resp <- Rfast::colsums(!is.na(data))

  # check the items which have all missing responses
  loc_allmiss <- which(n.resp == 0L)
  if (length(loc_allmiss) > 0L) {
    # IDs of the columns of the combined data set (id now holds the IDs of the new items only)
    allmiss_item <- x$id[loc_allmiss]
    memo2 <- paste0(paste0(allmiss_item, collapse = ", "), " has/have no item response data across all groups. \n")
    stop(memo2, call. = FALSE)
  }

  # stop when a response is not one of the score categories 0, ..., cats - 1
  check_resp_range(data = data, cats = x$cats, id = x$id)

  # stop when the responses to all fixed items show no variation (pooled over all groups)
  check_fix_info(data = data[, fix.loc, drop = FALSE], id = x$id[fix.loc])

  # warn about the items to be estimated whose maximum likelihood estimates do not exist,
  # using the responses pooled over all groups
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
  # create initial weights of prior ability distribution when it is not specified
  if (is.null(weights)) {
    # create quadrature points
    quadpt <- seq(-Quadrature[2], Quadrature[2], length.out = Quadrature[1])

    # create the data.frame containing the quadrature points and weights
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

  # a list containing the weights and densities for each group
  weights.gr <- replicate(n = ngroup, expr = weights, simplify = FALSE)
  names(weights.gr) <- group.name

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
  if (!is.null(x_new)) {
    # (no item other than the constrained 1PLM items when loc_else is empty)
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
    mmt_dist_old <- matrix(rep(c(group.mean, group.var), ngroup), ncol = ngroup)
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

  # create indices of samples (examinees) who belong to each group
  idx.tmp2 <- as.numeric(cumsum(nstd.gr))
  idx.tmp1 <- dplyr::lag(idx.tmp2, default = 0) + 1
  idx.std <- purrr::map2(.x = idx.tmp1, .y = idx.tmp2, ~ {
    .x:.y
  })
  names(idx.std) <- group.name

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
          data_all1 = data.new, weights = weights.gr, D = D, idx.std = idx.std
        )
      } else {
        estep <- Estep_fipc(
          elm_item1 = elm_item_new, elm_item2 = elm_item_all, idx.drm2 = idx.drm.all,
          idx.prm2 = idx.prm.all, data_drm2 = data_drm.all, data_prm2 = data_prm.all,
          data_all1 = data.new, weights = weights.gr, D = D, idx.std = idx.std
        )
      }
    } else {
      estep <- Estep_fipc(
        elm_item1 = elm_item_all, elm_item2 = elm_item_all, idx.drm2 = idx.drm.all,
        idx.prm2 = idx.prm.all, data_drm2 = data_drm.all, data_prm2 = data_prm.all,
        data_all1 = data.all, weights = weights.gr, D = D, idx.std = idx.std
      )
      estep$elm_item <- NULL
    }

    # implement M-step
    mstep <- Mstep(
      estep = estep, id = id, cats = cats, model = model, quadpt = quadpt, n.quad = n.quad,
      D = D, cols.item = cols.item, loc_1p_const = loc_1p_const, loc_else = loc_else, idx4est = idx4est,
      n.1PLM = n.1PLM, EmpHist = EmpHist, weights = weights.gr, fix.a.1pl = fix.a.1pl, fix.a.gpcm = fix.a.gpcm, fix.g = fix.g,
      a.val.1pl = a.val.1pl, a.val.gpcm = a.val.gpcm, g.val = g.val, use.aprior = use.aprior, use.bprior = use.bprior,
      use.gprior = use.gprior, aprior = aprior, bprior = bprior, gprior = gprior, group.mean = group.mean,
      group.var = group.var, nstd = nstd.gr, control = control,
      iter = r, fipc = TRUE, reloc.par = param_loc$reloc.par, ref.group = ref.group,
      free.group = free.group, parbd = parbd
    )

    if (!is.null(x_new)) {
      # compute the difference between previous and updated item parameter estimates
      diff_par <- mstep$elm_item$pars - elm_item_new$pars
      # convergence statistic: largest absolute parameter change
      max.diff <- suppressWarnings(max(abs(diff_par), na.rm = TRUE))
    } else {
      # compute the mean and variance of the updated prior distribution
      mmt_dist_new <-
        as.matrix(purrr::map_dfc(.x = mstep$weights, ~ {
          cal_moment(node = .x$theta, weight = .x$weight)
        }))
      diff_par <- mmt_dist_new - mmt_dist_old
      # convergence statistic: largest absolute change in prior mean and variance
      max.diff <- suppressWarnings(max(abs(diff_par), na.rm = TRUE))
    }

    # loglikelihood value
    llike <- do.call(what = "sum", args = mstep$loglike)

    # print
    if (verbose) {
      cat("\r", paste0(
        "EM iteration: ", r, ", Loglike: ",
        format(round(llike, 4), nsmall = 4), ", Max-Change: ",
        format(round(max.diff, 6), nsmall = 5)
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
    weights.gr <- mstep$weights

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

  # conduct one more E-step to update the posterior distribution using the final item parameter estimates
  if (!is.null(x_new)) {
    estep <- Estep_fipc(
      elm_item1 = elm_item_new, elm_item2 = elm_item_all, idx.drm2 = idx.drm.all,
      idx.prm2 = idx.prm.all, data_drm2 = data_drm.all, data_prm2 = data_prm.all,
      data_all1 = data.new, weights = weights.gr, D = D, idx.std = idx.std
    )
  } else {
    estep <- Estep_fipc(
      elm_item1 = elm_item_all, elm_item2 = elm_item_all, idx.drm2 = idx.drm.all,
      idx.prm2 = idx.prm.all, data_drm2 = data_drm.all, data_prm2 = data_prm.all,
      data_all1 = data.all, weights = weights.gr, D = D, idx.std = idx.std
    )
  }

  # compute the final log of marginal likelihood
  likehd.gr <- purrr::map(.x = idx.std, ~ {
    estep$likehd[.x, ]
  })

  # log-scale offsets of the rows rescaled in likelihood(), by group
  offset.gr <- purrr::map(.x = idx.std, ~ {
    estep$offset[.x]
  })
  llike.gr <- purrr::pmap(
    .l = list(likehd.gr, weights.gr, offset.gr),
    .f = function(lk, wt, os) {
      sum(log(lk %*% matrix(wt[, 2]))) + sum(os)
    }
  )
  llike <- do.call(what = "sum", args = llike.gr)

  # compute the mean and variance of the estimated density distributions
  pop_moments <-
    purrr::map(.x = weights.gr, ~ {
      cal_moment(node = quadpt, weight = .x[, 2])
    }) %>%
    purrr::map(~ {
      c(.x, sigma = sqrt(as.numeric(.x[2])))
    })
  if (length(ref.group) > 0) {
    pop_moments[ref.group] <-
      purrr::map(
        .x = 1:length(ref.group),
        ~ {
          c(mu = group.mean, sigma2 = group.var, sigma = sqrt(group.var))
        }
      )
  }

  ## ---------------------------------------------------------------------
  # estimates the information matrix and standard errors of item parameter estimates
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

    # compute the information matrix of item parameter estimates using the cross-product method
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
      # see est_mg_em()
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
    # 1) for the only new items (see est_irt_em() for the loc.par layout).
    se_df <- loc.par <- param_loc$loc.par
    mask  <- !is.na(loc.par)
    if (se) {
      se_df[mask] <- se_par[loc.par[mask]]
    } else {
      se_df[mask] <- NA_real_
    }

    # 2) for all items of the combined form
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

    # deploy the standard errors on the location of matrix as the item parameter estimates
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

  # divide the parameter estimation results to each group
  full_all_df_gr <-
    purrr::map(
      .x = item.id,
      .f = function(x) {
        df.tmp <-
          full_all_df %>%
          dplyr::filter(.data$id %in% x)
        df.tmp <- df.tmp[match(x, df.tmp$id), ]
        rownames(df.tmp) <- 1:nrow(df.tmp)
        df.tmp
      }
    )
  names(full_all_df_gr) <- group.name
  par_df_gr <-
    purrr::map(
      .x = item.id,
      .f = function(x) {
        df.tmp <-
          x_all %>%
          dplyr::filter(.data$id %in% x)
        df.tmp <- df.tmp[match(x, df.tmp$id), ]
        rownames(df.tmp) <- 1:nrow(df.tmp)
        df.tmp
      }
    )
  names(par_df_gr) <- group.name
  se_df_gr <-
    purrr::map(
      .x = item.id,
      .f = function(x) {
        df.tmp <-
          se_all_df %>%
          dplyr::filter(.data$id %in% x)
        df.tmp <- df.tmp[match(x, df.tmp$id), ]
        rownames(df.tmp) <- 1:nrow(df.tmp)
        df.tmp
      }
    )
  names(se_df_gr) <- group.name

  # population density parameters
  group.par <- vector("list", ngroup)
  for (i in 1:ngroup) {
    if (i %in% ref.group) {
      moments.se <- rep(NA_real_, 3)
      group.par.tmp <- data.frame(rbind(pop_moments[[i]], moments.se))
    } else {
      moments.est <- pop_moments[[i]]
      se.mu <- moments.est[3] / sqrt(nstd.gr[i])
      se.sigma2 <- moments.est[2] * sqrt(2 / (nstd.gr[i] - 1))
      se.sigma <- (1 / (2 * moments.est[3])) * se.sigma2 # using Delta method
      moments.se <- c(se.mu, se.sigma2, se.sigma)
      group.par.tmp <- data.frame(rbind(moments.est, moments.se))
    }
    colnames(group.par.tmp) <- c("mu", "sigma2", "sigma")
    rownames(group.par.tmp) <- c("estimates", "se")
    group.par[[i]] <- group.par.tmp
  }
  names(group.par) <- group.name

  # prior information
  if (use.aprior) aprior.dist <- aprior else aprior.dist <- NULL
  if (use.bprior) bprior.dist <- bprior else bprior.dist <- NULL
  if (use.gprior) gprior.dist <- gprior else gprior.dist <- NULL

  # count the free parameters: the item parameters to be estimated, the mean and variance of
  # each free group with a normal distribution, and the weights of each empirical histogram,
  # whose number of free values is the number of quadrature points minus one for a free group
  # (sum to one) and minus three for a reference group (sum to one, mean, and variance fixed)
  npar.item <- if (!is.null(x_new)) length(param_loc$reloc.par) else 0L
  npar.dist <- if (EmpHist) 0L else length(free.group) * 2L
  npar.weight <- if (EmpHist) length(free.group) * (n.quad - 1L) + length(ref.group) * (n.quad - 3L) else 0L
  npar.est <- npar.item + npar.dist + npar.weight
  npar.detail <- c(item = npar.item, dist = npar.dist, weight = npar.weight)
  neg2llke <- -2 * llike
  aic <- 2 * npar.est + neg2llke
  bic <- npar.est * log(nstd) + neg2llke

  ## ---------------------------------------------------------------
  # check end time
  end.time <- Sys.time()

  # record total computation time
  est_time3 <- round(as.numeric(difftime(end.time, start.time, units = "secs")), 2)

  rst <- list(
    estimates = list(overall = full_all_df, group = full_all_df_gr),
    par.est = list(overall = x_all, group = par_df_gr),
    se.est = list(overall = se_all_df, group = se_df_gr),
    pos.par = loc_all_df, covariance = cov_mat, loglikelihood = list(overall = llike, group = llike.gr), aic = aic, bic = bic,
    group.par = group.par, weights = weights.gr, posterior.dist = post_dist, data = list(overall = data_all, group = data.gr),
    scale.D = D, ncase = list(overall = nstd, group = nstd.gr), nitem = list(overall = nitem.all, group = nitem.gr),
    Etol = Etol, MaxE = MaxE, aprior = aprior.dist, bprior = bprior.dist, gprior = gprior.dist, npar.est = npar.est,
    niter = r, maxpar.diff = max.diff, EMtime = est_time1, SEtime = est_time2, TotalTime = est_time3, test.1 = memo3,
    test.2 = memo4, var.note = memo5, fipc = TRUE, fipc.method = fipc.method, fix.loc = list(overall = fix.loc, group = fix.loc.gr),
    npar.detail = npar.detail
  )

  if (verbose) {
    cat("Estimation is finished in", est_time3, "seconds.", "\n")
  }
  return(rst)
}


# This function stops when the arguments of est_mg() that are given per group do not
# have one element for each group in 'data'
check_mg_args <- function(x, data, group.name, model, cats, item.id, fipc, fipc.method) {
  # 'data' must be a list with one response data set per group
  if (!is.list(data) || is.data.frame(data)) {
    stop("'data' must be a list containing one response data matrix for each group.", call. = FALSE)
  }

  # number of groups and number of items (columns) in each group's response data
  ngroup <- length(data)
  nitem.gr <- vapply(data, NCOL, numeric(1))

  # one unique group name per group
  if (!is.null(group.name) && (length(group.name) != ngroup || anyDuplicated(group.name) > 0L)) {
    stop("'group.name' must contain ", ngroup, " unique names, one for each element of 'data'.", call. = FALSE)
  }

  # the FIPC method must be one of the two available methods
  if (fipc && !(length(fipc.method) == 1L && fipc.method %in% c("MEM", "OEM"))) {
    stop("'fipc.method' must be either \"MEM\" or \"OEM\".", call. = FALSE)
  }

  # the item metadata is a list with one data frame per group, one row per response column
  if (!is.null(x)) {
    if (!is.list(x) || is.data.frame(x) || length(x) != ngroup) {
      stop("'x' must be a list containing the item metadata of each of the ", ngroup, " groups.", call. = FALSE)
    }
    # groups whose item metadata and response data have different numbers of items
    bad.x <- which(vapply(x, NROW, numeric(1)) != nitem.gr)
    if (length(bad.x) > 0L) {
      stop("The number of items in 'x' and 'data' must be the same in each group. ",
           "They differ in group(s) ", paste(bad.x, collapse = ", "), ".", call. = FALSE)
    }
  } else if (!fipc) {
    # model and cats are recycled across the groups only when a single element is given
    for (arg in list(list(model, "model"), list(cats, "cats"))) {
      if (!is.null(arg[[1]]) && !length(arg[[1]]) %in% c(1L, ngroup)) {
        stop("'", arg[[2]], "' must be a list with one element for each of the ", ngroup, " groups.", call. = FALSE)
      }
    }
  }

  # the item IDs are a list with one ID per item of each group
  if (!is.null(item.id)) {
    if (!is.list(item.id) || length(item.id) != ngroup) {
      stop("'item.id' must be a list containing the item IDs of each of the ", ngroup, " groups.", call. = FALSE)
    }
    # groups whose item IDs do not have one value per item
    bad.id <- which(lengths(item.id) != nitem.gr)
    if (length(bad.id) > 0L) {
      stop("Each element of 'item.id' must have one ID for each item (column) of the group's response data. ",
           "Check group(s) ", paste(bad.id, collapse = ", "), ".", call. = FALSE)
    }
  }

  invisible(NULL)
}


# This function converts 'free.group' to group indices and finds the reference groups whose
# latent distributions are fixed; it stops when 'free.group' names a group that does not
# exist or, without FIPC, leaves no reference group to fix the scale
set_mg_groups <- function(free.group, group.name, fipc) {
  # number of groups
  ngroup <- length(group.name)

  if (is.null(free.group)) {
    # all groups are reference groups
    ref.group <- 1:ngroup
  } else {
    if (is.character(free.group)) {
      # stop when a name in 'free.group' is not one of the group names
      if (!all(free.group %in% group.name)) {
        stop("The group(s) ", paste(setdiff(free.group, group.name), collapse = ", "),
             " in 'free.group' do(does) not exist in the group.name argument.", call. = FALSE)
      }
      free.group <- which(group.name %in% free.group)
    } else {
      # stop when a value in 'free.group' is not one of the group indices
      if (!all(free.group %in% 1:ngroup)) {
        stop("The values in 'free.group' must be group indices between 1 and ", ngroup, ".", call. = FALSE)
      }
      # count each free group once
      free.group <- unique(free.group)
    }
    ref.group <- c(1:ngroup)[-free.group]
  }

  # without fixed items, at least one reference group is needed to fix the scale
  if (!fipc && length(ref.group) == 0L) {
    stop("At least one group must be left out of 'free.group' so that its latent distribution ",
         "fixes the scale when fipc = FALSE.", call. = FALSE)
  }

  list(free.group = free.group, ref.group = ref.group)
}


# This function evaluates an expression and gives the warning about 'DRM' items of
# confirm_df() once, however many groups raise it (other warnings pass through)
warn_drm_once <- function(expr) {
  drm.memo <- NULL
  out <- withCallingHandlers(expr, warning = function(w) {
    # hold back the 'DRM' warning
    if (grepl("'DRM' items are treated as '3PLM' items", conditionMessage(w), fixed = TRUE)) {
      drm.memo <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    }
  })
  # give the held warning once
  if (!is.null(drm.memo)) {
    warning(drm.memo, call. = FALSE)
  }
  out
}


# This function stops when an item ID is repeated within a group or has different numbers
# of score categories or different models across groups, and warns when a fixed item of
# FIPC has different parameter values across groups (those of the first group are used)
check_mg_items <- function(x.gr, fix.item = NULL) {
  # item IDs repeated within a group
  dup.id <- unique(unlist(lapply(x.gr, function(z) z$id[duplicated(z$id)])))
  if (length(dup.id) > 0L) {
    stop("Item IDs must be unique within each group. Repeated item ID(s): ",
         paste(dup.id, collapse = ", "), ".", call. = FALSE)
  }

  # stack the item metadata of all groups
  x.all <- dplyr::bind_rows(x.gr)

  # item IDs whose numbers of score categories or models differ across groups
  n.def <- tapply(paste(x.all$cats, x.all$model), x.all$id, function(k) length(unique(k)))
  bad.id <- names(n.def)[n.def > 1L]
  if (length(bad.id) > 0L) {
    stop("Items with the same ID must have the same number of score categories and the same model ",
         "in all groups. Check item(s) ", paste(bad.id, collapse = ", "), ".", call. = FALSE)
  }

  # fixed items whose parameter values differ across groups
  if (length(fix.item) > 0L) {
    x.fix <- x.all[x.all$id %in% fix.item, , drop = FALSE]
    # number of parameters of each item: a and b (and g for 3PLM), or a and cats - 1 thresholds
    npar <- ifelse(x.fix$cats == 2, ifelse(x.fix$model == "3PLM", 3, 2), x.fix$cats)
    # parameter values of each item row
    par.val <- lapply(seq_len(nrow(x.fix)), function(i) {
      as.numeric(unlist(x.fix[i, 3 + seq_len(npar[i])]))
    })
    # items whose values differ across groups by more than a numerical tolerance
    id.all <- sort(unique(x.fix$id))
    differ <- vapply(id.all, function(id) {
      loc <- which(x.fix$id == id)
      ref <- par.val[[loc[1]]]
      any(vapply(par.val[loc[-1]], function(p) {
        d <- abs(p - ref)
        # two missing values are equal
        d[is.na(p) & is.na(ref)] <- 0
        anyNA(d) || max(d) > 1e-10
      }, logical(1)))
    }, logical(1))
    diff.id <- id.all[differ]
    if (length(diff.id) > 0L) {
      warning("The fixed item(s) ", paste(diff.id, collapse = ", "), " have different parameter values ",
              "across groups; the values in the first group that includes each item were used.", call. = FALSE)
    }
  }

  invisible(NULL)
}
