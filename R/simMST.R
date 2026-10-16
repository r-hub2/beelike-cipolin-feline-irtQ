#' Simulated 1-3-3 MST Panel Data
#'
#' A simulated Multistage-Adaptive Test (MST) dataset representing a 1-3-3 panel
#' structure, assembled using module-level target test information functions
#' (TIFs) and design constraints similar to those used in the simulation study
#' by Lim et al. (2021).
#'
#' @usage simMST
#'
#' @format A list containing five internal objects:
#' \describe{
#'   \item{item_bank}{A data frame of item metadata for the 56 items assigned
#'   to the panel (one row per item; columns follow the \pkg{irtQ} item
#'   metadata format: \code{id}, \code{cats}, \code{model}, \code{par.1},
#'   \code{par.2}, \code{par.3}).}
#'
#'   \item{module}{A binary integer matrix (56 x 7) that maps items in the
#'   item bank to MST modules. This structure specifies the item-to-module
#'   assignment, similar to the `modules` argument in the `randomMST()`
#'   function from the \pkg{mstR} package (Magis et al., 2017).}
#'
#'   \item{route_map}{A binary square integer matrix (7 x 7) that defines the
#'   MST transition structure, showing module pathways across stages. Module 1
#'   routes to modules 2 to 4, module 2 to modules 5 and 6, module 3 to
#'   modules 5 to 7, and module 4 to modules 6 and 7. This corresponds to the
#'   `transMatrix` argument in the `randomMST()` function in the \pkg{mstR}
#'   package.}
#'
#'   \item{cut_score}{A list of two numeric vectors, `stage.2` and `stage.3`,
#'   giving the routing cut scores into stages 2 and 3. They are the
#'   TIF-crossing points in the `cut_score` element returned by
#'   `find_cut(x = simMST$item_bank, module = simMST$module,`
#'   `route_map = simMST$route_map, D = 1.702)`.}
#'
#'   \item{theta}{A numeric vector of 81 ability values from -4 to 4 in steps
#'   of 0.1, used to evaluate the measurement precision of the panel across the
#'   latent trait continuum.}
#' }
#'
#'   This 1-3-3 MST panel includes 7 modules across 3 stages, drawn from
#'   56 unique items (8 per module) with no item shared across modules. Each
#'   module contains 8 dichotomously scored items calibrated under the IRT
#'   3-parameter logistic (3PL) model. The item parameters are on the D = 1.702
#'   scale, so pass \code{D = 1.702} to \code{\link{find_cut}()},
#'   \code{\link{reval_mst}()}, and \code{\link{run_mst}()} when using this
#'   panel.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @references Lim, H., Davey, T., & Wells, C. S. (2021). A recursion-based
#'   analytical approach to evaluate the performance of MST. *Journal of
#'   Educational Measurement, 58*(2), 154-178. \doi{10.1111/jedm.12276}.
#'
#'   Magis, D., Yan, D., & von Davier, A. A. (2017). *Computerized adaptive and
#'   multistage testing with R: Using packages catR and mstR*. Springer.
#'   \doi{10.1007/978-3-319-69218-0}.
#'
#' @examples
#' # structure of the data
#' str(simMST, max.level = 1)
#'
#' # item metadata of the first three items
#' head(simMST$item_bank, 3)
#'
"simMST"
