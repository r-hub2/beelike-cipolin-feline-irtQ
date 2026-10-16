#' Recursion-Based Analytical Evaluation of an MST Panel
#'
#' This function evaluates the measurement precision and bias in
#' Multistage-Adaptive Test (MST) panels using a recursion-based evaluation
#' method introduced by Lim et al. (2021). This function computes conditional
#' biases and standard errors of measurement (CSEMs) across a range of IRT
#' ability levels, facilitating efficient and accurate MST panel assessments
#' without extensive simulations.
#'
#' @param x A data frame containing the item metadata of the item bank (e.g.,
#'   item parameters, number of categories, and IRT models). Each row is one
#'   item, and the rows correspond to the rows of `module`. See
#'   [irtQ::shape_df()] for creating item metadata.
#' @param D A scaling constant used in IRT models to make the logistic function
#'   closely approximate the normal ogive function. A value of 1.702 is commonly
#'   used for this purpose. Default is 1.
#' @param route_map A binary square matrix defining the MST structure. A value
#'   of 1 at row *i* and column *j* means that test takers can be routed from
#'   module *i* to module *j*. This is the same convention as the `transMatrix`
#'   argument of `randomMST()` in the \pkg{mstR} package (Magis et al., 2017).
#'   The stages are found from the transitions (see [irtQ::panel_info()]), and
#'   the first stage must contain a single routing module.
#' @param module A binary matrix that assigns the items in `x` to modules. It
#'   has one row per item (in the order of `x`) and one column per module, and a
#'   value of 1 at row *j* and column *m* means that item *j* belongs to module
#'   *m*. This is the same convention as the `modules` argument of `randomMST()`
#'   in the \pkg{mstR} package (Magis et al., 2017).
#' @param cut_score A list with one numeric vector per stage transition (the
#'   number of stages minus one). Element *s* holds the cut scores on the
#'   ability (theta) metric for routing from stage *s* to stage *s* + 1, finite
#'   and in strictly ascending order, with one fewer value than the number of
#'   modules in stage *s* + 1; cut score *k* separates the *k*-th and the
#'   (*k* + 1)-th module of that stage in the order of the module indices. A
#'   `NULL` element is read as an empty vector for a stage with a single
#'   module. For example, in a 1-3-3 panel,
#'   `cut_score = list(c(-0.5, 0.5), c(-0.6, 0.6))`. See **Details**.
#' @param theta A vector of ability levels (theta) at which the MST panel's
#'   performance is assessed. This allows for the evaluation of measurement
#'   precision and bias across a continuum of ability levels. The default is
#'   `theta = seq(-5, 5, 1)`.
#' @param intpol Logical. If `TRUE`, inverse TCC scoring assigns ability
#'   estimates to the sum scores that cannot be mapped through the TCC (sum
#'   scores less than or equal to the sum of the guessing parameters, and the
#'   maximum possible sum score), as in [irtQ::est_score()]. The low scores are
#'   mapped by linear interpolation (Lim et al., 2021), and the maximum sum
#'   score receives the second value of `range.tcc`. With `intpol = FALSE`,
#'   some sum scores have no ability estimate, and the function stops with an
#'   error. Default is `TRUE`.
#' @param range.tcc A numeric vector of length two giving the ability estimates
#'   assigned to the lowest and the maximum possible sum scores when
#'   `intpol = TRUE` (see [irtQ::est_score()]). The range must be wide enough
#'   to cover the ability estimates of the sum scores, or the function stops
#'   with an error. Default is `c(-7, 7)`.
#' @param tol A positive number giving the tolerance of the bisection search
#'   used by inverse TCC scoring; the search stops when the search interval is
#'   no wider than `tol`. Default is 1e-4.
#'
#' @details The [irtQ::reval_mst()] function evaluates an MST panel by
#'   implementing a recursion-based method to assess measurement precision
#'   across IRT ability levels. This approach, detailed in Lim et al. (2021),
#'   enables the computation of conditional biases and CSEMs efficiently,
#'   bypassing the need for extensive simulations traditionally required for MST
#'   evaluation.
#'
#'   The recursion of Lim et al. (2021) is built on inverse TCC ability
#'   estimates. At each stage, the sum score accumulated over all modules
#'   administered so far is converted to an ability estimate by inverse TCC
#'   scoring, this estimate is compared with the cut scores to route the test
#'   taker to the next module, and the final ability estimate is also the
#'   inverse TCC estimate of the total sum score. Accordingly, the function
#'   supports inverse TCC scoring with cut-score routing only. To evaluate
#'   other scoring or routing methods, use [irtQ::run_mst()], which runs a
#'   Monte Carlo simulation. With `route_method = NULL`, a `cut_score` list, and
#'   `route_score = list(method = "INV.TCC")`, [irtQ::run_mst()] follows the
#'   same design that this function evaluates analytically.
#'
#'   The first stage must have a single routing module, whatever its module
#'   index, and the function stops with an error when stage 1 has more than one
#'   module. [irtQ::run_mst()] also runs panels with several stage-1 modules, so
#'   the two functions differ in the panels they support.
#'
#'   All modules in the same stage must have the same maximum sum score (the
#'   sum of the maximum item scores), for example the same number of items when
#'   all items are dichotomous. The function stops with an error when this
#'   condition is not met.
#'
#'   The `module` argument, used in conjunction with the item bank metadata `x`,
#'   systematically organizes items into modules for MST panel evaluation. Each
#'   row of `x` corresponds to an item, detailing its characteristics like score
#'   categories and IRT model. The `module` matrix, structured with the same
#'   number of rows as `x` and columns representing modules, indicates item
#'   assignments with 1s. This precise mapping enables the [irtQ::reval_mst()]
#'   function to evaluate the MST panel's performance by analyzing how items
#'   within each module contribute to measurement precision and bias, reflecting
#'   the tailored progression logic inherent in MST designs.
#'
#'   The `route_map` argument defines the possible transitions between
#'   modules, following the `transMatrix` argument of `randomMST()` in the
#'   \pkg{mstR} package (Magis et al., 2017). A value of 1 at row *i* and column
#'   *j* means that test takers can move from module *i* directly to module *j*.
#'   [irtQ::panel_info()] uses this matrix to enumerate all pathways, and
#'   [irtQ::reval_mst()] evaluates the routing along these pathways
#'   analytically, without simulation.
#'
#'   To further detail the `cut_score` argument with an illustration: In a 1-3-3
#'   MST configuration, the list `cut_score = list(c(-0.5, 0.5), c(-0.6, 0.6))`
#'   operates as a decision guide at each stage. Initially, all test takers
#'   start in the routing module. Upon completion, the inverse TCC ability
#'   estimate of the sum score over all modules taken so far determines the next
#'   module: estimates below -0.5 route to the first module of the next stage,
#'   estimates from -0.5 up to but not including 0.5 to the second, and
#'   estimates of 0.5 or above to the third. When a module can reach only some
#'   modules of the next stage, only the cut scores that separate the reachable
#'   modules are used. For example, a test taker in a stage-2 module that can
#'   reach only the second and third stage-3 modules is routed by the second
#'   cut score alone.
#'
#' @return A list with the following seven elements:
#'
#' \item{panel.info}{The output of [irtQ::panel_info()]: `config`, `pathway`,
#' `n.module`, and `n.stage`.}
#'
#' \item{item.by.mod}{A list of item metadata data frames, one per module, named
#' `m.1`, `m.2`, etc. The item metadata are checked and completed as in
#' [irtQ::run_mst()], so these are plain data frames with upper-case model names
#' and without parameter columns that are all `NA`, whatever the class of `x`.
#' The same holds for `item.by.path`.}
#'
#' \item{item.by.path}{A list with one element per stage (`stage.1`, `stage.2`,
#' ...). Each element is a list of item metadata data frames, one per distinct
#' partial pathway up to that stage (`path.1`, `path.2`, ...), holding the items
#' of all modules on the partial pathway.}
#'
#' \item{eq.theta}{A list with one element per stage. Each element is a matrix
#' with one row per possible cumulative sum score (0 to the maximum) and one
#' column per partial pathway up to that stage (`path.1`, `path.2`, ...),
#' holding the inverse TCC ability estimates. The columns of the last stage
#' correspond to the complete pathways.}
#'
#' \item{cdist.by.mod}{A list with one element per value of `theta` (named by
#' the values). Each element is a matrix with one row per sum score and one
#' column per module (`m.1`, `m.2`, ...), holding the conditional sum score
#' distribution of each module; rows beyond a module's maximum score are 0.}
#'
#' \item{jdist.by.path}{A list with one element per stage. Each element is a
#' list with one element per value of `theta`, holding a matrix with one row per
#' cumulative sum score and one column per partial pathway, which gives the
#' joint probability of the sum score and the pathway.}
#'
#' \item{eval.tb}{A data frame with the columns `theta` (true ability), `mu` and
#' `sigma2` (conditional mean and variance of the final ability estimate), `bias`
#' (`mu - theta`), and `csem` (square root of `sigma2`).}
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [irtQ::run_mst()], [irtQ::panel_info()], [irtQ::find_cut()],
#'   [irtQ::shape_df()], [irtQ::est_score()]
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
#' \donttest{
#' ## ------------------------------------------------------------------------------
#' ## Evaluation of a 1-3-3 MST panel using simMST data.
#' ## This panel was assembled using module-level target TIFs and design
#' ## constraints similar to those used in Lim et al.'s (2021) simulation study
#' ## (it is not the identical dataset from that study).
#' ## Details:
#' ##    (a) Panel configuration: 1-3-3 MST panel
#' ##    (b) Test length: 24 items (each module contains 8 items across all stages)
#' ##    (c) IRT model: 3-parameter logistic model (3PLM)
#' ## ------------------------------------------------------------------------------
#' # Load the necessary library
#' library(dplyr)
#' library(tidyr)
#' library(ggplot2)
#'
#' # Import item bank metadata
#' x <- simMST$item_bank
#'
#' # Import module information
#' module <- simMST$module
#'
#' # Import routing map
#' route_map <- simMST$route_map
#'
#' # Import cut scores for routing to subsequent modules
#' cut_score <- simMST$cut_score
#'
#' # Import ability levels (theta) for evaluating measurement precision
#' theta <- simMST$theta
#'
#' # Evaluate MST panel using the reval_mst() function
#' eval <-
#'   reval_mst(x,
#'     D = 1.702, route_map = route_map, module = module,
#'     cut_score = cut_score, theta = theta, range.tcc = c(-7, 7)
#'   )
#'
#' # Review evaluation results
#' # The evaluation result table below includes conditional biases and
#' # standard errors of measurement (CSEMs) across ability levels
#' print(eval$eval.tb)
#'
#' # Generate plots for biases and CSEMs
#' p_eval <-
#'   eval$eval.tb %>%
#'   dplyr::select(theta, bias, csem) %>%
#'   tidyr::pivot_longer(
#'     cols = c(bias, csem),
#'     names_to = "criterion", values_to = "value"
#'   ) %>%
#'   ggplot2::ggplot(mapping = ggplot2::aes(x = theta, y = value)) +
#'   ggplot2::geom_point(mapping = ggplot2::aes(shape = criterion), size = 3) +
#'   ggplot2::geom_line(
#'     mapping = ggplot2::aes(
#'       color = criterion,
#'       linetype = criterion
#'     ),
#'     linewidth = 1.5
#'   ) +
#'   ggplot2::labs(x = expression(theta), y = NULL) +
#'   ggplot2::theme_bw() +
#'   ggplot2::theme(legend.key.width = unit(1.5, "cm"))
#' print(p_eval)
#' }
#'
#' @export
#' @import dplyr
reval_mst <- function(x,
                      D = 1,
                      route_map,
                      module,
                      cut_score,
                      theta = seq(-5, 5, 1),
                      intpol = TRUE,
                      range.tcc = c(-7, 7),
                      tol = 1e-4) {

  # Validate and normalize the item metadata (e.g., a missing guessing
  # parameter of 1PLM and 2PLM items is set to 0)
  x <- confirm_df(x)

  ## -----------------------------------------------------
  # Extract all the panel information from the route map
  ## -----------------------------------------------------
  # Run the panel_info() function
  panel_data <- panel_info(route_map)

  # Pathways allowed
  pathway <- panel_data$pathway

  # Number of modules for each stage
  n.mod <- panel_data$n.module

  # Total number of modules across all stages
  tn.mod <- sum(n.mod)

  # Stop when the module matrix does not have one column per module of the panel
  if (ncol(module) != tn.mod) {
    stop("'module' must have one column per module in 'route_map'.", call. = FALSE)
  }

  # Number of stages
  n.stg <- panel_data$n.stage

  # Stop when stage 1 has more than one module
  if (n.mod[1] != 1L) {
    stop("The first stage must have a single routing module.", call. = FALSE)
  }

  # Stop when cut_score does not have one vector per stage transition
  if (length(cut_score) != (n.stg - 1)) {
    stop(sprintf(
      "'cut_score' must be a list of length %d (one vector per stage transition).",
      n.stg - 1), call. = FALSE)
  }

  # Stop when a cut score vector does not separate the modules of the next stage
  for (s in seq_len(n.stg - 1)) {
    # treat a NULL element as an empty vector of cut scores
    if (is.null(cut_score[[s]])) cut_score[s] <- list(numeric(0))
    cut_s <- cut_score[[s]]
    if (!is.numeric(cut_s) || length(cut_s) != (n.mod[s + 1] - 1) ||
        !all(is.finite(cut_s)) || is.unsorted(cut_s, strictly = TRUE)) {
      # name the empty vector when the next stage has a single module
      if (n.mod[s + 1] == 1) {
        stop(sprintf(paste0(
          "'cut_score[[%d]]' must be empty (numeric(0) or NULL) because ",
          "stage %d has a single module."), s, s + 1), call. = FALSE)
      }
      stop(sprintf(paste0(
        "'cut_score[[%d]]' must contain %d finite value(s) in strictly ascending ",
        "order, one fewer than the number of modules in stage %d."),
        s, n.mod[s + 1] - 1, s + 1), call. = FALSE)
    }
  }

  ## -----------------------------------------------------
  # Create List of item metadata for modules and paths
  # across the stages
  ## -----------------------------------------------------
  # List containing item metadata of all modules
  meta_mod <-
    purrr::map(
      .x = 1:tn.mod,
      .f = ~ {
        x[module[, .x] == 1, ] %>%
          tibble::remove_rownames()
      }
    )
  names(meta_mod) <- paste0("m.", 1:tn.mod)

  # Check that all modules in the same stage have the same maximum sum score
  for (s in 1:n.stg) {
    # maximum sum scores of the modules in stage s
    max_sum_stg <-
      purrr::map_dbl(
        .x = meta_mod[panel_data$config[[s]]],
        .f = ~ sum(.x$cats - 1)
      )
    if (length(unique(max_sum_stg)) > 1L) {
      stop(
        sprintf(
          "All modules in stage %d must have the same maximum sum score (found: %s).",
          s, paste(max_sum_stg, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  # List containing item metadata for all possible (sub) pathways at each stage
  meta_path <-
    purrr::map(
      .x = 1:n.stg,
      .f = ~ {
        # Unique pathways up to the current stage
        path_uni <- unique(pathway[, 1:.x, drop = FALSE])

        # Item metadata for the unique pathways
        meta_uni <-
          purrr::map(
            .x = 1:nrow(path_uni),
            .f = ~ {
              meta_mod[path_uni[.x, ]] %>%
                dplyr::bind_rows()
            }
          )

        # Add the names and return the list
        names(meta_uni) <- paste0("path.", 1:length(meta_uni))
        return(meta_uni)
      }
    )
  names(meta_path) <- paste0("stage.", 1:n.stg)

  # Count the cumulative numbers of all possible (sub) pathways by each stage
  cum_path <- purrr::map_dbl(.x = meta_path,
                             .f = ~ {length(.x)})

  ## -----------------------------------------------------
  # Estimate the IRT theta scores using the Inverse-TCC method
  ## -----------------------------------------------------
  # Compute the IRT thetas corresponding to all observed sum scores using
  # the Inverse-TCC scoring method for all (sub) pathways
  eq_theta <-
    purrr::map(
      .x = meta_path,
      .f = ~ {
        # Implement the Inverse-TCC scoring for the current pathway
        tmp_meta <- .x
        tmp_eqtheta <-
          purrr::map(
            .x = tmp_meta,
            .f = ~ {
              inv_tcc_nr(
                x = .x, D = D,
                tol = tol, intpol = intpol,
                range.tcc = range.tcc
              )$est.theta
            }
          ) %>%
          do.call(what = "cbind")
      }
    )

  # Stop when a sum score has no inverse TCC estimate, which happens when
  # intpol = FALSE or range.tcc does not cover the estimates
  if (anyNA(unlist(eq_theta))) {
    stop(paste0(
      "Inverse TCC estimates are not available for some sum scores. ",
      "Use 'intpol = TRUE' and a 'range.tcc' wide enough to cover the estimates."),
      call. = FALSE)
  }

  ## -----------------------------------------------------
  # Compute the conditional distributions of modules at each theta
  ## -----------------------------------------------------
  # List of the conditional sum score distributions of each module,
  # one column per ability level
  cdist_by_mod <-
    purrr::map(
      .x = meta_mod,
      .f = ~ {
        lwrc(x = .x, theta = theta, D = D)
      }
    )

  # Rearrange the list by ability level: each element holds the
  # conditional distributions of all modules at one ability level
  cdist_by_th <-
    purrr::map(
      .x = 1:length(theta),
      .f = ~ {
        num <- .x
        tmp_tb <-
          purrr::map(
            .x = cdist_by_mod,
            .f = ~ {
              .x[, num]
            }
          ) %>%
          bind.fill(type = "cbind", fill = 0L)
        colnames(tmp_tb) <- paste0("m.", 1:tn.mod)
        tmp_tb
      }
    )
  names(cdist_by_th) <- theta

  ## -----------------------------------------------------
  # Compute the joint conditional distributions at each theta
  ## -----------------------------------------------------
  # List of possible observed scores across all modules
  score_mod <-
    purrr::map(
      .x = meta_mod[pathway[1, ]],
      .f = ~ {
        0:sum(.x$cats - 1)
      }
    )
  names(score_mod) <- paste0("stage.", 1:n.stg)

  # Cumulative maximum sum scores across the stages
  score_maxcum <- cumsum(sapply(X = score_mod, FUN = "max"))

  # Empty list to include the cut scores used to route a test taker to next module across all stages
  cut4route <- vector("list", n.stg - 1)
  names(cut4route) <- paste0("stage.", 2:n.stg)

  # Empty list to include path labels across all stages
  path_label <- vector("list", n.stg - 1)
  names(path_label) <- paste0("stage.", 2:n.stg)

  # Empty list to include probabilities of the joint distributions across the stages
  # Each element of the list will contain the joint distributions across all ability levels
  joint_dist <- vector("list", n.stg)
  names(joint_dist) <- paste0("stage.", 1:n.stg)

  # For the routing module in the first stage,
  # assign the conditional distributions of the routing modules across the ability levels
  joint_dist[[1]] <-
    purrr::map(
      .x = cdist_by_th,
      .f = ~ {
        .x[, panel_data$config[[1]], drop = FALSE] %>%
          unname()
      }
    )
  names(joint_dist[[1]]) <- round(theta, 10)

  # Calculate the probabilities of the conditional joint distributions
  # given the ability levels across all stages
  for (s in 2:n.stg) {
    # The numbers of all possible (sub) pathways at the current stage
    n_path <- cum_path[s - 1]

    # Empty list to contain the cut scores used to route a test taker to next module
    cut4route[[s - 1]] <- vector("list", n_path)

    # A data frame with the partial pathway up to the previous stage
    # (first column) and the module of the current stage (second column)
    path_pre <-
      data.frame(pathway[, 1:s]) %>%
      tidyr::unite(col = "path", 1:(s - 1), sep = "_")

    # Unique pathways by the current stage
    uni_path_pre <- unique(path_pre$path)

    # Empty list to contain the conditional joint distributions at the current stage
    jdist_temp <- vector("list", n_path)

    # Loop over all possible (sub) pathways
    for (i in 1:n_path) {
      # Find cut scores for routing a next module
      # Firstly, find all unique next modules that will be assigned at the current pathway
      nmod_next <-
        path_pre %>%
        dplyr::filter(.data$path %in% uni_path_pre[i]) %>%
        dplyr::pull(2) %>%
        unique()

      # Secondly, find the positions of the reachable modules among the modules
      # of the next stage in ascending order of module index
      idx_next <- match(nmod_next, panel_data$config[[s]])

      # Lastly, find the cut scores to be used to assign the next modules
      cut4route[[s - 1]][[i]] <- cut_score[[s - 1]][dplyr::lag(idx_next)[-1]]

      # Find the path labels to determine which module will be assigned to each
      # ability level based on the cut scores
      path_label[[s - 1]][[i]] <-
        give_path(
          score = unlist(eq_theta[[s - 1]][, i]),
          cut_sc = cut4route[[s - 1]][[i]]
        )$path

      # Cumulative maximum sum scores by the current stages
      max_score <- score_maxcum[s]

      # Number of all unique next modules
      n_route <- length(nmod_next)

      # All possible observed scores by the current stage
      score_cur <- 0:score_maxcum[s - 1]

      # Possible observed score for the next module
      score_next <- score_mod[[s]]

      # The conditional joint distributions given the ability levels by the current stage
      jdist_cur <-
        purrr::map(
          .x = joint_dist[[s - 1]],
          .f = ~ {
            .x[, i, drop = FALSE]
          }
        )

      # The conditional distribution of the next modules given the ability levels
      cdist_next <-
        purrr::map(
          .x = cdist_by_th,
          .f = ~ {
            .x[, nmod_next, drop = FALSE]
          }
        )

      # Column bind of the above two distributions given the ability levels
      cbind_dist <-
        purrr::map(
          .x = 1:length(theta),
          .f = ~ {
            bind.fill(
              List = list(jdist_cur[[.x]], cdist_next[[.x]]),
              type = "cbind",
              fill = 0L
            )
          }
        )

      # Compute the joint conditional distributions
      jdist_temp[[i]] <-
        jdist(
          cbind_dist = cbind_dist, path_lb = path_label[[s - 1]][[i]],
          n_route = n_route, max_score = max_score, score_cur = score_cur,
          score_next = score_next, theta = theta
        )
    }

    # Add to the list
    joint_dist[[s]] <-
      purrr::map(
        .x = 1:length(theta),
        .f = ~ {
          do.call("cbind", lapply(jdist_temp, "[[", .x))
        }
      )
    names(joint_dist[[s]]) <- round(theta, 10)
  }

  ## -----------------------------------------------------
  # Evaluate the measurement precision of the MST panel
  ## -----------------------------------------------------
  # Calculate the means and variances of the conditional joint distributions
  # across all ability levels
  cmoment <-
    lapply(
      X = joint_dist[[s]],
      FUN = cal_moment,
      node = eq_theta[[s]]
    ) %>%
    do.call(what = "rbind")

  # Create a table containing the moments and
  # conditional biases and standard errors of measurement (SEMs)
  eval_tb <-
    data.frame(theta = theta, cmoment) %>%
    dplyr::mutate(
      bias = .data$mu - .data$theta,
      csem = sqrt(.data$sigma2)
    ) %>%
    tibble::remove_rownames()


  ## -----------------------------------------------------
  # Return the results
  ## -----------------------------------------------------
  rst <- list(
    panel.info = panel_data,
    item.by.mod = meta_mod,
    item.by.path = meta_path,
    eq.theta = eq_theta,
    cdist.by.mod = cdist_by_th,
    jdist.by.path = joint_dist,
    eval.tb = eval_tb
  )
  rst
}


# This function computes the conditional joint distributions of the current pathways
# given the ability levels
jdist <- function(cbind_dist, path_lb, n_route, max_score, score_cur,
                  score_next, theta) {
  # Empty list across all ability levels
  cj_dist_all <- vector("list", length(theta))

  # Loop over the ability levels
  for (j in 1:length(theta)) {
    # Matrix to contain the conditional joint distribution at a given ability level
    cj_dist <- array(0, c(max_score + 1, n_route))
    rownames(cj_dist) <- 0:max_score

    # Loop over the observed scores of the (sub) pathway up to the current stage
    for (i in 1:length(score_cur)) {
      # Compute the conditional joint distribution
      cj_dist[i:(i + max(score_next)), path_lb[i]] <-
        cj_dist[i:(i + max(score_next)), path_lb[i]] +
        (cbind_dist[[j]][1:length(score_next), (path_lb[i] + 1)] * cbind_dist[[j]][i, 1])
    }

    # Add to the list
    cj_dist_all[[j]] <- cj_dist
  }

  # Return the results
  cj_dist_all
}

# This function determines the next path that each test taker with
# a given ability level needs to take
give_path <- function(score, cut_sc) {
  # Number of theta scores
  nscore <- length(score)

  # Number of categories (or labels) to which each theta will be given
  ncats <- length(cut_sc) + 1

  # A vector of labels to be given to each test taker
  labels <- 1:ncats

  # Create pathway vectors; intervals are closed on the left, [c_k, c_(k+1)),
  # so a score equal to a cut score goes to the higher module, and the last
  # interval includes Inf
  path <-
    cut(x = score, breaks = c(-Inf, cut_sc, Inf), labels = labels,
        right = FALSE, include.lowest = TRUE) %>%
    as.numeric()

  # Return the results
  rst <- list(path = path, score = score, cut.sc = cut_sc)
  rst
}
