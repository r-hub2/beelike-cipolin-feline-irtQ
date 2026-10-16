#' Extract Panel Structure from an MST Route Map
#'
#' @description
#' Given a binary transition matrix that defines the Multistage-Adaptive Test
#' (MST) structure, this function extracts key panel information: the module
#' configuration by stage, all valid test pathways, the number of modules per
#' stage, and the total number of stages.
#'
#' @param route_map A binary square matrix defining the MST transition structure,
#'   where a value of 1 at row \emph{i}, column \emph{j} indicates that a test
#'   taker can be routed from module \emph{i} to module \emph{j}. This is
#'   equivalent to the \code{transMatrix} argument in the \code{randomMST()}
#'   function of the \pkg{mstR} package (Magis et al., 2017). See
#'   \code{\link{reval_mst}} for details on how to construct this matrix.
#'
#' @return A named list with four elements:
#' \describe{
#'   \item{\code{config}}{A named list of length equal to the number of stages.
#'     Each element is an integer vector of module indices belonging to that
#'     stage. Names are \code{"stage.1"}, \code{"stage.2"}, etc.}
#'   \item{\code{pathway}}{An integer matrix of all valid test pathways, with
#'     one row per pathway and one column per stage. Row names are
#'     \code{"path.1"}, \code{"path.2"}, etc., and column names are
#'     \code{"stage.1"}, \code{"stage.2"}, etc. Each entry is the index of the
#'     module administered at that stage.}
#'   \item{\code{n.module}}{A named integer vector giving the number of modules
#'     in each stage.}
#'   \item{\code{n.stage}}{An integer scalar giving the total number of stages.}
#' }
#'
#' @details
#' The function identifies stage-1 modules as those with no incoming transitions
#' (column sum of \code{route_map} equals zero). It then traverses the route map
#' stage by stage, collecting modules reachable from the previous stage. All
#' logically possible pathways are enumerated via \code{expand.grid()}, and
#' pathways that violate the transition constraints in \code{route_map} are
#' removed. The remaining pathways are sorted and returned as a matrix.
#'
#' The route map must be a square matrix of 0s and 1s, and the modules do not
#' need to be numbered in the order of the stages. The function stops with an
#' error when the route map defines a single stage, contains a cycle, has a
#' module with no transition to or from another module, has a module that
#' cannot be reached from a stage-1 module or that belongs to more than one
#' stage, or has a module before the last stage without a transition to a
#' module of the next stage. A logical route map is read as a 0/1 matrix.
#'
#' @references
#'   Magis, D., Yan, D., & von Davier, A. A. (2017). *Computerized adaptive and
#'   multistage testing with R: Using packages catR and mstR*. Springer.
#'   \doi{10.1007/978-3-319-69218-0}.
#'
#' @seealso \code{\link{reval_mst}}, \code{\link{run_mst}}
#'
#' @examples
#' # Use the 1-3-3 simMST panel route map
#' route_map <- simMST$route_map
#' pinfo <- panel_info(route_map)
#'
#' # Stage configuration: which modules belong to which stage
#' pinfo$config
#'
#' # All valid pathways (rows = pathways, columns = stages)
#' pinfo$pathway
#'
#' # Number of modules per stage
#' pinfo$n.module
#'
#' # Total number of stages
#' pinfo$n.stage
#'
#' @export
panel_info <- function(route_map) {

  # Transform the format of the route_map to a numeric matrix
  route_map <- as.matrix(route_map)

  # Treat a logical route map as a 0/1 matrix
  if (is.logical(route_map)) storage.mode(route_map) <- "integer"

  # Stop when route_map is not a square binary matrix
  if (!is.numeric(route_map) || nrow(route_map) != ncol(route_map) ||
      anyNA(route_map) || !all(route_map %in% c(0, 1))) {
    stop("'route_map' must be a square binary (0/1) matrix.", call. = FALSE)
  }

  # Count the total number of modules (columns = total modules)
  end_col <- ncol(route_map)

  # Find the routing module (the module(s) in the first stage):
  # A first-stage module has no incoming transitions -> column sum equals 0
  mod_stg1 <- unname(which(colSums(route_map) == 0))

  # Stop when a module has no transition to or from another module
  mod_iso <- which(rowSums(route_map) == 0 & colSums(route_map) == 0)
  if (end_col > 1L && length(mod_iso) > 0L) {
    stop(sprintf("Module(s) %s in 'route_map' have no transition to or from another module.",
                 paste(mod_iso, collapse = ", ")), call. = FALSE)
  }

  # Stop when every module has an incoming transition, which happens when the map has a cycle
  if (length(mod_stg1) == 0L) {
    stop("'route_map' must not contain a cycle.", call. = FALSE)
  }

  # Traverse the route map stage by stage to build the config list
  config <- list()
  config[[1]] <- mod_stg1
  col_bef <- mod_stg1
  repeat {
    # Modules in the next stage: columns with at least one "1" in the current
    # stage's rows
    col_aft <- unname(which(colSums(route_map[col_bef, , drop = FALSE] == 1) > 0))
    # Stop when the current stage has no outgoing transitions (final stage)
    if (length(col_aft) == 0L) break
    # Stop when the stages never end, which happens when the map has a cycle
    if (length(config) >= end_col) {
      stop("'route_map' must not contain a cycle.", call. = FALSE)
    }
    config[[length(config) + 1L]] <- col_aft
    col_bef <- col_aft
  }

  # Stop when the route map has no transition, which gives a single stage
  if (length(config) < 2L) {
    stop("'route_map' must define at least two stages.", call. = FALSE)
  }

  # Stop when a module is unreachable or belongs to more than one stage
  if (!identical(sort(unlist(config)), seq_len(end_col))) {
    stop("Each module in 'route_map' must belong to exactly one stage.",
         call. = FALSE)
  }

  # Stop when a module before the final stage has no outgoing transition
  mod_nonfinal <- unlist(config[-length(config)])
  if (any(rowSums(route_map[mod_nonfinal, , drop = FALSE]) == 0)) {
    stop("Every module before the final stage must route to a module of the next stage.",
         call. = FALSE)
  }

  # Name stages and strip names from module vectors
  names(config) <- paste0("stage.", 1:length(config))
  config <- lapply(X = config, FUN = "unname")

  # Count stages and modules per stage
  n.mod <- sapply(X = config, FUN = "length")
  n.stg <- length(n.mod)

  # Enumerate all logically possible pathways (before applying constraints)
  path_all <- expand.grid(config)

  # Remove pathways that violate transition constraints in route_map
  remain <- path_all
  for (s in 1:(n.stg - 1)) {
    for (m in 1:n.mod[s]) {
      tmp.1 <- config[[s]][m]                        # current stage module
      tmp.2 <- which(route_map[tmp.1, ] == 1)        # allowed next modules
      # Delete rows where stage-s module is tmp.1 but next stage module is NOT
      # in the allowed set
      delete <- (remain[, s] == tmp.1 & !(remain[, s + 1] %in% tmp.2))
      remain <- remain[!delete, ]
    }
  }

  # Sort pathways and label rows
  remain <-
    dplyr::arrange(remain, dplyr::pick(dplyr::everything())) %>%
    as.matrix()
  rownames(remain) <- paste0("path.", 1:nrow(remain))

  # Return the panel structure
  rst <- list(config = config, pathway = remain, n.module = n.mod, n.stage = n.stg)
  rst
}
