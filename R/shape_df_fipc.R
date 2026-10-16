#' Combine Fixed and New Item Metadata for Fixed Item Parameter Calibration
#' (FIPC)
#'
#' This function merges existing fixed-item metadata with automatically
#' generated metadata for new items, producing a single data frame ordered by
#' specified test positions, to facilitate fixed item parameter calibration
#' using [irtQ::est_irt()].
#'
#' @param x A data.frame of metadata for items whose parameters remain fixed
#'   (e.g., output from [irtQ::shape_df()]).
#' @param fix.loc An integer vector specifying the row positions in the final
#'   output where the fixed items are placed, in the row order of `x`. This
#'   argument must be specified.
#' @param item.id A character vector of IDs for the new items whose parameters
#'   are estimated. If `NULL`, default IDs (`"V1"`, `"V2"`, ...) are assigned.
#' @param cats An integer vector giving the number of score categories of each
#'   new item, in the order of `item.id`. A single value is recycled across all
#'   new items.
#' @param model A character vector of IRT model names for the new items:
#'   `"1PLM"`, `"2PLM"`, `"3PLM"`, or `"DRM"` for dichotomous items, and
#'   `"GRM"` or `"GPCM"` for polytomous items. A single value is recycled
#'   across all new items.
#'
#' @details First, prepare a metadata frame `x` that contains only the fixed
#'   items, created by [irtQ::shape_df()] or imported from external software
#'   (e.g., with [irtQ::bring.flexmirt()]). It must include the columns `id`,
#'   `cats`, `model`, and the parameter columns (`par.1`, `par.2`, ...). Then
#'   use `fix.loc` to give the row positions of these items in the final test
#'   form. The length of `fix.loc` must equal the number of rows in `x`, and
#'   the i-th element of `fix.loc` is the position of the i-th row of `x`.
#'
#'   Next, provide `item.id`, `cats`, and `model` for the new items, whose
#'   number equals the total form length minus the length of `fix.loc`. A
#'   single value of `cats` or `model` is used for all new items.
#'   If `item.id` is `NULL`, default IDs (`"V1"`, `"V2"`, ...) are assigned.
#'
#'
#' @return A data frame of item metadata for all items (fixed and new), ordered
#'   by test position, with columns `id`, `cats`, `model`, and `par.1`,
#'   `par.2`, .... New items receive default parameters as in [irtQ::shape_df()]
#'   with `default.par = TRUE`. Fixed `"DRM"` items are relabeled `"3PLM"`.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [irtQ::shape_df()]
#'
#' @examples
#' ## Import the flexMIRT parameter output file
#' prm_file <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
#' x_fixed <- bring.flexmirt(file = prm_file, "par")$Group1$full_df
#'
#' ## Define positions of fixed items in the test form
#' fixed_pos <- c(1:40, 43:57)
#'
#' ## Specify IDs, models, and category counts for new items
#' new_ids <- paste0("NI", 1:6)
#' new_models <- c("3PLM", "1PLM", "2PLM", "GRM", "GRM", "GPCM")
#' new_cats <- c(2, 2, 2, 4, 5, 6)
#'
#' ## Generate combined metadata for FIPC
#' shape_df_fipc(x = x_fixed, fix.loc = fixed_pos, item.id = new_ids,
#'   cats = new_cats, model = new_models)
#'
#' @export
shape_df_fipc <- function(x, fix.loc = NULL, item.id = NULL, cats, model) {

  # Validate and standardize the fixed-item metadata
  x_fix <- confirm_df(x)

  # Number of new items: the length of item.id when given, else the longer of cats and model
  n_new <- if (!is.null(item.id)) length(item.id) else max(length(cats), length(model))

  # Repeat a single value of cats or model for every new item
  if (length(cats) == 1L) cats <- rep(cats, n_new)
  if (length(model) == 1L) model <- rep(model, n_new)

  # Stop when cats or model do not match the number of new items
  if (length(cats) != n_new || length(model) != n_new) {
    stop("The lengths of `cats` and `model` must be 1 or equal to the number of new items.", call. = FALSE)
  }

  # Generate default metadata for the new items
  x_new <- shape_df(item.id = item.id, cats = cats, model = model, default.par = TRUE)

  # Merge fixed and new metadata into a single data frame
  x_all <- dplyr::bind_rows(x_fix, x_new)

  # Determine the total number of items
  nitem <- nrow(x_all)

  # Stop unless fix.loc gives one distinct valid position for every fixed item
  if (is.null(fix.loc) || length(fix.loc) != nrow(x_fix) ||
      anyDuplicated(fix.loc) > 0L || !all(fix.loc %in% seq_len(nitem))) {
    stop(
      "`fix.loc` must contain one distinct position in 1:", nitem,
      " for each row of `x`.",
      call. = FALSE
    )
  }

  # Identify row positions reserved for new items
  nfix.loc <- setdiff(seq_len(nitem), fix.loc)

  # Extract the original fixed-item rows
  x_fix2 <- x_all[seq_len(nrow(x_fix)), ]

  # Extract the new-item rows aligned to every column of the merged frame
  x_new2 <- x_all[nrow(x_fix) + seq_len(nrow(x_new)), ]

  # Place fixed items and new items at their specified positions
  x_all[fix.loc, ] <- x_fix2
  x_all[nfix.loc, ] <- x_new2

  # Return the combined and ordered metadata
  x_all

}

