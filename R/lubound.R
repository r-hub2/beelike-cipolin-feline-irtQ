# This function sets the lower and upper bounds of the item parameters during the estimation.
# The slopes of dichotomous items are unbounded, the slopes of GRM and GPCM items have a lower
# bound of 1e-8, and the guessing parameters lie in [1e-8, 1 - 1e-8].
lubound <- function(model, cats, n.1PLM, idx4est = NULL, fix.a.1pl, fix.g, fix.a.gpcm) {
  bound.par <- list()

  if (!is.null(idx4est$drm.slc)) {
    lower <- c(-Inf, rep(-Inf, n.1PLM))
    upper <- c(Inf, rep(Inf, n.1PLM))
    bound.par$drm.slc <- list(lower = lower, upper = upper)
  } else {
    bound.par$drm.slc <- NULL
  }

  if (!is.null(idx4est$drm.else)) {
    bound.par$drm.else <- list()
    for (i in 1:length(idx4est$drm.else)) {
      if (model[idx4est$drm.else[i]] == "1PLM") {
        lower <- -Inf
        upper <- Inf
        bound.par$drm.else[[i]] <- list(lower = lower, upper = upper)
      } else if (model[idx4est$drm.else[i]] == "2PLM") {
        lower <- c(-Inf, -Inf)
        upper <- c(Inf, Inf)
        bound.par$drm.else[[i]] <- list(lower = lower, upper = upper)
      } else if (model[idx4est$drm.else[i]] == "3PLM" & fix.g) {
        lower <- c(-Inf, -Inf)
        upper <- c(Inf, Inf)
        bound.par$drm.else[[i]] <- list(lower = lower, upper = upper)
      } else if (model[idx4est$drm.else[i]] == "3PLM" & !fix.g) {
        lower <- c(-Inf, -Inf, 1e-8)
        upper <- c(Inf, Inf, 1 - (1e-8))
        bound.par$drm.else[[i]] <- list(lower = lower, upper = upper)
      }
    }
  } else {
    bound.par$drm.else <- NULL
  }

  if (!is.null(idx4est$prm)) {
    bound.par$prm <- list()
    for (i in 1:length(idx4est$prm)) {
      if (model[idx4est$prm[i]] == "GRM") {
        lower <- c(1e-8, rep(-Inf, (cats[idx4est$prm[i]] - 1)))
        upper <- c(Inf, rep(Inf, (cats[idx4est$prm[i]] - 1)))
        bound.par$prm[[i]] <- list(lower = lower, upper = upper)
      } else if (model[idx4est$prm[i]] == "GPCM" & fix.a.gpcm) {
        lower <- rep(-Inf, (cats[idx4est$prm[i]] - 1))
        upper <- rep(Inf, (cats[idx4est$prm[i]] - 1))
        bound.par$prm[[i]] <- list(lower = lower, upper = upper)
      } else if (model[idx4est$prm[i]] == "GPCM" & !fix.a.gpcm) {
        lower <- c(1e-8, rep(-Inf, (cats[idx4est$prm[i]] - 1)))
        upper <- c(Inf, rep(Inf, (cats[idx4est$prm[i]] - 1)))
        bound.par$prm[[i]] <- list(lower = lower, upper = upper)
      }
    }
  } else {
    bound.par$prm <- NULL
  }

  # return the results
  bound.par
}
