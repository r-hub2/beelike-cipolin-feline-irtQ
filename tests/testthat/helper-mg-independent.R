# Independent formulas and data for the multiple-group tests (built on helper-em-independent.R)

# per-examinee marginal log-likelihood of one group given its quadrature weights
indep_marg_i <- function(x, data, theta, w, D) {
  indep_lse(sweep(indep_llmat(x, data, theta, D), 2, log(w), "+"))
}

# multiple-group marginal log-likelihood by group: x is the overall item metadata, data a list of
# group response matrices whose column names are item IDs, and wts a list of (theta, weight) data frames
indep_mg_loglik <- function(x, data, wts, D) {
  vapply(seq_along(data), function(k) {
    xk <- x[match(colnames(data[[k]]), x$id), ]
    sum(indep_marg_i(xk, data[[k]], wts[[k]][, 1], wts[[k]][, 2], D))
  }, numeric(1))
}

# normalized normal weights on a grid
indep_normw <- function(theta, mu, s2) {
  w <- dnorm(theta, mu, sqrt(s2))
  w / sum(w)
}

# data of three groups whose forms share the common items C1 to C6 (Group 1 and 2 share C1 to C4,
# Group 2 and 3 share C5 and C6) and have unique items
make_mg_data <- function(seed = 11, n = c(300, 300, 300)) {
  set.seed(seed)
  ids <- c(paste0("C", 1:6), paste0("A", 1:4), paste0("B", 1:4), paste0("E", 1:4))
  mods <- c("2PLM", "3PLM", "2PLM", "GRM", "2PLM", "GPCM", "3PLM", "2PLM", "2PLM", "GRM", "2PLM", "3PLM", "2PLM",
            "2PLM", "2PLM", "GPCM", "3PLM", "2PLM")
  cats <- ifelse(mods %in% c("GRM", "GPCM"), 3, 2)
  a <- rlnorm(length(ids), 0, .25)
  b <- rnorm(length(ids))
  pd <- lapply(seq_along(ids), function(i) {
    if (cats[i] == 2) c(a[i], b[i], if (mods[i] == "3PLM") .15 else NA) else c(a[i], sort(rnorm(2)))
  })
  xall <- data.frame(id = ids, cats = cats, model = mods, par.1 = sapply(pd, `[`, 1), par.2 = sapply(pd, `[`, 2),
                     par.3 = sapply(pd, `[`, 3), stringsAsFactors = FALSE)
  forms <- list(c(paste0("C", 1:4), paste0("A", 1:4)), c(paste0("C", 1:6), paste0("B", 1:4)), c(paste0("C", 5:6), paste0("E", 1:4)))
  xs <- lapply(forms, function(f) {
    r <- xall[match(f, xall$id), ]
    rownames(r) <- NULL
    r
  })
  th <- list(rnorm(n[1], 0, 1), rnorm(n[2], 0.5, 1.2), rnorm(n[3], -0.3, 0.8))
  dat <- lapply(1:3, function(k) {
    d <- simdat(x = xs[[k]], theta = th[[k]], D = 1)
    colnames(d) <- xs[[k]]$id
    d
  })
  list(x = xs, data = dat, xall = xall)
}
