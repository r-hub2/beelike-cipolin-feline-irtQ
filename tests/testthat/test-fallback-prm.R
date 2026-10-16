# the search of the first EM cycles keeps the previous values of a polytomous item when it fails

test_that("a failed search in the first EM cycles keeps the previous values of a GRM item", {
  quadpt <- seq(-4, 4, length.out = 21)
  f_i <- rep(20, 21)
  r_i <- round(f_i * plogis(quadpt))
  r3 <- cbind(f_i - r_i, round(r_i / 2), r_i - round(r_i / 2))
  # a gradient that cannot be evaluated makes every gradient-based search fail
  fnp <- irtQ:::make_prm_optim_fns
  local_mocked_bindings(make_prm_optim_fns = function(...) {
    out <- fnp(...)
    out$gradient <- function(p) rep(NaN, length(p))
    out
  })
  for (it in 1:2) {
    est <- irtQ:::estimation2(
      f_i = f_i, r_i = r3, s_i = NULL, quadpt = quadpt, mod = "GRM", D = 1, n.quad = 21,
      control = list(eval.max = 50, iter.max = 20), startval = c(1, -0.3, 0.4),
      lower = c(0.01, -6, -6), upper = c(5, 6, 6), iter = it
    )
    expect_equal(est$pars, c(1, -0.3, 0.4))
    expect_equal(est$convergence, 1L)
  }
})
