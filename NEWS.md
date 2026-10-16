
<!-- NEWS.md is generated from NEWS.Rmd. Please edit that file -->

# irtQ 1.4.0

## Bug Fixes

- `ctt()` returned a standardized alpha (`alpha_std`) that was too
  large, and could exceed 1, when a test included a constant item. The
  constant item is now kept in the item count as in the raw alpha, and
  the value is 0 when only one item varies, as is the raw alpha. An
  undefined standardized alpha or mean discrimination is `NA` instead of
  `NaN`.

- `ctt()`, `ctt_distr()`, and `freq_score()` silently left out of their
  counts a score within rounding error of a whole number, such as
  `(0.1 + 0.2) * 10`. Such a score is now counted as that whole number.

- `freq_score()` tabulated a factor by its level codes (1, 2, 3, ...)
  instead of its labels. A factor is now converted through its labels.

- `ctt_distr()` listed an option twice in selected-response mode when
  the responses wrote it in different ways, such as `"1"` and `"01"`,
  which inflated the frequencies and percentages. The option is now
  listed once.

- `score_resp()` and `ctt_distr()` did not treat a response with
  surrounding spaces, such as `" 9"`, as missing when `missing = "9"`.
  Responses are now compared with `missing` after trimming spaces.

- `score_resp()` and `ctt_distr()` silently scored every response to an
  item as wrong when its key value contained a comma, such as `"1,5"`.
  They now stop with an error, because every item must have exactly one
  correct option.

- `ctt()` did not flag an item whose discrimination is undefined, such
  as an item with a constant score, so the item could pass the
  `crit.dis` check. Such an item is now flagged with
  `"discrimination undefined"`.

- `shape_df_fipc()` misplaced the parameter columns of the new items
  when the fixed items included a polytomous item and every new item was
  dichotomous, which changed a threshold of a fixed polytomous item. The
  new items now follow the columns of the fixed items. FIPC results
  built from such metadata should be rerun.

- The GPCM probabilities of `prm()`, `traceline()`, and `info()` were
  distorted when an exponent exceeded 700, and the values at the other
  thetas and items of the same call changed with them. Each row is now
  rescaled by its own maximum. Results in the usual range change by less
  than 1e-12, and `prm()` no longer returns a GPCM probability of
  exactly 0.

- `info()` returned `NaN` instead of 0 for a dichotomous item at a theta
  far from the item location.

- `traceline()` returned a vector instead of a one-row matrix in
  `prob.cats` for a polytomous item when a single theta was given, which
  made `plot()` fail.

- `drm()` checked each probability against the guessing parameter of
  another item, and the expected item score of `traceline()` included
  categories that the item does not have. Probabilities and expected
  scores could differ by about 1e-10, and item parameter estimates by
  about 1e-9.

- `traceline()` and `info()` stopped with an error for a two-category
  GRM or GPCM item, such as one made by `bring.mirt()`, and for a 3PLM
  or DRM item with an `NA` guessing parameter, for which `drm()`
  returned `NA`. An `NA` guessing parameter is now treated as 0, so a
  two-category GRM or GPCM item is computed as a 2PLM item, also when
  the test has polytomous items of the same model.

- `shape_df()` used a single `g` value for the first dichotomous item
  only and set the guessing parameters of the other items to 0. A single
  `g` value is now used for every dichotomous item, as in `simdat()`,
  and a `g` vector of any other wrong length stops with an error.

- `est_score()` with `ncore > 1` failed when a chunk of examinees had a
  single row, for example with 3 examinees and `ncore = 2`, and it did
  not close the cluster when a worker failed. Each chunk is now kept as
  a matrix, also when it has one row, and the cluster is always closed.

- `est_score()` with `method = "ML"`, `"WL"`, `"MLF"`, or `"MAP"` could
  stop at a point that was not a solution when the iterations alternated
  between two values, which happened mostly for low-scoring examinees on
  tests with 3PLM items. The estimate then depended on `max.iter`. When
  the iterations do not converge, the root of the score function is now
  found within the range of the visited values.

- `est_score()` with `method = "INV.TCC"` returned wrong standard errors
  when `intpol = FALSE` or when the interpolation was skipped with a
  warning (a `range.tcc` that is too narrow), because the estimates were
  paired with the wrong sum scores of the score distribution. The
  standard deviation is now taken over the sum scores that have an
  estimate. Results with `intpol = TRUE` and no warning are unchanged.

- `est_score()` with `method = "MLF"` gave an examinee with all
  responses missing an estimate near 0 and a standard error of 316.2
  without a warning. The examinee now gets `NA` and the warning, as with
  the other methods.

- `est_score()` with `method = "EAP"` returned `NaN`, with a warning
  about missing responses, for tests of about 1,000 items or more,
  because the product of the item probabilities underflowed. The
  posterior is now computed from log-likelihoods, and results for usual
  tests change by less than 1e-14.

- `llike_score()` with `method = "MLF"` stopped with an error under the
  default `fence.b = NULL`. The fences are now placed at -5 and 5, the
  default `range` of `est_score()`. The function also merged examinees
  whose rows had the same row name, and an examinee without any
  response, who stopped the function (`"ML"`, `"MAP"`) or got values
  from the fence items alone (`"MLF"`), now gets `NA`.

- `lwrc()` stopped with an error for a single `theta` when the test had
  polytomous items.

- `sx2_fit()` counted some summed score groups twice in the tables of
  polytomous items when the test was so short that the groups merged at
  the two ends overlapped (for example, two items with five categories).
  The observed and expected frequencies and the degrees of freedom were
  too large. Each summed score group is now counted once.

- `sx2_fit()` and `irtfit()` reported a p-value of 0 (or `NaN` with a
  warning) for an item without degrees of freedom left after collapsing,
  which read as a strong misfit. The critical value and p-value of such
  an item are now `NA`, with a warning that names the item.

- `sx2_fit()` named the columns of the observed frequency and proportion
  tables (`obs_freq` and `obs_prop`) with an `exp_freq.` prefix
  (`exp_freq.score.0`, ...) when the test had a single polytomous item.
  The columns are now named `score.0`, `score.1`, ... as for tests with
  several polytomous items. The values are unchanged.

- `plot.irtfit()` ignored `xlab.text` when `type = "both"` and always
  used theta as the x-axis title. The given title is now used, as for
  the other types.

- `plot.irtfit()` drew the standardized residuals as blue crosses,
  instead of red circles, when every plotted residual exceeded `overSR`.
  Each kind of point now keeps its color and shape.

- `irtfit()` stopped with an error for a logical response matrix, and
  `irtfit()` and `sx2_fit()` stopped with an error when the groups of an
  item were merged into a single row (for example, with a large
  `min.collapse`). These cases now work; an item left without degrees of
  freedom gets `NA` for the critical value and p-value of the chi-square
  statistic, with a warning.

- `sx2_fit()` for an `est_irt` object computed the expected frequencies
  with a standard normal distribution even when a different latent
  distribution was used or estimated in the calibration (for example,
  with `EmpHist = TRUE` or FIPC). When `weights`, `norm.prior`, and
  `nquad` are all omitted, the latent distribution stored in the object
  is now used. The statistics change slightly for a calibration with a
  standard normal distribution and can change enough to alter the
  decision at `alpha` when the distribution was estimated. Set
  `norm.prior = c(0, 1)` to reproduce the earlier results.

- `run_mst()` stopped with an obscure error for a panel whose routing
  module is not module 1 when `ini_mod` was `NULL`, because the starting
  module was drawn with `sample()` from a single module number. The
  examinees now start in the routing module.

- `reval_mst()` took the score distribution of the first stage from
  module 1 even when the routing module had another number, which gave
  wrong results without an error. The routing module is now taken from
  `route_map`, and a first stage with more than one module stops with an
  error. `run_mst()` still runs such panels.

- `cac_lee()` and `cac_rud()` used weights that did not sum to 1 as
  given, so the marginal indices could exceed 1. Such weights are now
  rescaled to sum to 1, with a warning. Weights that are missing,
  infinite, or negative, or that have a non-positive sum, now stop with
  an error.

- `plot.find_cut()` with `theta_range` still stretched the x-axis to the
  crossing lines outside the range. The x-axis now covers the requested
  range.

- `ripd()` stopped with an obscure error when an examinee had no ability
  estimate: a missing value in `score`, an examinee without any
  response, a focal group examinee left out by `min.resp` (so `min.resp`
  could not be used), or an examinee whose responses were all removed
  during purification. Such examinees are now excluded, and the function
  stops when no examinee has an ability estimate. A warning reports the
  examinees excluded because of `score` or `min.resp`, and one warning
  at the end of the purification reports the examinees excluded during
  it. `min.resp` now also applies to the first analysis when `score` is
  given, and an examinee with a missing value in `score` stays out of
  the purification. The purification stops with a warning and
  `complete = FALSE` when `min.resp`, which counts the items that
  remain, leaves a group without any ability estimate; `score` and the
  number of examinees reported as excluded then belong to the last
  iteration whose statistics are reported.

- `ripd()` and `pcd2()` stopped with an obscure error when purification
  left a single item. They now work. When every item is flagged, the
  purification stops with a warning and returns all items as flagged
  with `complete = FALSE`.

- `ripd()` and `pcd2()` silently accepted responses outside the score
  categories (for example, an unrecorded missing code such as -9, which
  gave chi-square values in the hundreds of thousands), a `group` of the
  wrong length (`ripd()`), and an `item.skip` that was not a set of item
  positions (an ID given to `pcd2()` added a spurious row and skipped
  nothing). They now stop with a clear message, a logical `item.skip` is
  read as item positions, also during purification, and character
  responses are read as numbers.

- `pcd2()` warned about excluded items even when `min.resp` excluded
  none, and worded the warning poorly. It now warns only when items are
  excluded.

- `ripd()` compared the p-values rounded to four decimals with `alpha`,
  so an item with a p-value just above `alpha` (for example, 0.05004 for
  `alpha = 0.05`) was flagged. The decision now uses the unrounded
  p-values, and the tables still show four decimals.

- `ripd()` computed `RIPD_RS` from arbitrary values when the covariance
  matrix of `RIPD_R` and `RIPD_S` was singular, which happens when all
  examinees who responded to an item have the same probability of a
  correct response (a tiny constant was added to every entry of the
  matrix). The statistic is now computed with a generalized inverse and
  the rank of the matrix as degrees of freedom, and a warning names the
  items. Because the degrees of freedom can then differ between items,
  the purification removes the flagged item with the smallest p-value of
  the `purify.by` statistic, compared on the log scale, instead of the
  item with the largest statistic. Results for regular data are
  unchanged.

- `rdif()` and `crdif()` stopped with an obscure error when an examinee
  had no ability estimate: a missing value in `score`, an examinee
  without any response, an examinee left out by `min.resp` (so
  `min.resp` could not exclude any examinee), or an examinee whose
  responses were all removed during purification. Such examinees are now
  excluded, and the functions stop when no examinee has an ability
  estimate. A warning reports the examinees excluded because of `score`
  or `min.resp`, and one warning at the end of the purification reports
  the examinees excluded during the purification.

- `rdif()` and `crdif()` stopped with an obscure error when purification
  left a single item, and `rdif()` also stopped with "non-conformable
  arguments" when purification removed the item with the most score
  categories from a test with polytomous items. They now work. When
  every item is flagged, the purification stops with a warning and
  returns all items as flagged with `complete = FALSE`. The purification
  also stops with a warning and `complete = FALSE` when `min.resp`,
  which counts the items that remain, leaves a group without any ability
  estimate; `score` and the number of examinees reported as excluded
  then belong to the last iteration whose statistics are reported.

- `rdif()` and `crdif()` silently accepted responses outside the score
  categories (for example, an unrecorded missing code such as -9, which
  flagged every item), a `group` of the wrong length, a `focal.name`
  that is not in `group`, an `item.skip` that was not a set of item
  positions (a logical `item.skip` skipped only the first item during
  purification), an `alpha` outside (0, 1), a `max.iter` that was not a
  positive whole number when `purify = TRUE`, and a `group` with more
  than two values (the values other than `focal.name` were pooled into
  the reference group). They now stop with a clear message; use
  `grdif()` to compare more than two groups. A logical `item.skip` is
  read as item positions, also during purification, and character
  responses are read as numbers.

- `rdif()` compared the p-values rounded to four decimals with `alpha`,
  so an item with a p-value just above `alpha` (for example, 0.05004 for
  `alpha = 0.05`) was flagged, and an item with a p-value below 0.00005
  was flagged even for an `alpha` of 1e-30. The decision now uses the
  unrounded p-values, and the tables still show four decimals.

- `rdif()` computed `RDIF_RS` from arbitrary values, including negative
  chi-square values, when the covariance matrix of `RDIF_R` and `RDIF_S`
  was singular, which happens when all examinees who responded to an
  item have the same probability of a correct response (a tiny constant
  was added to every entry of the matrix). The statistic is now computed
  with a generalized inverse and the rank of the matrix as degrees of
  freedom, and a warning names the items. Results for regular data are
  unchanged.

- `crdif()` compared `RDIF_R-CR` and `RDIF_RS-CR` of polytomous items
  with chi-square distributions with K and 2K degrees of freedom,
  although the category residuals of each examinee sum to zero, so the
  covariance matrices have ranks K - 1 and 2K - 1. The tests were
  therefore conservative (about 2% to 3% of the DIF-free items were
  flagged by `RDIF_R-CR` and `RDIF_RS-CR` at `alpha = 0.05`). The
  degrees of freedom are now K - 1, K, and 2K - 1, and the statistics
  are computed without adding a small constant to a singular covariance
  matrix, which changed them in the fourth decimal. Results for
  dichotomous items are unchanged.

- The purification of `crdif()` removed the item with the smallest
  p-value, which was the item listed first when several p-values were
  too small to be represented (they were 0). The purification of
  `rdif()` and `crdif()` now removes the flagged item with the smallest
  p-value of the `purify.by` statistic, compared on the log scale, so
  the choice also holds for an item of `rdif()` whose `RDIF_RS` has
  fewer degrees of freedom because of a singular covariance matrix.
  Results for regular data are unchanged.

- `rdif()` and `crdif()` ignored `min.resp` in the initial analysis when
  `score` was supplied. `min.resp` is now applied from the initial
  analysis, with or without `score`: the scores of the examinees with
  fewer than `min.resp` responses are set to missing and these examinees
  are excluded.

- `grdif()` stopped with the error "object 'chisq' not found" when a
  group had no response to an item. The three GRDIF statistics and their
  p-values of such an item are now `NaN`, and the item is not flagged.

- `grdif()` stopped with an obscure error when an examinee had a missing
  ability estimate, and `catsib()` left such an examinee out without a
  message. Both functions now exclude the examinee with a warning, and
  stop with a message when no examinee has an ability estimate.

- `grdif()` put the examinees whose group is missing into the reference
  group, while its post-hoc analyses left them out. They are now left
  out of every analysis, as in the other DIF functions.

- The purification of `grdif()` and `catsib()` stopped with an obscure
  error when one item was left, so it could not end when every item was
  flagged. It now stops with a warning and `complete = FALSE` when every
  item is flagged or when a group has no ability estimate left
  (`catsib()`: fewer than two within `range`, which happens when a
  single item is left and `method = "ML"`), and the remaining items keep
  the results of the last counted iteration, as in `rdif()` and
  `crdif()`. The flagged item with the smallest p-value is removed at
  each iteration, as in `rdif()` and `crdif()`; for regular data this is
  the item removed before. One warning reports the examinees who lost
  their ability estimates during the purification, the warning for an
  incomplete purification reads as that of `rdif()`, and `n.iter` of the
  `with_purify` results is the integer 0 when no item is flagged in the
  initial analysis.

- `grdif()` and `catsib()` silently computed wrong statistics for some
  invalid inputs, such as character responses, responses outside 0 and
  1, or a `group`, `score`, or `se` whose length does not match `data`,
  and stopped with an obscure message for others, such as a `focal.name`
  that is not found in `group` or `method = "EAP.SUM"` or `"INV.TCC"`.
  They now stop with a clear message, also for an `alpha` outside (0,
  1), a `max.iter` that is not a positive whole number when
  `purify = TRUE`, or an invalid `min.resp`. The `est_item` method of
  `grdif()` stops when `score` is passed, which it ignored before, and
  the `est_irt` and `est_item` methods stop when `data` or `D` is
  passed. `catsib()` also stops when `min.binsize` is not an integer of
  2 or more (with 1, the p-values of all items were missing), when
  `n.bin` is not the maximum and the minimum number of bins in this
  order, and when `max.del` is not from 0 to less than 1.

- A logical `item.skip` of `catsib()` skipped the wrong items during the
  purification, so a skipped item could be tested, flagged, and removed.
  It is now read as item positions, also during the purification.

- `grdif()` and `catsib()` compared the p-values rounded to four
  decimals with `alpha`, so an item whose p-value was within rounding of
  `alpha` could be flagged wrongly, and the purification could remove
  it. The unrounded p-values are now compared; the tables still show
  four decimals.

- `grdif()` computed `GRDIF_R`, `GRDIF_S`, and `GRDIF_RS` after adding
  small constants to the covariance matrix of an item when the matrix
  was singular or nearly singular, for example when all examinees of a
  group had the same ability estimate, and could stop with an error. The
  statistics now use a generalized inverse, the degrees of freedom are
  the rank of the matrix instead of G - 1, G - 1, and 2(G - 1), and one
  warning names the items. Results for regular data are unchanged.

- `catsib()` returned a statistic of 0, a standard error of 0, and a
  p-value of `NaN` for an item when the corrected ability estimates of
  the two groups did not overlap, and ended the purification as
  complete. The statistics of such an item are now `NA`, one warning
  names the items, and the item is not flagged. The purification stops
  with a warning and `complete = FALSE` when this happens to every
  remaining item.

- The estimation of GPCM items with large slopes in `est_irt()`,
  `est_mg()`, and `est_item()` rescaled all cumulative logits by a
  common factor once one exceeded 700 in absolute value, which distorted
  the objective function, and the gradient and the Hessian could become
  `NaN` even for logits well below 700. Each row is now shifted by its
  own maximum when a logit exceeds 100 in absolute value; results below
  that range are unchanged.

- `est_irt()` and `est_mg()` computed the standard errors of GRM and
  GPCM items with `D = 1` whatever `D` was. With `D` other than 1, such
  as 1.702, the standard errors of polytomous items were about 1.4 to 2
  times too large, those of the dichotomous items of the same test were
  also affected, and so were `covariance` and `test.2`. The estimates
  were not affected.

- With `fix.g = TRUE` and a `g.val` other than 0.2, `est_irt()` and
  `est_mg()` computed the standard errors of 3PLM items as if the
  guessing parameter were 0.2. The error grows with the distance between
  `g.val` and 0.2, for example about 9 percent with `g.val = 0.25` and
  about 25 percent with 0.05 or 0.35 for the `LSAT6` data. The estimates
  were not affected.

- When `D` was not 1, the Hessian of 1PLM items with a common estimated
  slope (`fix.a.1pl = FALSE`) had a wrong cross term between the slope
  and the difficulties. In `est_item()` the standard error of the slope
  came out too small and the optimization could stop short of the
  maximum; in `est_irt()` and `est_mg()` the M-step followed a different
  path, which changed the number of EM cycles and the estimates within
  the convergence tolerance. The Hessian is now correct, and results for
  `D = 1` are unchanged.

- `covirt()` underestimated the standard errors of polytomous items by
  about 10 to 25 percent, because it computed their information with
  Equation A8 of Li and Lissitz (2004), which treats the category
  indicators of a response as independent Bernoulli variables. It now
  uses the information of the multinomial response, whose standard
  errors agreed with the empirical standard deviations of the estimates
  in a simulation; see Details of `covirt()`. The results for
  dichotomous items do not change.

- `covirt()` stopped with an error for a GRM or GPCM item with two score
  categories. Such an item is now treated as a 2PLM item, and a GPCM
  item listed in `pcm.loc` as a 1PLM item with a fixed slope.

- `covirt()` printed a message instead of issuing a warning when the
  information matrix of an item was singular, and the matrix and the
  vector of standard errors filled with `NA` had no names. It now issues
  a warning, and the `NA` values carry the parameter names.

- In a test with many items, such as several hundred, the likelihood of
  an examinee could underflow to 0 at every quadrature point.
  `est_irt()` and `est_mg()` then stopped with "NA/NaN gradient
  evaluation", and `pcd2()` returned `NaN` for every item. The
  likelihood of such an examinee is now rescaled by its largest
  log-likelihood, and the marginal log-likelihood adds the scale back.
  Results for tests where no likelihood underflows are unchanged.
  `pcd2()` also no longer returns `NaN` when a quadrature point has no
  posterior mass.

- `est_irt()`, `est_mg()`, and `est_item()` could return the thresholds
  of a GRM item in a decreasing order, which is not a valid solution
  because the category probabilities then sum to more than 1 and the
  reported log-likelihood is inflated. It happened mostly in small
  samples or when a score category has few responses (for example, in 13
  of 40 simulated data sets with five categories and 300 examinees).
  When the unconstrained estimates are not in increasing order, the item
  is now estimated again with its thresholds kept in increasing order.
  Results of items whose thresholds were already ordered are unchanged.
  Fits with a GRM item whose thresholds were reversed should be rerun.

- When the second column of `weights` did not sum to 1, `est_irt()` and
  `est_mg()` reported a log-likelihood, AIC, and BIC that were off by
  `N` times the log of the sum, and a variance in `group.par` that was
  multiplied by the sum. With `EmpHist = TRUE` the empirical histogram
  was rescaled to that wrong variance, which changed the scale of the
  item parameters (by up to 1.8 in an example with a sum of 4). The
  weights are now rescaled to sum to 1 with the warning "The weights do
  not sum to 1 and were rescaled to sum to 1.", as in `cac_lee()` and
  `cac_rud()`; weights that are negative, not finite, or sum to zero
  stop with an error. Weights that sum to 1 give the same results as
  before.

- With `EmpHist = TRUE` and user `weights` whose number of quadrature
  points is not 49 (the default of `Quadrature[1]`), `est_irt()` and
  `est_mg()` stopped with an error such as "non-conformable arguments"
  when there were fewer points, and extrapolated the upper tail of the
  empirical histogram from the wrong point when there were more. The
  number of points is now taken from `weights`. The results for 49
  points are unchanged.

- `est_irt()` with `fipc = TRUE` counted a fixed item twice in the first
  E-step without notice when `fix.loc` listed it more than once, which
  changed the estimates of the items to calibrate (by up to 0.45 in an
  example with `fipc.method = "OEM"`). It now stops with an error.

- With `fipc = TRUE`, all items fixed, and user `weights`, the change in
  the first EM cycle was measured from the mean and variance of the
  default distribution instead of those of `weights`, which could shift
  the cycle at which the iteration stopped. It is now measured from the
  moments of `weights`.

- `npar.est`, `aic`, and `bic` of `est_irt()` and `est_mg()` did not
  count the weights of the empirical histogram (`EmpHist = TRUE`) as
  free parameters, and `est_irt()` with FIPC counted only the mean and
  variance of the latent distribution (2 values) even with
  `EmpHist = TRUE`. With 49 quadrature points, the AIC of a model with
  an empirical histogram was 92 too small (2 times the 46 weights),
  which made it look better than a model with a normal distribution. The
  count now includes the free weights: Q - 3 for a histogram rescaled to
  a fixed mean and variance (`est_irt()` without FIPC, and the reference
  groups of `est_mg()`) and Q - 1 for a histogram that is not rescaled
  (`est_irt()` with FIPC, and the free groups of `est_mg()`), where Q is
  the number of quadrature points. The values for a normal distribution
  do not change.

- `est_mg()` with `x = NULL` stopped with an error when `cats = NULL`,
  although the help page allows it when all models are dichotomous, and
  when `model` included `"DRM"`. Both now work, as in `est_irt()`.

- `est_mg()` counted a group that was listed twice in `free.group` twice
  in `npar.est`, which made `aic` and `bic` too large. It silently
  treated a group name in `free.group` that is not in `group.name` as a
  reference group, and it silently dropped a position in `fix.loc` that
  is outside the rows of `x`, so that a different set of items was
  fixed. A repeated group is now counted once, and the other two cases
  stop with an error.

- `est_mg()` ran with a response such as 0.5, which the likelihood and
  the frequency tables treated differently. Responses must now be whole
  numbers from 0 to the number of score categories minus 1, as in
  `est_irt()`.

- With `missing` supplied, the `group` element of `data` returned by
  `est_mg()` still held the missing-value code, while `overall` had
  `NA`. Both now have `NA`.

- `est_item()` stopped with an unclear error when `score` was a data
  frame, contained a missing value, or had a length different from the
  number of rows of `data`, and when only one item had responses. A data
  frame `score` now works, and the other cases stop with a clear message
  or work. The message about items that did not converge numbered the
  items after the items without responses had been removed; it now uses
  the positions in the supplied data.

- `est_mg()` with `fipc = TRUE` always stopped with an error when
  `item.id` was supplied: with "do not share common items" when the IDs
  were those of `x`, and with "At least one model name is mis-specified"
  when they were different (present in 1.3.1). A supplied `item.id` now
  replaces the IDs in `x`, as the help page says, and `fix.id` refers to
  the new IDs.

- FIPC with fixed items whose responses are all the same (for example,
  all correct) produced an error ("NA/NaN gradient evaluation") with
  `fipc.method = "MEM"`, degenerate estimates (a slope of 0.01 and a
  difficulty above 10) with `"OEM"`, and `NaN` group parameters after
  500 cycles when all items were fixed. Such fixed items carry no
  information to place the ability scale. `est_irt()` and `est_mg()` now
  stop with a message when the responses to all fixed items show no
  variation (`est_mg()` uses the responses pooled over the groups). A
  failed short search for starting values in the first two EM cycles is
  now skipped, and when the last search of an item in an M-step fails,
  the item keeps its previous values and is reported as not converged;
  both cases stopped the run before.

- For `est_mg()` with FIPC and all items fixed, `test.2` read "Solution
  is a possible local maximum." although no information matrix is
  computed. It now reads as in `est_irt()`: "No item parameters are
  estimated, so the information matrix is not computed."

- `est_item()` returned the log-likelihood plus the log prior densities
  as `loglikelihood` when prior distributions were used, for example for
  3PLM items with the default `use.gprior = TRUE`, although the help
  page and the printed "-2loglikelihood" call it the log-likelihood. It
  now returns the log-likelihood at the final estimates without the log
  prior densities, as `est_irt()` and `est_mg()` do. The estimates,
  standard errors, and covariances are unchanged, and so is
  `loglikelihood` when no prior is used.

- `bring.mirt()` returned `NA` parameters, or silently wrong ones, for
  some mixed-format fits (present in 1.3.1). The columns of the
  coefficient matrix of mirt follow the model of the first item, so a
  fit that began with a polytomous item gave `NA` for its dichotomous
  items, and a fit that used the graded or gpcm model for every item,
  including two-category items, gave wrong values (the category
  probabilities differed by up to 0.5). The parameters of each item are
  now read by the column names. The intercept `c` of the gpcmIRT model,
  which was dropped, is moved into the thresholds (`b - c / a`); it
  changes the thresholds only when `c` is not 0. A dichotomous item of
  an unsupported model in a mixed-format fit, such as a 4PL item, now
  stops with an error instead of losing its upper asymptote. A model
  fitted to multiple groups, which is not supported, now stops with a
  message that asks for one extracted group (`mirt::extract.group()`).

- `write.flexmirt()` wrapped the lines of the file at the width of the
  console (80 characters in `Rscript`), which broke the file:
  `bring.flexmirt()` could not read it, and flexMIRT would not either.
  This happened for the lines of items with several categories (GRM and
  GPCM), so the files written by the examples of the help page were
  affected. Each item and each group is now written as one line, and the
  file does not depend on the console width.

- `write.flexmirt()` wrote a wrong first step parameter for GPCM items
  with a number of categories other than five (the difference was up to
  0.69 for the four-category items of `simCAT_MX`), because the mean of
  the step parameters was computed with a constant 4 instead of the
  number of steps. With `rePar = TRUE`, it also left out the last
  scoring function value of GPCM items. Files written for GPCM items
  with five categories are unchanged. A GPCM item with two categories
  and `rePar = FALSE` stopped with a subscript error; it is now written
  as a two-parameter item. The file is also closed when the function
  stops with an error.

- `est_irt()` stopped with an out-of-range error for response columns
  stored as character, such as `"0"` and `"1"`, and `est_mg()` stopped
  with an error when only some of its groups had character columns
  (present in 1.3.1). Character and factor columns are now read by their
  values in `est_irt()` and `est_mg()`, as in `est_item()`. A
  non-numeric missing code, such as `missing = "."`, is recoded before
  the columns are read in `est_irt()`, as in `est_mg()`.

## Changes in Behavior

- `catsib()` has a new argument `rho.screen`, `TRUE` by default. When
  the ability estimates of a few examinees have very large standard
  errors, the estimated reliability (rho^2) of the regression correction
  can be far below the reliability of the other estimates, or negative,
  which shrinks the corrected abilities too much and inflates the Type I
  error rate. With `rho.screen = TRUE`, the examinees whose squared
  standard error is equal to or greater than the variance of the ability
  estimates of their group are left out of the computation of the group
  mean, the variance, and the error variance used in the regression
  correction; they are still matched and tested. This changes the
  default results of `catsib()`. The warning about rho^2 falling below
  0.05 is given less often, and once per call instead of at every
  purification iteration. `rho.screen = FALSE` computes rho^2 as earlier
  versions did.

- The post-hoc pairwise analyses of `grdif()` take the reference group
  as the reference group of every pair that includes it, and, for a pair
  of focal groups, the group listed earlier in `focal.name`. Before, the
  group that came later in the sorted order of the group values was the
  focal group of a pair. When the two orders differ, `group.pair` is
  labeled in the new order, the signs of `rdifr`, `z.rdifr`, `rdifs`,
  and `z.rdifs` are reversed, and `n.ref` and `n.foc` are swapped.
  `rdifrs`, the p-values, and the flags do not change.

- The default of `range` in `grdif()` is `c(-5, 5)`, as in `rdif()`,
  `crdif()`, `ripd()`, `catsib()`, and `est_score()`, instead of
  `c(-4, 4)`. The ability estimates of the examinees whose ML estimates
  reach a limit of the range (for example, all responses correct or all
  wrong) change, and so do the results that depend on them.

- `grdif()` supports dichotomous items only and stops with an error when
  `x` includes a polytomous item or any item of the "GRM" or "GPCM"
  model, even one with two score categories, because the GRDIF framework
  has been evaluated for dichotomous items (Lim et al., 2024). Before,
  polytomous items were analyzed as net DIF, and the default post-hoc
  analyses of mixed-format data could stop with an error.

- `grdif()` stops with an error when the examinees whose group is not in
  `focal.name` belong to more than one group, because they form the
  reference group. Before, the GRDIF statistics combined them into one
  reference group while the post-hoc analyses compared them separately.

- `grdif()` and `catsib()` apply `min.resp` from the initial analysis
  also when `score` is supplied: the scores of the examinees with fewer
  than `min.resp` responses are set to missing and these examinees are
  excluded, and a warning reports their number.

- `est_irt()` stops with a clear message that names the items when a
  response is not a whole number from 0 to the number of score
  categories minus 1. Before, it stopped with "subscript out of bounds".
  With `fipc = TRUE` it also stops when `x` and `data` have different
  numbers of items, when `fix.loc` is empty or points outside the rows
  of `x`, and when an ID in `fix.id` is not found in the item IDs.
  Before, an unmatched ID in `fix.id` was dropped without notice.

- A supplied `item.id` replaces the IDs in `x` in `est_irt()` also when
  `fipc = FALSE`, as it already did with `fipc = TRUE`. Before, it was
  ignored when `x` was supplied. It must have one value per item. The
  same holds for `est_mg()` (one vector of IDs per group, also without
  FIPC) and `est_item()`.

- For `est_irt()` with FIPC and all items fixed, `test.2` reads "No item
  parameters are estimated, so the information matrix is not computed.",
  instead of the message about a possible local maximum.

- `est_irt()`, `est_mg()`, and `est_item()` issue one warning, with the
  item IDs, when an item to be estimated has a score category without
  responses or all its responses are the same. The maximum likelihood
  estimate does not exist for such an item, and the functions returned
  extreme or starting values with standard errors of 99999 without
  notice. The estimates are not changed. `est_mg()` judges an item on
  its responses pooled over the groups, and `est_item()` judges the
  items that have responses.

- `bring.mirt()` labels a two-category item fitted with the graded,
  gpcm, or gpcmIRT model as `"DRM"`, with a guessing parameter of 0,
  because a two-category item is a dichotomous item in irtQ. Such an
  item was labeled `"GRM"` or `"GPCM"`, with wrong or missing parameters
  (see Bug Fixes).

- The `type` argument of `bring.flexmirt()`, `bring.bilog()`, and
  `bring.parscale()` is `"par"` when omitted, as the usage line
  suggests, and its case is still ignored. Omitting it stopped with "the
  condition has length \> 1" before. An unknown type stops with the
  available choices.

- `bring.parscale()` stops with an error for the model codes 2, 4, 6,
  and 8 of PARSCALE when the file has polytomous items. The thresholds
  of such items were all `NA`.

## Minor Improvements

- `est_mg()` checks its inputs and stops with a clear message when
  `group.name` does not have one unique name per group, when `x`,
  `model`, `cats`, or `item.id` does not match the groups or the columns
  of `data`, when an item ID is repeated within a group, when items with
  the same ID have different models or numbers of score categories in
  different groups, when `free.group` names or numbers a group that does
  not exist or leaves no reference group without FIPC, when `fix.loc` is
  not a list with one element per group or has repeated positions, when
  an ID in `fix.id` is not found, and when `fipc.method` is not `"MEM"`
  or `"OEM"` (also in `est_irt()`).

- `est_mg()` with FIPC warns when a fixed item has different parameter
  values in different groups. The values of the first group that
  includes the item are used, as before. Differences of at most 1e-10
  are ignored.

- `run_flexmirt()` stops with a clear message on systems other than
  Windows, where flexMIRT does not run, and when "C:/Program Files"
  holds more than one flexMIRT folder (give `dir.flex` then).
  `write.flexmirt()` stops with a clear message when `file` is `NULL`.

- `getirt()` can extract the matched `call` of `est_irt()`, `est_mg()`,
  and `est_item()` objects with `what = "call"`.

- `ctt()` records the `correct` setting in its `crit` element, and the
  `summary()` report names the item-total correlation used for flagging.

- `ctt()` and `ctt_distr()` stop with a clear message when `crit.p`,
  `crit.dis`, `crit.distractor`, or `total` is not valid, when fewer
  than two examinees remain, or when `ctt_distr()` receives numeric and
  non-numeric keys together.

- Inputs that were silently misread now stop with a clear message.
  Functions that take item metadata reject parameter columns stored as
  text, a missing or too small `cats`, a dichotomous model with `cats`
  other than 2, and a polytomous item with more thresholds than
  `cats - 1`, including any value in `par.3` of a two-category GRM or
  GPCM item. `shape_df()` and `simdat()` reject parameter vectors that
  do not match the items, `shape_df_fipc()` checks `fix.loc`,
  `gen.weight()` checks `dist` and `theta`, and `prm()` rejects an
  unknown `pr.model`.

- `prm()` uses the GRM when `pr.model` is omitted, `simdat()` repeats a
  single `pr.model` for all polytomous items, `shape_df()` accepts
  `par.drm` without `g`, and item metadata without a `par.3` column is
  accepted when every item is 1PLM or 2PLM.

- `plot.info()` and `plot.traceline()` show item IDs as given
  (`"item-2"` was shown as `item.2`), and the panel titles of
  `plot.traceline()` read `Score: 0` instead of `Score: resp.0`.

- The warning for `"DRM"` items reads "All 'DRM' items are treated as
  '3PLM' items." and no longer mentions item parameter estimation,
  because it also appears in functions such as `info()` and `simdat()`.
  `est_mg()` gives it once, not once for each group.

- `est_score()` stops with a clear message when `data` has a number of
  columns other than the number of items (it scored the wrong items
  before), when a response is not an integer within the categories of
  its item, and when `method`, `stval.opt`, `tol`, or `max.iter` is not
  valid (a `max.iter` below 1 or not a whole number removed the
  iteration limit). For an `est_irt` object, supplying `data` or `D` is
  now an error instead of being ignored.

- `est_score()` now has a method for `est_item` objects, as its
  documentation already stated.

- `irtfit()` and `sx2_fit()` stop with a clear message when `data` has a
  number of columns other than the number of items or a response is not
  an integer within the categories of its item (such a response was
  counted in the wrong cells before), and they read character and factor
  responses as numbers instead of level codes. `irtfit()` also stops for
  a `score` of the wrong length, scores without variation, an invalid
  `range.score`, and a `loc.theta` other than `"average"` or `"middle"`
  (ignoring case; any other value was treated as `"average"`), leaves
  out examinees with a missing ability estimate instead of failing, and
  works when only one item remains after items with fewer than two
  responses are excluded. `plot.irtfit()` stops for a `type` other than
  `"both"`, `"icc"`, or `"sr"` and for an `item.loc` that is not an
  evaluated item. `sx2_fit()` stops for a test with a single item, which
  has no S-X2 statistic. The warnings for excluded items and for missing
  responses in `sx2_fit()` name the items, and the error for a response
  outside the categories names the item and its column.

- `bisection()` warns when the function values at `lb` and `ub` have the
  same sign, instead of returning a bound without notice, and it stops
  after exactly `max.it` iterations instead of `max.it + 1`. Its help
  page names the last returned element `delta`, as in the code.

- `est_score()` with `method = "INV.TCC"` stops with a clear message
  when the test characteristic curve does not reach a sum score within
  theta from -20 to 20, for example for items with very small slopes or
  guessing parameters that sum to nearly the maximum score.

- `ctt()` computes alpha with each item removed faster for large tests.

- `cac_lee()` and `cac_rud()` leave out examinees with a missing ability
  estimate (or standard error) with a warning. Before, `cac_lee()` and
  `cac_rud()` stopped with an obscure error for a missing ability
  estimate, and `cac_rud()` returned wrong indices for a missing
  standard error. Both stop for cut scores that are not in strictly
  ascending order, and `cac_rud()` stops for a standard error that is
  not positive. `cac_lee()` no longer stops when a performance level
  contains no integer summed score, and it is much faster with ability
  estimates.

- `panel_info()` stops with a clear message for a route map that does
  not define a panel (a module without transitions, a cycle, or a module
  outside the stages) instead of running forever, and it finds every
  stage when the modules are not numbered in the order of the stages.
  `panel_info()`, `run_mst()`, `reval_mst()`, and `find_cut()` read a
  logical route map as a 0/1 matrix.

- `reval_mst()` accepts item metadata with a missing guessing parameter,
  such as 2PLM items from `shape_df()`. It checks and completes the item
  metadata as `run_mst()` does, so `item.by.mod` and `item.by.path` hold
  plain data frames with upper-case model names and without parameter
  columns that are all `NA`. `reval_mst()` and `run_mst()` stop with a
  clear message when `cut_score` does not match the panel or `module`
  has a number of columns other than the number of modules, they read a
  `NULL` element of `cut_score` for a stage with a single module as an
  empty vector, and they assign the cut scores to the modules of a stage
  in the order of the module indices. `reval_mst()` also stops with a
  clear message when a sum score has no inverse TCC estimate (with
  `intpol = FALSE` or a `range.tcc` that does not cover the estimates);
  before, it stopped with an obscure error or returned `NA` for every
  row of `eval.tb`.

- `run_mst()` stops with a clear message when `tol` or `max.iter` in
  `route_score` or `final_score` is not valid (a `max.iter` of 0 or 2.5
  never ended), and it warns once about examinees without any observed
  response to the administered items.

- `find_cut()` no longer warns about two proper crossings when a
  crossing lies exactly on a grid point.

- `ripd()` results have a `print()` method that prints the flagged items
  and the RIPD statistics of the analyses without and with purification;
  `what` selects `"no_purify"` or `"with_purify"`. Before, `print()`
  listed the whole object, including the ability estimates of every
  examinee.

- The `est_irt` and `est_item` methods of `rdif()` and `crdif()` stop
  with a clear message when `data` or `D` is passed, because both are
  taken from the object, and the `est_item` methods do so for `score` as
  well. Before, they were silently ignored or caused an obscure error.

- The purification of `rdif()`, `crdif()`, `ripd()`, `grdif()`, and
  `catsib()` is faster: after the first iteration, only the examinees
  who responded to the removed item are rescored when `ncore = 1` and
  `method` is a pattern scoring method. Results are unchanged.

- The purification of `rdif()`, `crdif()`, `ripd()`, `grdif()`, and
  `catsib()` gives the warning of `est_score()` about `ncore > 1` once
  instead of at every iteration.

- `rdif()`, `crdif()`, and `ripd()` stop with a clear message for
  `method = "EAP.SUM"` or `"INV.TCC"` when they estimate the abilities
  (`score = NULL` or `purify = TRUE`), because these methods give no
  ability estimate for each response pattern. Before, `rdif()` and
  `crdif()` stopped with an obscure message. `ripd()` and `pcd2()` also
  stop with a clear message when `min.resp` is not a single nonnegative
  number, and `ripd()` stops when `score` is not numeric, when `alpha`
  is not a single number between 0 and 1 (an `alpha` of 1 or more, or a
  text value, flagged every item), or when `max.iter` is not a positive
  whole number and `purify = TRUE`. The `est_irt` and `est_item` methods
  of `ripd()` stop with a clear message when `data` or `D` is passed,
  and the `est_item` method when `score` is passed, as those of `rdif()`
  and `crdif()` do.

- `pcd2()` stops with a clear message when `max.iter` is not a positive
  whole number and `purify = TRUE`, and when `crit.val` is not `NULL` or
  a single nonnegative number. Before, it accepted `max.iter = 2.5`, and
  a text `crit.val` flagged no item without a message.

- The `print()` methods of `grdif()` and `catsib()` leave the
  significance symbol of an item with a missing p-value blank, as those
  of `rdif()` and `crdif()` do for an item skipped with `item.skip`.
  Before, it showed `?`, which the legend does not explain.

- The `print()` methods of `rdif()` and `crdif()` leave the significance
  symbol of an item skipped with `item.skip` blank. Before, it showed
  `?`, which the legend does not explain.

- The `dif_item` elements `rdifr` and `rdifs` of the `no_purify` results
  of `rdif()` are integer vectors like `rdifrs`, and `n.iter` of the
  `with_purify` results of `rdif()` and `crdif()` is the integer 0 when
  no item is flagged in the initial analysis. `npar.est` of `est_mg()`
  and `n.iter` of the `with_purify` results of `catsib()` are integers
  as well. The values are unchanged.

- When the Hessian of an item is singular during the estimation of
  `est_irt()`, `est_mg()`, and `est_item()`, it is now repaired by
  adding a ridge to its diagonal instead of random numbers. The results
  are deterministic, and the random number stream of the user is no
  longer used by the estimation.

- The help pages of `est_irt()`, `est_mg()`, `est_item()`, and
  `getirt()` now describe how the standard errors are computed: the
  information of the item parameter priors is added to the cross-product
  information matrix of `est_irt()` and `est_mg()`, and `est_item()`
  uses the observed information with the abilities treated as known. The
  help pages of `est_irt()`, `est_mg()`, and `getirt()` also state that
  `test.2` checks whether the information matrix is singular, and those
  of `est_irt()`, `est_mg()`, and `est_item()` state that the slope
  prior `aprior` is placed on `D` times the slope.

- `est_irt()` and `est_mg()` return a new last element, `npar.detail`,
  an integer vector with the numbers of free item parameters (`item`),
  of means and variances of the latent distributions (`dist`), and of
  weights of the empirical histograms (`weight`); their sum is
  `npar.est`. The `summary()` output shows these three numbers under
  "Number of free parameters", and `getirt()` extracts it with
  `what = "npar.detail"`. Objects created by earlier versions have no
  `npar.detail`, and their summary shows the total only.

- Corrected and extended the help page of `est_irt()`. The returned
  `loglikelihood` is the sum over examinees of the marginal
  log-likelihood of the response pattern, without the prior densities
  (it is not a sum over items); `ncase` and `posterior.dist` refer to
  examinees, not response patterns; `data` is a matrix; and the returned
  `call` is listed. The help page also describes the response coding (0
  to the number of categories minus 1), the role of `x`, `item.id`,
  `weights`, `group.mean`, and `group.var` with and without FIPC, the EM
  algorithm, the empirical histogram, the two FIPC methods (including
  that the returned distribution of `fipc.method = "OEM"` is updated
  once afterward and that the standard errors in `group.par` for FIPC
  are approximate), the counting of the free parameters, the bounds of
  the slopes, and the items without a maximum likelihood estimate.
  `Etol` and `MaxE` are described as the criterion and the maximum
  number of EM cycles, and the argument tables, the lists in the
  argument descriptions, and the examples were corrected. The same
  corrections apply to `getirt()` and to the arguments shared with other
  help pages.

- Corrected and extended the help pages of `est_mg()`, `est_item()`, and
  `getirt()` in the same way. The help page of `est_mg()` now explains
  `free.group`, the roles of `weights`, `group.mean`, `group.var`, and
  `EmpHist` for groups in and not in `free.group`, the MG-FIPC algorithm
  (Kim & Kolen, 2016), and the structure of the returned
  `loglikelihood`, `group.par`, `posterior.dist`, `data`, `ncase`, and
  `nitem`. The help page of `est_item()` states that the estimates are
  Bayes modal when priors are used, that the items are estimated
  separately, and what `group.par`, `data`, `deleted.item`, and
  `n.response` contain. The `README` and the reference of Kim and
  Kolen (2016) were corrected.

- The help pages of `est_irt()` and `est_mg()` state that the EM cycles
  stop on the largest absolute change in the estimated item parameters
  (with all items fixed in FIPC, in the means and variances of the
  latent distributions). The estimated means and variances of the latent
  distributions and the weights of empirical histograms are not part of
  this criterion, so they can be less converged than the item
  parameters, and a smaller `Etol` gives more precise values; the
  default `Etol` of `est_mg()` (`1e-3`) is larger than that of
  `est_irt()` (`1e-4`). The help pages also state that, with
  `EmpHist = TRUE`, `group.par` reports the fixed mean and variance of a
  rescaled histogram, which can differ slightly from the mean and
  variance of the returned `weights`, and that the standard errors of
  the item parameters treat the latent distributions as known.

- Corrected the help pages of `bring.flexmirt()` and the other import
  functions (the returned elements, the arguments `type`, `rePar`, and
  `rePar.gpc`, the model labels, and the item types of mirt that are
  supported), `write.flexmirt()` (`norm.pop` holds the mean and the
  variance, metadata from irtQ needs `rePar = FALSE`, and the scaling
  constant), `run_flexmirt()` (Windows only), `summary()` (the returned
  object), `sx2_fit()` (missing responses are replaced with 0), and the
  data sets `LSAT6`, `simCAT_DC`, `simCAT_MX`, `simIPD`, `simMG`, and
  `simMST`. The package overview gives `D = 1` as the default and lists
  the functions for item parameter drift, MST panels, and the import and
  export of files. The help pages of `rdif()`, `crdif()`, `grdif()`,
  `catsib()`, `est_irt()`, `est_mg()`, and `est_item()` state that the
  `digits` argument of their print methods is not used.

- The print methods are now documented on the help pages of the
  functions that create the objects (`rdif()`, `crdif()`, `grdif()`,
  `catsib()`, `irtfit()`, `run_mst()`, and `find_cut()`), and those of
  the `est_irt()`, `est_mg()`, and `est_item()` objects and their
  summaries on the help page of `summary()`, including the use of their
  `digits` argument.

# irtQ 1.3.1

## Bug Fixes

- Fixed the routing in `run_mst()`. The ability estimate used to choose
  the next module was computed from the responses to the current module
  only. It is now computed from the responses to all modules
  administered so far, as intended. This applies to every routing method
  (`"bmat"`, `"mfi"`, and cut scores) and every `route_score` method, so
  simulated paths and results can differ from earlier versions, in
  particular from stage 3 onward.

- Fixed cut-score routing in `run_mst()` for route maps in which a
  module reaches only some modules of the next stage. The examinee could
  be sent to a module that did not match the cut scores. Only the cut
  scores that separate the reachable modules are now used, which is the
  rule applied by `reval_mst()`. Results are unchanged when every module
  of a stage can be reached.

- Fixed the final estimate of `run_mst()` for responses with missing
  values (`response` argument). The ML, WL, MLF, MAP, and EAP estimates
  used the parameters of other items instead of those of the observed
  items.

- `score_resp()` and `ctt_distr()` now match a data frame `key` to the
  item columns by item number when its `item` column is character or
  factor. Before, keys for tests with 10 or more items could be assigned
  to the wrong items.

- `shape_df()` with `default.par = TRUE` now repeats a single value of
  `cats` or `model` for all items. Before, a single `cats` value gave
  item IDs of "V1" for every item and a guessing parameter of 0 instead
  of 0.2 for 3PLM items.

- `irtfit()` now computes the observed category proportions directly
  from the frequencies. Before, they came from
  `janitor::adorn_percentages()`, which can include the total column in
  the denominator in some environments and halve the proportions, the
  residuals, and `overSR.prop`. Results are unchanged when the
  proportions were computed correctly.

- `catsib()` no longer overwrites `score` with `se` when `se` is given
  as a matrix or data frame.

- `cac_rud()` no longer stops with a dimnames error when a performance
  level has no examinees, and the label of the total row of `marginal`
  is now "marginal", as in `cac_lee()`.

- `simdat()` treats NA values in `g.drm` as zeros, as the item metadata
  input does. Before, the responses to those items were all NA.

- The Wald confidence intervals in `plot.irtfit()` now use the two-sided
  critical value, `qnorm(1 - alpha / 2)`. Before, they used
  `qnorm(1 - alpha)` and were 90% intervals at the default
  `alpha = 0.05`, while the Wilson intervals were 95% intervals.

## Minor Improvements

- `reval_mst()` now stops with an informative message when the modules
  in a stage differ in maximum sum score.

- In `run_mst()` and `reval_mst()`, an ability estimate equal to a cut
  score is now routed to the higher module, as in the classification
  rule of `cac_lee()` and `cac_rud()`. Results change only when an
  estimate equals a cut score exactly.

- `ctt()`, `ctt_distr()` (scored-category mode), and the CTT helper
  functions now stop with an informative error when an item score is not
  a whole number between 0 and `cats - 1`. When `cats` is inferred,
  every item has at least two categories, so an item that every examinee
  scores 0 on gets a difficulty of 0 and is flagged.

- `ctt()` reports the listwise deletion of incomplete rows once, and the
  item table has default row names.

- `plot.find_cut()` no longer passes an unused `inherit.aes` argument to
  `geom_vline()`, which caused warnings with some ggplot2 versions.

- `catsib()` now stops with an informative message when `score` is
  supplied without `se`. Before, it failed with a "missing value where
  TRUE/FALSE needed" error.

## Documentation

- Stated in `?reval_mst` that the recursion is based on inverse TCC
  estimates, that only inverse TCC scoring with cut-score routing is
  supported, and that modules in a stage must have the same maximum sum
  score; corrected the documented default of `theta` and the description
  of the returned list.

- In `?run_mst`, described the cumulative routing estimate and the
  relation to `reval_mst()`, and stated that `"EAP.SUM"` and `"INV.TCC"`
  count a missing response as 0 in the final sum score.

- Corrected the interpretation of alpha with the item removed in `?ctt`.
  The `ctt()` examples and the CTT article now use simulated IRT data.

- Corrected the description of the expected frequencies in `pcd2()` and
  stated that `crit.val = NULL` flags no items.

- Reworded the motivation of `ripd()`.

- In `?simMST`, noted that the item parameters are on the D = 1.702
  scale. In `?find_cut`, clarified that the cut scores in `simMST` were
  obtained with `find_cut()` and softened the statement about path
  reversals.

- Corrected the probability matrix of the first example in `?lwrc`, the
  GPCM formula note in `?irtQ`, and the class of `prob.cats` in
  `?traceline`. Corrected statements in the README, the vignette
  overview, and the articles (shrinkage of MAP and EAP, `range.score` in
  `irtfit()`, `fix.id` in `est_mg()`, the effect of `EmpHist` in FIPC,
  the fixed-slope 1PLM, fixed guessing, the CATSIB regression
  correction, the purification procedure, and the usage notes of the
  utility functions). The DIF article now simulates item difficulties in
  a narrower range so that the pooled calibrations converge.

# irtQ 1.3.0

## New Features

- New classical test theory (CTT) functions: `ctt()` for item- and
  test-level statistics, with `print()` and `summary()` methods;
  `freq_score()` for a total-score frequency table; `ctt_distr()` for
  option or score-category distributions with point-biserial
  correlations; and `score_resp()` for scoring raw selected-response
  data with an answer key. Item-total and option-total correlations are
  reported both with and without the item.

- Added a vignette, "Introduction to irtQ" (`vignette("irtQ")`), that
  walks through calibration, scoring, and model-fit evaluation.

## Bug Fixes

- Fixed the EM convergence check in `est_irt()` and `est_mg()`. It used
  the largest signed change instead of the largest absolute change, so
  the EM could stop before the `Etol` criterion was met, most often when
  some parameters were held constant (for example, the guessing
  parameter of 1PLM and 2PLM items, or slopes fixed by `fix.a.1pl` or
  `fix.a.gpcm`). Refitting a model can change the estimates and take
  more EM cycles. A fit from an earlier version that reported
  `maxpar.diff` equal to 0 most likely stopped early.

- `est_irt()` and `est_mg()` no longer warn that the convergence
  criteria are not satisfied when `fipc.method = "OEM"` or when the EM
  converges at the `MaxE`-th cycle. In these cases, `test.2` can now
  report a possible local maximum.

- `est_item()` now uses Beta(5, 16) as the default prior for the
  guessing parameter, as in `est_irt()` and `est_mg()`.

- `getirt()` now returns `posterior.dist` and `scale.D` for objects from
  `est_irt()` and `est_mg()`, as documented.

## Minor Improvements

- `find_cut()` now uses `D = 1` by default, like the other functions in
  the package (previously 1.702). Supply `D = 1.702` to reproduce
  earlier results.

- `shape_df_fipc()` now accepts a single value of `cats` or `model` for
  all new items.

- The mirt package moved from Imports to Suggests. Only `bring.mirt()`
  uses it, and it now stops with an informative message when mirt is not
  installed.

## Documentation

- Added examples to `print.ctt()`, `summary.ctt()`,
  `print.summary.ctt()`, and the data sets.

- Corrected the descriptions of `Etol`, `test.1`, `test.2`, and
  `maxpar.diff`, the return class of `est_mg()`, `dif_item` in the
  purified results of `rdif()`, `crdif()`, and `grdif()`, and the
  component names of `simCAT_DC` and `simCAT_MX`.

- Corrected bibliographic details, added DOIs, cited the polytomous RDIF
  extensions in `rdif()` and `crdif()`, and listed every reference cited
  in the help pages on the package help page (`?irtQ`).

# irtQ 1.2.0

## New Features

- Added a new function, `find_cut()`, which identifies TIF-crossing
  routing cut scores for MST panels. For each adjacent module pair
  within a stage, the function locates the theta where the two modules'
  test information functions (TIFs) intersect. A warning is issued when
  the mean difficulty order of modules within a stage differs from their
  input index order.

- Added a new S3 method, `plot.find_cut()`, which visualizes TIF curves
  and routing cut scores stage by stage using **ggplot2** facets.
  Proper, anomalous, and unselected cut scores are distinguished by line
  type. A `layout` argument (`"vertical"` / `"horizontal"`) controls the
  facet orientation.

- Added a new function, `run_mst()`, which simulates MST administrations
  for a given panel structure and returns response data along with
  ability and routing information for each simulated examinee.

- Added a new exported dataset, `simIPD`, which contains simulated CAT
  response data for illustrating IPD detection with `ripd()` and
  `pcd2()`. The dataset represents one replication of a CAT simulation
  (N = 3,000; test length = 30; 360-item 3PLM pool) in which 5% of items
  (18 items) had both $a$ and $b$ parameters decreased by 0.5.

## Minor Improvements

- `est_irt()`, `est_item()`, and `est_mg()` now accept a partial
  `control` list. Users can specify only the arguments they wish to
  override (e.g., `control = list(iter.max = 500)`); unspecified
  arguments fall back to their defaults via `modifyList()`.

- Reorganized the **pkgdown** reference page: `ripd()` and `pcd2()` are
  now grouped under a new *Item Parameter Drift (IPD)* section, and
  `reval_mst()`, `panel_info()`, `find_cut()`, `plot.find_cut()`, and
  `run_mst()` are grouped under a new *Multistage-Adaptive Test (MST)*
  section.

## Documentation

- Expanded the documentation for `ripd()` by adding a `@details` section
  that covers the theoretical background of the RIPD framework,
  asymptotic distributions of the three RIPD statistics, a drift-type
  diagnostic guide, the CAT-specific three-step workflow, and the
  purification procedure. Also added a `\donttest{}` example
  demonstrating a complete CAT-based IPD detection workflow using the
  `simIPD` dataset. Updated `@references` to include Lim & Choe (2023)
  and replaced the previous conference paper citation with the journal
  reference (Lim & Han, 2026).

- Added a `\donttest{}` example to `pcd2()` demonstrating CAT-based IPD
  detection using the `simIPD` dataset, including the bootstrap critical
  value procedure described in Lim & Han (2026).

- Updated the *MST Panel Evaluation and Simulation* article
  (`vignettes/articles/mst-panel-evaluation.Rmd`) to introduce
  `run_mst()` and extend the existing `reval_mst()` content with routing
  and scoring examples, a `find_cut()`-based principled cut score
  derivation, a side-by-side routing method comparison, and a Monte
  Carlo-vs-analytical validation example (Example 6).

## Bug Fixes

- Rebuilt the `simMST` dataset: the previous version had 9 items
  duplicated across non-adjacent modules because the original assembly
  only enforced no-overlap within a single routing pathway. The new
  version enforces a global no-overlap constraint across all 7 modules
  (56 unique items total) and adds a mean(*b*) per-module band
  constraint; cut scores were regenerated via `find_cut()` on the
  rebuilt modules.

# irtQ 1.1.0

## Major Improvements

- Improved the speed and reduced memory usage of item parameter
  estimation and standard error computation in `est_irt()`,
  `est_item()`, and `est_mg()`.
- Improved the computational speed of `est_score()` by up to 52% for
  dichotomous items (N = 10,000) and up to 27% for mixed-format tests,
  through a series of optimizations.
- Improved the computational speed of `sx2_fit()` substantially by
  replacing the O(J^2) Lord-Wingersky recursion with a forward-backward
  pass (up to 11x faster for mixed-format tests with J = 55 items) and
  vectorizing internal helper functions `expFreq()`, `obsFreq()`, and
  the PRM category-collapsing routine.

## New Features

- Added a unit test suite using the **testthat** 3rd edition
  (`testthat >= 3.0.0`). Tests cover core functions including `drm()`,
  `prm()`, `est_irt()`, `est_score()`, `est_mg()`, `rdif()`, `crdif()`,
  and `catsib()`, with the relevant tests across dichotomous,
  polytomous, and mixed-format item scenarios (355 tests total).
- Added a new function, `ripd()`, which implements the Residual-based
  Item Parameter Drift (RIPD) detection framework. The function computes
  three RIPD statistics: $RIPD_R$, $RIPD_S$, and $RIPD_{RS}$ (one for
  each item). $RIPD_R$ captures uniform item parameter drift (IPD) via
  differences in mean raw residuals between groups, $RIPD_S$ captures
  nonuniform IPD via differences in mean squared residuals, and
  $RIPD_{RS}$ is a combined chi-square-based statistic sensitive to both
  types of drift. An optional purification procedure is also supported.

## New Articles

- Launched the irtQ documentation website at
  <https://hwangQ.github.io/irtQ/>, built with **pkgdown**. The site
  includes a full function reference index and the following vignettes
  covering the complete irtQ workflow:
  - *Getting Started with irtQ*: an end-to-end overview of the package
    workflow.
  - *Item Parameter Estimation*: detailed guidance on `est_irt()`,
    `est_item()`, and `est_mg()`.
  - *Ability Estimation*: scoring methods available in `est_score()`.
  - *Model-Data Fit Evaluation*: using `irtfit()` and `sx2_fit()` to
    assess model fit.
  - *DIF Detection*: applying `rdif()`, `grdif()`, and `catsib()` to
    detect item bias.
  - *Classification Accuracy and Consistency*: computing indices via
    `cac_lee()` and `cac_rud()`.
  - *Utility Functions*: usage of `info()`, `traceline()`, `lwrc()`,
    `simdat()`, and related helpers.
  - *Evaluating MST Panels with `reval_mst()`*: measurement precision
    and bias evaluation for multistage adaptive tests.

## Bug Fixes

- Fixed a minor bug in `sx2_fit()` that caused incorrect cell collapsing
  between two adjacent score categories for polytomous items when
  computing the S-$X^2$ item fit statistic.
- Fixed minor bugs in `est_score()` and `info()`.
- Resolved an issue in `catsib()` where the function failed when all
  responses were missing (NA) in either the reference or focal group.
- Updated `cac_rud()` to include the `x` argument, allowing users to
  pass item metadata data frames directly.
- Revised default `control` parameters in `est_irt()`, `est_item()`, and
  `est_mg()`, and updated the documentation accordingly.
- Fixed a minor bug in `est_score()` function in terms of Newton-Raphson
  method.
- Fixed multiple stability issues in `catsib()`:
  - The final bin exclusion step in `catsib_item()` used a hardcoded
    threshold of 3 instead of the user-supplied `min.binsize` argument,
    causing inconsistent bin filtering behavior.
  - The reliability estimate `rho2` in `catsib_one()` was not clamped to
    $[0, 1]$, so when `errvar > sigma2` (e.g., very few items or
    purification cascade), a negative `rho2` reversed the regression
    correction direction, inflating the Type I error rate.
  - When `errvar >= sigma2` during purification, `rho2` collapsed to 0,
    causing all corrected scores to converge to the group mean. With
    group mean differences (impact), this produced empty bin data frames
    and an invalid purification result. A minimum floor of 0.05 is now
    enforced for `rho2` to preserve score spread.
- Fixed a critical bug in `covirt()` where the guessing parameter
  (`par[,3]`) was incorrectly passed as the difficulty parameter
  (`par[,2]`) to the `integrand()` function for DRM items. This caused
  the gradient computation to receive `c = b`, which zeroed out the
  $\partial P/\partial a$ and $\partial P/\partial b$ gradient
  components and produced a singular Fisher information matrix.
  Additionally, `NA` values in `par.3` for 1PLM and 2PLM items are now
  substituted with 0 prior to gradient evaluation to prevent `NA`
  propagation.

# irtQ 1.0.0

- The documentation for the `irtQ` package has been revised to reflect
  updates to function behavior, fix typos, and provide more relevant and
  detailed information for existing functions.

- A new function, `crdif()`, has been added. This function computes
  three statistics from the residual-based DIF detection framework using
  categorical residuals (RDIF-CR). It allows for the detection of global
  DIF, particularly in polytomously scored items.

- A new function, `shape_df_fipc()`, has been introduced. This function
  merges fixed-item metadata with automatically generated metadata for
  new items and produces a single data frame ordered by test position.
  It is designed to support fixed item parameter calibration (FIPC) via
  the `est_irt()` function.

- The `plot()` method has been enhanced to support the display of all
  item characteristic curves for a given item in a single panel.

- The `pcd2()` function has been updated to include a purification
  procedure.

- The `rdif()` and `catsib()` functions now include an `item.skip`
  argument. This allows users to specify a numeric vector of item
  indices to exclude from the DIF analysis.

# irtQ 0.2.1

- Enhanced functionality of the `bind.fill()` function by adding a new
  argument `fill`. The value in the argument is used to fill in missing
  data when aligning datasets.

- Fixed a bug within the `est_irt()` function that was previously unable
  to implement the fixed item parameter calibration (FIPC) when only
  freely estimating a single item given that all other items are fixed.

- Added a new function, `reval_mst()`, which evaluates the measurement
  precision and bias in Multistage-Adaptive Test (MST) panels using a
  recursion-based evaluation method introduced by Lim et al. (2020).

- Added a new function, `pcd2()`, which computes the Pseudo-count
  $D^{2}$ statistics (Cappaert et al., 2018; Stone, 2000) to detect item
  parameter drift.

# irtQ 0.2.0

- Introduced Warm's (1989) Weighted Likelihood (WL) estimation method to
  the `est_score()` function. This WL scoring method can now be utilized
  by setting `method = "WL"`.

- Enhanced the speed of ability parameter estimation in the
  `est_score()` function when using the ML, MLF, or MAP methods for the
  `method` argument. The updated version performs approximately three
  times faster than its predecessor.

- Addressed a bug within the `est_score()` function that was previously
  unable to accurately compute scores when only a single item data was
  provided. This issue was occurring with the EAP.SUM and INV.TCC
  estimation methods.

- Added two new functions for computing classification accuracy and
  consistency: `cac_rud()` and `cac_lee()`.

  - `cac_rud`: This function implements Rudner's (2001, 2005) method for
    computing classification accuracy and consistency. It takes cut
    scores, ability estimates, standard errors, and optional weights as
    inputs and returns a list containing a confusion matrix, marginal
    and conditional classification accuracy and consistency indices, the
    probability of being assigned to each level category, and the cut
    scores used in the analysis.
  - `cac_lee`: This function implements Lee's (2010) method for
    computing classification accuracy and consistency. It takes a data
    frame containing item metadata, cut scores, optional ability
    estimates, optional weights, a scaling factor, and a logical value
    indicating the cut score metric as inputs. It returns a list similar
    to `cac_rud`.

- Added a new function, `llike_score()`, which computes the
  loglikelihood of ability parameters given the item parameters and
  response data.

- Enhanced functionality of the `rdif()` and `grdif()` functions: Both
  now support the graded response model (GRM) and generalized partial
  credit model (GPCM).

- Fixed an issue in the `grdif()` function that inaccurately calculated
  the GRDIF statistics when group membership was specified in a
  non-standard way. Specifically, the problem arose when 0 wasn't used
  as the reference group and consecutive numbers (e.g., 1, 2, 3) weren't
  used to represent focal groups in the `group` argument.

# irtQ 0.1.1

- Resolved the misalignment issue of standard errors in the output of
  the `est_irt()` function when `fix.a.1pl = TRUE` is specified and the
  items are calibrated using the 1PLM.

- Added a new function, `grdif()`, to perform differential item
  functioning (DIF) analysis across multiple groups. This function
  calculates three generalized IRT residual DIF (GRDIF) statistics. For
  more information about the function and its usage, please refer to the
  accompanying documentation.

- Fixed several typos in the manual documentation

# irtQ 0.1.0

- Initial release on CRAN

- The `irtQ` package is a successor of the `irtplay` package which was
  retracted from R CRAN due to the intellectual property (IP) violation.
  All issues of the IP violation have been clearly resolved in the
  `irtQ` package.

- Most of the functions the `irtQ` package are identical in appearance
  and functionality to those of `irtplay` package except a few functions
  (e.g., `shape_df()`, `est_score()`). However, the computing speed of
  several functions (e.g., `est_irt()`, `est_score()`, `lwrc()`) in the
  `irtQ` package are faster than the previous ones in the `irtplay`
  package. Read the documentation carefully prior to using the
  functions.
