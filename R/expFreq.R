# This function returns a contingency table of the expected frequencies
# to be used to compute S-X2 fit statistic
expFreq <- function(t.score, cats, prob.cats, lkhd_noitem, lkhd, wts, score.freq) {
  # number of possible scores for the test without this item
  t_noitem <- t.score - cats + 2L

  # compute joint[k, j] = sum_theta w_theta * P(cat=j-1|theta) * P(rest_score=k-1|theta)
  # one BLAS DGEMM computes all K column sums at once
  joint <- crossprod(lkhd_noitem, prob.cats * wts[, 2])  # (t_noitem x cats)

  # fill the (t.score+1) x cats staircase matrix in one vectorized assignment:
  # category j occupies rows j:(j + t_noitem - 1) of column j
  row_idx <- rep(seq_len(t_noitem), cats) + rep(seq_len(cats) - 1L, each = t_noitem)
  col_idx <- rep(seq_len(cats), each = t_noitem)
  tmp1 <- array(0, c(t.score + 1, cats))
  tmp1[cbind(row_idx, col_idx)] <- joint

  # divide each row by the marginal score distribution (one BLAS DGEMV);
  # R recycles denom column-wise so each row i is divided by denom[i]
  denom <- as.vector(crossprod(lkhd, wts[, 2]))  # length t.score+1
  tmp2  <- tmp1 / denom

  colnames(tmp2) <- paste0("score.", 0:(cats - 1))
  rownames(tmp2) <- paste0("score.", 0:t.score)

  # convert the conditional proportions to expected frequencies
  tmp2 <- score.freq * tmp2

  # drop the extreme summed scores and merge the sparse score groups at both ends
  merge_end_scores(tmp2, cats)
}

# This function returns a contingency table of the observed frequencies
# to be used to compute S-X2 fit statistic
obsFreq <- function(rawscore, response, t.score, cats) {
  # Encode (rawscore, response) as a row-major linear index into the
  # (t.score+1) x cats matrix, then rebuild with a single tabulate() call;
  # avoids factor allocation and table() overhead called J times
  lin_idx <- rawscore * cats + response + 1L
  tmp2 <- matrix(tabulate(lin_idx, nbins = (t.score + 1L) * cats),
                 nrow = t.score + 1L, ncol = cats, byrow = TRUE)

  colnames(tmp2) <- paste0("score.", 0:(cats - 1))
  rownames(tmp2) <- paste0("score.", 0:t.score)

  # drop the extreme summed scores and merge the sparse score groups at both ends
  merge_end_scores(tmp2, cats)
}

# This function drops the lowest and highest summed scores and pools the first and
# last (cats - 1) remaining summed scores without counting any score twice
merge_end_scores <- function(tab, cats) {
  # remove the lowest and highest summed scores and keep the matrix form
  tab <- tab[-c(1, nrow(tab)), , drop = FALSE]
  n <- nrow(tab)

  # find the rows pooled into the first score group
  first.idx <- seq_len(min(cats - 1, n))

  # find the rows pooled into the last score group, excluding rows of the first group
  last.idx <- setdiff(seq(max(1, n - cats + 2), n), first.idx)

  # find the rows kept as separate score groups
  mid.idx <- setdiff(seq_len(n), c(first.idx, last.idx))

  # sum the frequencies within the first group
  row.first <- colSums(tab[first.idx, , drop = FALSE])

  # combine the first group with the middle groups
  out <- rbind(row.first, tab[mid.idx, , drop = FALSE])
  rownames(out) <- c(rownames(tab)[max(first.idx)], rownames(tab)[mid.idx])

  # add the last group when it holds any row
  if (length(last.idx) > 0L) {
    row.end <- colSums(tab[last.idx, , drop = FALSE])
    out <- rbind(out, row.end)
    rownames(out)[nrow(out)] <- rownames(tab)[min(last.idx)]
  }

  data.frame(out)
}
