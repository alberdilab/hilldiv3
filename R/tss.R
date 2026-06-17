#' Total Sum Scaling normalisation
#'
#' Normalise a numeric vector or count matrix so that each sample (column) sums
#' to one. Columns that sum to zero are returned as all-zero (the `0/0 = NaN`
#' case is mapped to `0`).
#'
#' @param abund A numeric vector or a matrix/data.frame of counts with taxa
#'   (OTUs/ASVs/MAGs) in rows and samples in columns.
#'
#' @return A normalised object of the same shape as `abund` (vector in, vector
#'   out; matrix/data.frame in, matrix out).
#'
#' @examples
#' tss(c(a = 1, b = 3))
#' tss(matrix(c(1, 0, 3, 0, 0, 2), nrow = 3))
#' @export
tss <- function(abund) {
  if (is.null(dim(abund))) {
    total <- sum(abund)
    if (total == 0) return(abund * 0)
    return(abund / total)
  }
  m <- as.matrix(abund)
  sweep_totals <- colSums(m)
  out <- sweep(m, 2, sweep_totals, FUN = "/")
  out[is.nan(out)] <- 0
  out
}
