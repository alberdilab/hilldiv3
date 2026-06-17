#' Hill numbers computation
#'
#' Compute neutral, phylogenetic or functional Hill numbers (alpha diversity)
#' from a single sample or a count table. The diversity type is inferred from
#' the inputs: counts only -> neutral; counts + `tree` -> phylogenetic;
#' counts + `dist` -> functional.
#'
#' @param data Counts: a numeric vector (one sample), a matrix/data.frame
#'   (taxa x samples), a `phyloseq` object or a `TreeSummarizedExperiment`.
#' @param q Numeric vector of diversity orders (>= 0). Defaults to
#'   `c(0, 1, 2)` (richness, Shannon, Simpson).
#' @param tree A phylogenetic tree of class `phylo` whose tip labels match the
#'   taxa in `data`.
#' @param dist A functional distance matrix (or `dist`) over the taxa.
#' @param tau Optional functional distance threshold. Defaults to `max(dist)`.
#' @param type Diversity type: `"auto"` (default) infers it from the inputs
#'   (counts only -> neutral, `+tree` -> phylogenetic, `+dist` -> functional);
#'   an explicit `"neutral"`, `"phylogenetic"` or `"functional"` asserts the
#'   type and is validated against the inputs (e.g. `"phylogenetic"` requires a
#'   `tree`; `"neutral"` ignores any tree/dist carried by the object).
#' @param reference Reference tree depth for *phylogenetic* Hill numbers
#'   (ignored for neutral and functional types). `"pool"` (default) reads every
#'   sample at one common depth `T = mean(T_j)`, so values share a comparable
#'   axis (matching hilldiv2's `multi` behaviour); `"sample"` reads each sample
#'   at its own depth `T_j` (effective lineages at that sample's depth). The two
#'   coincide on ultrametric trees. This reference depth is intentionally *not*
#'   offered by [hillpart()]: in a partition `T` is fixed at the mean per-sample
#'   depth of Chiu et al. (2014), the unique value for which `gamma / alpha`
#'   is a valid decomposition with `beta` in `[1, N]`.
#' @param out Output shape: `"tibble"` (default) returns a long-format
#'   `data.frame` with columns `q`, `sample`, `value` and `print()`/`plot()`
#'   methods; `"matrix"` returns the legacy matrix (orders in rows, samples in
#'   columns).
#'
#' @return A long-format `data.frame` of class `hill_diversity` (default), or a
#'   matrix of Hill numbers with diversity orders in rows (`q0`, `q1`, ...) and
#'   samples in columns when `out = "matrix"`.
#'
#' @references
#' Chao, A., Chiu, C.-H. & Jost, L. (2010). Phylogenetic diversity measures
#' based on Hill numbers. Phil. Trans. R. Soc. B, 365, 3599-3609.\cr\cr
#' Alberdi, A. & Gilbert, M.T.P. (2019). A guide to the application of Hill
#' numbers to DNA-based diversity analyses. Mol. Ecol. Resour., 19, 804-817.
#'
#' @seealso [hillpart()], [hilldiss()], [hillprof()]
#' @examples
#' counts <- matrix(c(10, 0, 5, 2, 8, 1), nrow = 3,
#'                  dimnames = list(c("t1", "t2", "t3"), c("s1", "s2")))
#' hilldiv(counts)
#' hilldiv(counts, q = c(0, 1, 2))
#' plot(hilldiv(counts, q = c(0, 1, 2)))
#' @export
hilldiv <- function(data, q = c(0, 1, 2), tree = NULL, dist = NULL, tau = NULL,
                    type = c("auto", "neutral", "phylogenetic", "functional"),
                    reference = c("pool", "sample"),
                    out = c("tibble", "matrix")) {
  out <- match.arg(out)
  type <- match.arg(type)
  reference <- match.arg(reference)
  x <- prep_data(as_hill_input(data, tree = tree, dist = dist), q, type)
  type <- attr(x, "type")
  cli::cli_inform("Computing {type} Hill numbers of {.val {paste0('q', q)}}.")
  mat <- hill_alpha(x$counts, q = q, type = type,
                    tree = x$tree, dist = x$dist, tau = tau,
                    reference = reference)
  if (out == "matrix") return(mat)
  new_hill_result(.hill_longify(mat, q, "sample"), "hill_diversity", type)
}
