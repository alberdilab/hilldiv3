#' Hill numbers computation
#'
#' Compute neutral, phylogenetic and/or functional Hill numbers (alpha
#' diversity) from a single sample or a count table. By default the computation
#' is *cumulative*: every diversity type whose inputs are present is returned.
#' Counts are always available, so neutral is always computed; a `tree` adds
#' phylogenetic and a `dist` adds functional. Supplying both a `tree` and a
#' `dist` therefore returns neutral, phylogenetic and functional side by side in
#' a single tibble (with a `type` column). Use `type` to restrict the output to
#' a subset.
#'
#' @param data Counts: a numeric vector (one sample), a matrix/data.frame
#'   (taxa x samples), a `phyloseq` object or a `TreeSummarizedExperiment`.
#' @param q Numeric vector of diversity orders (>= 0). Defaults to
#'   `c(0, 1, 2)` (richness, Shannon, Simpson).
#' @param tree A phylogenetic tree of class `phylo` whose tip labels match the
#'   taxa in `data`.
#' @param dist A functional distance matrix (or `dist`) over the taxa.
#' @param tau Optional functional distance threshold. Defaults to `max(dist)`.
#' @param type Diversity type(s) to compute. `"auto"` (default) returns every
#'   type whose inputs are present (always neutral, plus phylogenetic with a
#'   `tree` and functional with a `dist`). Pass an explicit type, or a character
#'   vector of types, to restrict the output -- e.g. `"neutral"` ignores any
#'   tree/dist carried by the object, and `c("neutral", "phylogenetic")` drops
#'   functional even when a `dist` is supplied. A requested type that lacks its
#'   input (e.g. `"phylogenetic"` without a `tree`) is an error.
#' @param reference Reference tree depth for *phylogenetic* Hill numbers
#'   (ignored for neutral and functional types). `"pool"` (default) reads every
#'   sample at one common depth `T = mean(T_j)`, so values share a comparable
#'   axis across samples; `"sample"` reads each sample at its own depth `T_j`
#'   (effective lineages at that sample's depth). The two
#'   coincide on ultrametric trees. This reference depth is intentionally *not*
#'   offered by [hillpart()]: in a partition `T` is fixed at the mean per-sample
#'   depth of Chiu et al. (2014), the unique value for which `gamma / alpha`
#'   is a valid decomposition with `beta` in `[1, N]`.
#' @param out Output shape: `"tibble"` (default) returns a long-format
#'   `data.frame` with columns `q`, `sample`, `value` (plus a `type` column when
#'   more than one type is computed) and `print()`/`plot()` methods; `"matrix"`
#'   returns a matrix with samples in rows and diversity orders (`q0`, `q1`,
#'   ...) in columns, or, when more than one type is computed, a named list of
#'   such matrices (one per type).
#'
#' @return A long-format `data.frame` of class `hill_diversity` (default). With
#'   `out = "matrix"`, a matrix of Hill numbers (samples in rows, diversity
#'   orders `q0`, `q1`, ... in columns) for a single type, or a named list of
#'   such matrices when several types are computed.
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
#'
#' # Supplying both a tree and a distance matrix returns neutral, phylogenetic
#' # and functional diversity together, distinguished by a `type` column.
#' tree <- ape::read.tree(text = "((t1:1,t2:1):1,t3:2);")
#' dist <- as.matrix(stats::dist(c(t1 = 0, t2 = 1, t3 = 4)))
#' hilldiv(counts, tree = tree, dist = dist)
#'
#' # Restrict the output with `type` (a scalar or a vector):
#' hilldiv(counts, tree = tree, dist = dist, type = c("neutral", "functional"))
#' @export
hilldiv <- function(data, q = c(0, 1, 2), tree = NULL, dist = NULL, tau = NULL,
                    type = c("auto", "neutral", "phylogenetic", "functional"),
                    reference = c("pool", "sample"),
                    out = c("tibble", "matrix")) {
  out <- match.arg(out)
  type <- match.arg(type, several.ok = TRUE)
  reference <- match.arg(reference)

  xin <- as_hill_input(data, tree = tree, dist = dist)
  types <- .hilldiv_types(xin, type)
  n_taxa <- nrow(xin$counts)
  n_samples <- ncol(xin$counts)
  cli::cli_inform(c(
    "Computing {.val {types}} Hill numbers of {.val {paste0('q', q)}}.",
    "i" = "{n_taxa} {?taxon/taxa} across {n_samples} sample{?s}."
  ))

  # Compute each type independently. Strip the irrelevant reference first so the
  # shared single-type prep_data() aligns the counts for this type (and so the
  # "both tree and dist" guard in hill_type() never fires).
  mats <- lapply(types, function(ty) {
    xt <- xin
    if (ty != "phylogenetic") xt$tree <- NULL
    if (ty != "functional") xt$dist <- NULL
    x <- prep_data(xt, q, ty)
    hill_alpha(x$counts, q = q, type = ty, tree = x$tree, dist = x$dist,
               tau = tau, reference = reference)
  })
  names(mats) <- types

  if (out == "matrix") {
    if (length(types) == 1L) return(t(mats[[1L]]))
    return(lapply(mats, t))
  }

  if (length(types) == 1L) {
    df <- .hill_longify(mats[[1L]], q, "sample")
    return(new_hill_result(df, "hill_diversity", types))
  }

  parts <- Map(function(mat, ty) {
    d <- .hill_longify(mat, q, "sample")
    d$type <- ty
    d[c("q", "sample", "type", "value")]
  }, mats, types)
  df <- do.call(rbind, parts)
  rownames(df) <- NULL
  new_hill_result(df, "hill_diversity", paste(types, collapse = ", "))
}

# Which diversity types should hilldiv() compute? Unlike the single-type
# resolve_type() shared by the other hill* functions, hilldiv() is cumulative:
# "auto" returns every type whose inputs are present, and an explicit `type`
# (scalar or vector) selects a subset in the canonical neutral/phylo/functional
# order. Validation that a requested type has its input is deferred to the
# per-type prep_data() call in the loop above, which reuses resolve_type()'s
# error messages.
.hilldiv_types <- function(x, type) {
  if ("auto" %in% type) {
    return(c("neutral",
             if (!is.null(x$tree)) "phylogenetic",
             if (!is.null(x$dist)) "functional"))
  }
  types <- intersect(c("neutral", "phylogenetic", "functional"), type)
  if (length(types) == 0L) {
    cli::cli_abort("No valid diversity {.arg type} requested.")
  }
  types
}
