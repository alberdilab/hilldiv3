#' Coerce inputs into a uniform internal representation
#'
#' Adapter layer that turns the various objects users may supply (matrix, data
#' frame, tibble, `phyloseq`, `TreeSummarizedExperiment`) into a single internal
#' `hill_input` structure. Everything downstream (validation, engine) operates
#' on this uniform object.
#'
#' @param data A counts object: numeric matrix/data.frame (taxa x samples),
#'   a numeric vector (single sample), a `phyloseq` object, or a
#'   `TreeSummarizedExperiment`.
#' @param tree Optional `phylo` tree. Ignored if `data` already carries a tree.
#' @param dist Optional distance matrix/`dist` of functional distances.
#'
#' @return An object of class `hill_input`: a list with `counts` (matrix,
#'   taxa x samples), `tree`, `dist` and `metadata` (or `NULL`).
#' @keywords internal
#' @noRd
as_hill_input <- function(data, tree = NULL, dist = NULL) {
  UseMethod("as_hill_input")
}

#' @noRd
as_hill_input.default <- function(data, tree = NULL, dist = NULL) {
  if (is.null(dim(data))) {
    # Single-sample vector -> one-column matrix.
    counts <- matrix(data, ncol = 1, dimnames = list(names(data), "sample1"))
  } else {
    counts <- as.matrix(data)
  }
  new_hill_input(counts, tree = tree, dist = dist)
}

#' @noRd
as_hill_input.phyloseq <- function(data, tree = NULL, dist = NULL) {
  rlang::check_installed("phyloseq")
  otu <- methods::as(phyloseq::otu_table(data), "matrix")
  if (!phyloseq::taxa_are_rows(data)) otu <- t(otu)
  if (is.null(tree) && !is.null(phyloseq::phy_tree(data, errorIfNULL = FALSE))) {
    tree <- phyloseq::phy_tree(data)
  }
  meta <- tryCatch(methods::as(phyloseq::sample_data(data), "data.frame"),
                   error = function(e) NULL)
  new_hill_input(otu, tree = tree, dist = dist, metadata = meta)
}

#' @noRd
as_hill_input.TreeSummarizedExperiment <- function(data, tree = NULL, dist = NULL) {
  rlang::check_installed("TreeSummarizedExperiment")
  counts <- as.matrix(SummarizedExperiment::assay(data))
  if (is.null(tree)) {
    tree <- TreeSummarizedExperiment::rowTree(data)
  }
  meta <- as.data.frame(SummarizedExperiment::colData(data))
  new_hill_input(counts, tree = tree, dist = dist, metadata = meta)
}

# Constructor.
new_hill_input <- function(counts, tree = NULL, dist = NULL, metadata = NULL) {
  if (!is.null(dist) && !is.matrix(dist)) dist <- as.matrix(dist)
  structure(
    list(counts = counts, tree = tree, dist = dist, metadata = metadata),
    class = "hill_input"
  )
}

# Decide the diversity type from which inputs are present.
hill_type <- function(x) {
  has_tree <- !is.null(x$tree)
  has_dist <- !is.null(x$dist)
  if (has_tree && has_dist) {
    cli::cli_abort(c(
      "Cannot use phylogenetic and functional information at once.",
      "i" = "Supply either {.arg tree} or {.arg dist}, not both."
    ))
  }
  if (has_tree) return("phylogenetic")
  if (has_dist) return("functional")
  "neutral"
}
