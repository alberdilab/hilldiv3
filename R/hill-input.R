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
  phy <- phyloseq::phy_tree(data, errorIfNULL = FALSE)
  if (is.null(tree) && !is.null(phy)) {
    tree <- phy
  }
  meta <- tryCatch(methods::as(phyloseq::sample_data(data), "data.frame"),
                   error = function(e) NULL)
  new_hill_input(otu, tree = tree, dist = dist, metadata = meta)
}

#' @noRd
as_hill_input.TreeSummarizedExperiment <- function(data, tree = NULL,
                                                    dist = NULL) {
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

# Resolve the requested diversity type against the supplied inputs. `"auto"`
# (the default) falls back to input-based detection; an explicit type is
# validated -- "phylogenetic"/"functional" require the matching input, and any
# request drops the inputs that type does not use (so e.g. `type = "neutral"`
# deliberately ignores a tree carried along on a phyloseq object). Returns the
# possibly-trimmed `hill_input` and the resolved type.
resolve_type <- function(x, type = "auto") {
  type <- match.arg(type, c("auto", "neutral", "phylogenetic", "functional"))
  detected <- hill_type(x)
  if (type == "auto") {
    return(list(x = x, type = detected))
  }
  if (type == "phylogenetic" && is.null(x$tree)) {
    cli::cli_abort(c(
      "{.code type = \"phylogenetic\"} needs a {.arg tree}.",
      "i" = "Supply a {.cls phylo} tree, or use {.code type = \"auto\"}."
    ))
  }
  if (type == "functional" && is.null(x$dist)) {
    cli::cli_abort(c(
      "{.code type = \"functional\"} needs a {.arg dist} matrix.",
      "i" = "Supply a distance matrix, or use {.code type = \"auto\"}."
    ))
  }
  if (type != "phylogenetic") x$tree <- NULL
  if (type != "functional") x$dist <- NULL
  list(x = x, type = type)
}
