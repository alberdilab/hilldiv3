#' Validate and align a hill_input
#'
#' Single validation/alignment layer. Checks that q values are valid, that taxa
#' names agree between the counts and the tree/distance matrix using
#' [setequal()] (not the buggy elementwise comparison used in hilldiv2), and
#' **reorders** the counts so their row order matches the tree tips / distance
#' rows. This is what guarantees downstream alignment.
#'
#' @param x A `hill_input` (see [as_hill_input()]).
#' @param q Numeric vector of diversity orders.
#' @param type Requested diversity type, one of `"auto"` (detect from inputs),
#'   `"neutral"`, `"phylogenetic"` or `"functional"`. Resolved and validated by
#'   [resolve_type()].
#'
#' @return The validated, aligned `hill_input`, with a `type` attribute.
#' @keywords internal
#' @noRd
prep_data <- function(x, q, type = "auto") {
  if (any(q < 0)) {
    cli::cli_abort("Diversity orders {.arg q} must be >= 0.")
  }
  if (!is.numeric(x$counts)) {
    cli::cli_abort("Count data must be numeric.")
  }

  resolved <- resolve_type(x, type)
  x <- resolved$x
  type <- resolved$type
  taxa <- rownames(x$counts)

  if (type == "phylogenetic") {
    if (is.null(taxa)) {
      cli::cli_abort("Count data must have row names to match the tree tips.")
    }
    if (!setequal(taxa, x$tree$tip.label)) {
      .abort_mismatch(taxa, x$tree$tip.label, "tree tips")
    }
    x$counts <- x$counts[x$tree$tip.label, , drop = FALSE]
  }

  if (type == "functional") {
    dnames <- rownames(x$dist)
    if (is.null(taxa) || is.null(dnames)) {
      cli::cli_abort("Both count data and the distance matrix need names.")
    }
    if (!setequal(taxa, dnames)) {
      .abort_mismatch(taxa, dnames, "distance matrix")
    }
    x$counts <- x$counts[dnames, , drop = FALSE]
    x$dist <- x$dist[dnames, dnames, drop = FALSE]
  }

  attr(x, "type") <- type
  x
}

# Build a helpful mismatch error message.
.abort_mismatch <- function(a, b, what) {
  only_a <- setdiff(a, b)
  only_b <- setdiff(b, a)
  cli::cli_abort(c(
    "Taxa names in the count data and the {what} do not match.",
    "i" = if (length(only_a)) "Only in counts: {.val {utils::head(only_a, 5)}}",
    "i" = if (length(only_b)) "Only in {what}: {.val {utils::head(only_b, 5)}}"
  ))
}

#' Match and align a count table to a tree or distance matrix
#'
#' Subsets and reorders a count table so that its taxa match those of a
#' phylogenetic tree or a functional distance matrix, dropping taxa absent from
#' the reference. This realises the `match_data()` helper that hilldiv2's
#' documentation referred to but never provided.
#'
#' @param data A count matrix/data.frame (taxa x samples) with row names.
#' @param tree A `phylo` tree (optional).
#' @param dist A distance matrix (optional).
#'
#' @return The count matrix restricted to and ordered by the shared taxa.
#' @examples
#' counts <- matrix(1:6, nrow = 3,
#'                  dimnames = list(c("t1", "t2", "t3"), c("s1", "s2")))
#' tree <- ape::read.tree(text = "((t1:1,t2:1):1,t4:2);")
#' match_data(counts, tree = tree)
#' @export
match_data <- function(data, tree = NULL, dist = NULL) {
  if (!is.null(tree) && !is.null(dist)) {
    cli::cli_abort("Supply either {.arg tree} or {.arg dist}, not both.")
  }
  if (is.null(tree) && is.null(dist)) {
    cli::cli_abort("Provide a {.arg tree} or {.arg dist} to match {.arg data}
                    against.")
  }

  counts <- as_hill_input(data)$counts
  taxa <- rownames(counts)
  if (is.null(taxa)) {
    cli::cli_abort("{.arg data} must have row names (taxa) to match.")
  }

  if (!is.null(tree)) {
    ref <- tree$tip.label
    what <- "tree tips"
  } else {
    ref <- rownames(as.matrix(dist))
    what <- "distance matrix"
    if (is.null(ref)) {
      cli::cli_abort("{.arg dist} must have names to match against.")
    }
  }

  shared <- intersect(ref, taxa)            # keep the reference ordering
  if (length(shared) == 0) {
    cli::cli_abort("No taxa in common between {.arg data} and the {what}.")
  }

  n_drop_data <- length(setdiff(taxa, ref))
  n_drop_ref <- length(setdiff(ref, taxa))
  if (n_drop_data > 0) {
    cli::cli_inform("Dropped {n_drop_data} taxon{?s} from {.arg data} not in the
                     {what}.")
  }
  if (n_drop_ref > 0) {
    cli::cli_inform("{n_drop_ref} taxon{?s} in the {what} {?has/have} no counts;
                     prune {?it/them} before downstream analysis.")
  }

  counts[shared, , drop = FALSE]
}
