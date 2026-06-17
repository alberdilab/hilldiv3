#' Branch-abundance kernel for phylogenetic Hill numbers
#'
#' Internal helper that, for a phylogenetic tree, computes for every branch
#' (edge) the summed relative abundance of the tips descending from it. This is
#' the `a_i` quantity in the phylogenetic Hill-number framework of Chao et al.
#' (2010). It uses a single post-order traversal in O(edges), replacing the
#' O(edges x tips) `geiger::tips()` pattern used in hilldiv2.
#'
#' @param tree A phylogenetic tree of class `phylo`.
#' @param p A numeric matrix of *relative* abundances (taxa in rows, samples in
#'   columns) whose rownames match `tree$tip.label`. Columns are assumed to be
#'   already normalised (see [tss()]).
#'
#' @return A list with:
#'   \describe{
#'     \item{Li}{Numeric vector of branch lengths, one per edge.}
#'     \item{ai}{Numeric matrix (edges x samples) of descendant-tip abundance
#'       sums for each branch.}
#'   }
#'
#' @details The accumulation works on a samples-by-nodes layout so that each
#'   node's per-sample vector is contiguous in R's column-major storage; the
#'   per-edge update is then a unit-stride column add, which is markedly faster
#'   than the row-wise (strided) access a nodes-by-samples layout would force.
#' @keywords internal
#' @noRd
branch_abundance <- function(tree, p) {
  stopifnot(inherits(tree, "phylo"))
  p <- as.matrix(p)
  # Align tips to the matrix rows.
  p <- p[tree$tip.label, , drop = FALSE]

  n_tip <- length(tree$tip.label)
  n_total <- n_tip + tree$Nnode

  # node_abund[sample, node] = summed abundance of tips below `node`. The layout
  # is samples-by-nodes (each node's per-sample vector contiguous) so the
  # accumulation below is a unit-stride column add. Tip columns start at the
  # observed abundances; internal-node columns start at zero.
  node_abund <- matrix(0, nrow = ncol(p), ncol = n_total)
  node_abund[, seq_len(n_tip)] <- t(p)

  # Post-order: children before parents. ape stores edges parent->child; ape's
  # postorder reorder lists them so each child is summed before its parent is
  # needed.
  tree_po <- ape::reorder.phylo(tree, order = "postorder")
  parent <- tree_po$edge[, 1]
  child <- tree_po$edge[, 2]
  for (e in seq_along(parent)) {
    node_abund[, parent[e]] <- node_abund[, parent[e]] + node_abund[, child[e]]
  }

  # a_i for each edge = abundance below its child node; transpose back to the
  # edges-by-samples shape the engines expect (original edge order).
  ai <- t(node_abund[, tree$edge[, 2], drop = FALSE])
  list(Li = tree$edge.length, ai = ai)
}
