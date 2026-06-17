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
#' @keywords internal
#' @noRd
branch_abundance <- function(tree, p) {
  stopifnot(inherits(tree, "phylo"))
  p <- as.matrix(p)
  # Align tips to the matrix rows.
  p <- p[tree$tip.label, , drop = FALSE]

  n_tip <- length(tree$tip.label)
  n_node <- tree$Nnode
  n_total <- n_tip + n_node
  n_samp <- ncol(p)

  # node_abund[node, sample] = summed abundance of tips below `node`.
  node_abund <- matrix(0, nrow = n_total, ncol = n_samp)
  node_abund[seq_len(n_tip), ] <- p

  # Post-order: children before parents. ape stores edges parent->child; we
  # accumulate child contributions into parents from tips upward.
  # Process nodes in decreasing order of depth via edge ordering.
  edge <- tree$edge
  # Traverse edges in reverse of a root-to-tip order so children are summed
  # before their parent is needed. ape's postorder reorder gives this directly.
  tree_po <- ape::reorder.phylo(tree, order = "postorder")
  for (e in seq_len(nrow(tree_po$edge))) {
    parent <- tree_po$edge[e, 1]
    child <- tree_po$edge[e, 2]
    node_abund[parent, ] <- node_abund[parent, ] + node_abund[child, ]
  }

  # a_i for each edge = abundance below its child node.
  ai <- node_abund[edge[, 2], , drop = FALSE]
  list(Li = tree$edge.length, ai = ai)
}
