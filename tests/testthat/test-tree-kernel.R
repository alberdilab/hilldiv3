# branch_abundance() must return, for every edge, the summed abundance of the
# tips descending from that edge's child node (the a_i of Chao et al. 2010),
# in the tree's original edge order and shaped edges x samples.

test_that("branch abundances match hand-computed descendant sums", {
  # ((t1,t2),(t3,t4)); edges in ape's storage order.
  tree <- ape::read.tree(text = "((t1:1,t2:1):1,(t3:1,t4:2):1.5);")
  p <- matrix(c(0.4, 0.3, 0.2, 0.1,
                0.1, 0.1, 0.5, 0.3),
              nrow = 4,
              dimnames = list(c("t1", "t2", "t3", "t4"), c("s1", "s2")))

  ba <- branch_abundance(tree, p)

  expect_identical(ba$Li, tree$edge.length)
  expect_equal(dim(ba$ai), c(nrow(tree$edge), ncol(p)))

  # Reference: for each edge, sum the abundances of the descendant tips.
  tips_below <- function(node) {
    if (node <= length(tree$tip.label)) return(tree$tip.label[node])
    kids <- tree$edge[tree$edge[, 1] == node, 2]
    unlist(lapply(kids, tips_below))
  }
  ref <- t(vapply(tree$edge[, 2], function(node) {
    colSums(p[tips_below(node), , drop = FALSE])
  }, numeric(ncol(p))))
  expect_equal(unname(ba$ai), unname(ref))
})

test_that("root edges carry the full per-sample abundance", {
  # A node feeding the two basal clades sees every tip, so summed a_i over the
  # two root-child edges equals the column totals.
  tree <- ape::rtree(30)
  p <- matrix(runif(30 * 3), nrow = 30,
              dimnames = list(tree$tip.label, c("s1", "s2", "s3")))
  ba <- branch_abundance(tree, p)

  root <- length(tree$tip.label) + 1L
  root_edges <- which(tree$edge[, 1] == root)
  expect_equal(unname(colSums(ba$ai[root_edges, , drop = FALSE])),
               unname(colSums(p)))
})
