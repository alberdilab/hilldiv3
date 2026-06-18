# Non-ultrametric so that pool (common T) and sample (per-sample T) differ.
phylo_tree <- function() {
  ape::read.tree(text = "((t1:1,t2:1):1,(t3:1,t4:2):1.5);")
}

phylo_counts <- function() {
  matrix(c(5, 3, 1, 1, 2, 2, 6, 2), ncol = 2,
         dimnames = list(c("t1", "t2", "t3", "t4"), c("s1", "s2")))
}

test_that("q = 0 phylogenetic Hill number is Faith's PD / T", {
  tree <- phylo_tree()
  counts <- matrix(c(5, 3, 1, 1), ncol = 1,
                   dimnames = list(c("t1", "t2", "t3", "t4"), "s1"))
  out <- suppressMessages(
    hilldiv(counts, q = 0, tree = tree, reference = "sample", out = "matrix"))
  # All tips present -> PD is the whole tree length; T is the abundance-weighted
  # mean depth, here equal to tree height since the tree is rooted at depth.
  ba <- hilldiv3:::branch_abundance(tree, tss(counts))
  Tj <- sum(ba$Li * ba$ai[, 1])
  expect_equal(unname(out[1, "q0"]), sum(tree$edge.length) / Tj)
})

test_that("sample reference agrees with the partition engine at N = 1", {
  tree <- phylo_tree()
  counts <- phylo_counts()
  q <- c(0, 1, 2)
  for (j in seq_len(ncol(counts))) {
    cj <- counts[, j, drop = FALSE]
    prep <- hilldiv3:::part_prep(cj, "phylogenetic", tree = tree)
    pe <- hilldiv3:::part_eval(prep, 1, q)
    Tj <- sum(prep$Li * prep$aij[, 1])
    hd <- suppressMessages(
      hilldiv(cj, q = q, tree = tree, reference = "sample", out = "matrix"))
    expect_equal(unname(hd[1, ]), unname(pe[, "alpha"] / Tj))
  }
})

test_that("pool reference reads all samples at one common depth", {
  tree <- phylo_tree()
  counts <- phylo_counts()
  out <- suppressMessages(
    hilldiv(counts, q = 0, tree = tree, reference = "pool", out = "matrix"))
  # Both samples contain every tip, so at q = 0 the pooled-depth phylo Hill
  # number (PD / T_pool) is identical across samples.
  expect_equal(unname(out["s1", "q0"]), unname(out["s2", "q0"]))
})

test_that("pool and sample references coincide on an ultrametric tree", {
  ut <- ape::read.tree(text = "((t1:2,t2:2):1,(t3:1.5,t4:1.5):1.5);")
  counts <- phylo_counts()
  q <- c(0, 1, 2)
  a <- suppressMessages(
    hilldiv(counts, q = q, tree = ut, reference = "sample", out = "matrix"))
  b <- suppressMessages(
    hilldiv(counts, q = q, tree = ut, reference = "pool", out = "matrix"))
  expect_equal(a, b)
})

test_that("pool is the default reference", {
  tree <- phylo_tree()
  counts <- phylo_counts()
  q <- c(0, 1, 2)
  default <- suppressMessages(
    hilldiv(counts, q = q, tree = tree, out = "matrix"))
  pool <- suppressMessages(
    hilldiv(counts, q = q, tree = tree, reference = "pool", out = "matrix"))
  expect_equal(default, pool)
})
