q_set <- c(0, 0.5, 1, 2)

# A fixed, explicit tree so expectations are deterministic.
tree <- ape::read.tree(text = "((t1:1,t2:1):1,(t3:1,t4:2):1.5);")
counts <- matrix(
  c(8, 0, 4, 2,
    0, 5, 3, 1,
    2, 2, 0, 6),
  nrow = 4,
  dimnames = list(c("t1", "t2", "t3", "t4"), c("s1", "s2", "s3"))
)

test_that("phylogenetic partition on a unit-length star tree equals neutral", {
  # On a star phylogeny with unit branch lengths every tip is its own branch,
  # so phylogenetic Hill partitioning collapses onto the neutral case.
  star <- ape::stree(4, type = "star")
  star$tip.label <- rownames(counts)
  star$edge.length <- rep(1, nrow(star$edge))

  phylo <- suppressMessages(hillpart(counts, q = q_set, tree = star,
                                     out = "matrix"))
  neutral <- suppressMessages(hillpart(counts, q = q_set, out = "matrix"))
  expect_equal(phylo, neutral)
})

test_that("phylogenetic partition matches Chiu et al. (2014) structure", {
  part <- suppressMessages(hillpart(counts, q = q_set, tree = tree,
                                    out = "matrix"))
  N <- ncol(counts)

  # beta lies in [1, N]; gamma >= alpha.
  expect_true(all(part[, "beta"] >= 1 - 1e-9))
  expect_true(all(part[, "beta"] <= N + 1e-9))
  expect_true(all(part[, "gamma"] >= part[, "alpha"] - 1e-9))

  # At q = 0, gamma is Faith's PD of the pooled assemblage: the total length of
  # all branches with any descendant abundance (here every branch is present).
  expect_equal(unname(part["q0", "gamma"]), sum(tree$edge.length))
})

test_that("identical samples give no phylogenetic turnover (beta = 1)", {
  one <- counts[, 1, drop = FALSE]
  twin <- cbind(s1 = one[, 1], s2 = one[, 1])
  rownames(twin) <- rownames(counts)
  part <- suppressMessages(hillpart(twin, q = q_set, tree = tree,
                                    out = "matrix"))
  expect_equal(unname(part[, "beta"]), rep(1, length(q_set)))
  expect_equal(unname(part[, "alpha"]), unname(part[, "gamma"]))
})

test_that("q = 1 phylogenetic partition is the limit of nearby q", {
  lim <- suppressMessages(hillpart(counts, q = 1, tree = tree, out = "matrix"))
  near <- suppressMessages(hillpart(counts, q = 1 + 1e-6, tree = tree,
                                    out = "matrix"))
  expect_equal(unname(lim), unname(near), tolerance = 1e-4)
})

test_that("hilldiss/hillsim wire through the phylogenetic path", {
  d <- suppressMessages(hilldiss(counts, q = q_set, tree = tree,
                                 out = "matrix"))
  s <- suppressMessages(hillsim(counts, q = q_set, tree = tree,
                                out = "matrix"))
  expect_equal(d + s, matrix(1, nrow(d), ncol(d), dimnames = dimnames(d)))
  expect_true(all(d >= -1e-9 & d <= 1 + 1e-9))

  # beta = 1 (identical samples) implies zero dissimilarity on every metric.
  twin <- cbind(counts[, 1], counts[, 1])
  rownames(twin) <- rownames(counts)
  d0 <- suppressMessages(hilldiss(twin, q = q_set, tree = tree,
                                  out = "matrix"))
  expect_equal(unname(d0), matrix(0, length(q_set), 4))
})

test_that("phylogenetic partitioning requires a tree", {
  expect_error(hilldiv3:::part_prep(counts, "phylogenetic", tree = NULL))
})
