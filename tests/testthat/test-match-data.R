counts <- matrix(
  c(10, 0, 5, 2,
    8, 1, 3, 4),
  nrow = 4,
  dimnames = list(c("t1", "t2", "t3", "t4"), c("s1", "s2"))
)
# Tree shares t1, t2, t3 with the counts and adds t5; t4 is counts-only.
tree <- ape::read.tree(text = "((t1:1,t3:1):1,(t2:1,t5:2):1.5);")

test_that("match_data subsets and reorders counts to the tree tips", {
  md <- suppressMessages(match_data(counts, tree = tree))
  shared <- intersect(tree$tip.label, rownames(counts))
  expect_equal(rownames(md), shared)
  expect_equal(md, counts[shared, , drop = FALSE])
  # the result is now usable with hilldiv once the tree is pruned to match
  pruned <- ape::drop.tip(tree, setdiff(tree$tip.label, rownames(md)))
  expect_silent(suppressMessages(hilldiv(md, q = 0, tree = pruned)))
})

test_that("match_data reports drops from both sides", {
  expect_message(match_data(counts, tree = tree), "Dropped 1 taxon")
  expect_message(match_data(counts, tree = tree), "no counts")
})

test_that("match_data works against a distance matrix", {
  dist <- matrix(0, 3, 3,
                 dimnames = list(c("t3", "t1", "t9"), c("t3", "t1", "t9")))
  md <- suppressMessages(match_data(counts, dist = dist))
  expect_equal(rownames(md), c("t3", "t1"))
})

test_that("match_data errors on bad inputs", {
  expect_error(match_data(counts))
  expect_error(match_data(counts, tree = tree, dist = as.matrix(dist(counts))))
  no_names <- matrix(1:4, 2)
  expect_error(suppressMessages(match_data(no_names, tree = tree)))
  disjoint <- ape::read.tree(text = "(z1:1,z2:1);")
  expect_error(suppressMessages(match_data(counts, tree = disjoint)))
})
