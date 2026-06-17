counts <- matrix(
  c(10, 0, 5, 2,
    8, 1, 3, 4,
    0, 6, 2, 7),
  nrow = 4,
  dimnames = list(c("t1", "t2", "t3", "t4"), c("s1", "s2", "s3"))
)
tree <- ape::read.tree(text = "((t1:1,t2:1):1,(t3:1,t4:2):1.5);")
dist <- as.matrix(stats::dist(matrix(c(1, 2, 5, 6, 0, 1, 4, 5), nrow = 4,
                                     dimnames = list(rownames(counts), NULL))))

test_that("type = 'auto' reproduces input-based detection", {
  auto <- suppressMessages(hilldiv(counts, q = c(0, 1), tree = tree,
                                   type = "auto", out = "matrix"))
  detected <- suppressMessages(hilldiv(counts, q = c(0, 1), tree = tree,
                                       out = "matrix"))
  expect_equal(auto, detected)
})

test_that("explicit type matching the inputs equals auto", {
  ph <- suppressMessages(hilldiv(counts, q = c(0, 1), tree = tree,
                                 type = "phylogenetic", out = "matrix"))
  auto <- suppressMessages(hilldiv(counts, q = c(0, 1), tree = tree,
                                   out = "matrix"))
  expect_equal(ph, auto)
})

test_that("type = 'neutral' ignores a supplied tree", {
  forced <- suppressMessages(hilldiv(counts, q = c(0, 1), tree = tree,
                                     type = "neutral", out = "matrix"))
  plain <- suppressMessages(hilldiv(counts, q = c(0, 1), out = "matrix"))
  expect_equal(forced, plain)
})

test_that("requesting phylogenetic without a tree errors", {
  expect_error(
    suppressMessages(hilldiv(counts, q = 0, type = "phylogenetic")),
    "needs a"
  )
})

test_that("requesting functional without a dist errors", {
  expect_error(
    suppressMessages(hilldiv(counts, q = 0, type = "functional")),
    "needs a"
  )
})

test_that("type threads through partitioning and dissimilarity", {
  p <- suppressMessages(hillpart(counts, q = c(0, 1), tree = tree,
                                 type = "phylogenetic", out = "matrix"))
  expect_equal(colnames(p), c("alpha", "gamma", "beta"))

  d <- suppressMessages(hilldiss(counts, q = c(0, 1), dist = dist,
                                 type = "functional", out = "matrix"))
  expect_true(all(d >= -1e-9 & d <= 1 + 1e-9))
})
