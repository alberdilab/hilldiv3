# Golden-value tests: cross-check the engine against independent references
# (vegan, hand-computed constants) and pin down edge-case behaviour.

test_that("neutral Hill numbers match vegan", {
  skip_if_not_installed("vegan")
  counts <- t(matrix(c(10, 0, 5, 2, 8, 1, 3, 4, 0, 6, 2, 7), nrow = 3,
                     dimnames = list(c("t1", "t2", "t3"),
                                     c("s1", "s2", "s3", "s4"))))
  # vegan works on samples x species; hilldiv on species x samples.
  comm <- counts
  got <- suppressMessages(hilldiv(t(comm), q = c(0, 1, 2), out = "matrix"))

  expect_equal(unname(got["q0", ]), unname(vegan::specnumber(comm)))
  expect_equal(unname(got["q1", ]),
               unname(exp(vegan::diversity(comm, "shannon"))))
  expect_equal(unname(got["q2", ]),
               unname(vegan::diversity(comm, "invsimpson")))
})

test_that("a single taxon has diversity 1 at every order", {
  one <- matrix(c(5, 9), nrow = 1, dimnames = list("t1", c("s1", "s2")))
  got <- suppressMessages(hilldiv(one, q = c(0, 1, 2), out = "matrix"))
  expect_equal(unname(got), matrix(1, 3, 2))
})

test_that("an all-zero sample yields zero diversity", {
  counts <- matrix(c(0, 0, 0, 4, 2, 6), nrow = 3,
                   dimnames = list(c("t1", "t2", "t3"), c("empty", "s2")))
  got <- suppressMessages(hilldiv(counts, q = c(0, 1, 2), out = "matrix"))
  expect_equal(unname(got[, "empty"]), c(0, 0, 0))
  expect_true(all(got[, "s2"] > 0))
})

test_that("phylogenetic q0 gamma equals Faith's PD of the pooled tree", {
  tree <- ape::read.tree(text = "((t1:1,t2:1):1,(t3:1,t4:2):1.5);")
  counts <- matrix(c(8, 0, 4, 2, 0, 5, 3, 1), nrow = 4,
                   dimnames = list(c("t1", "t2", "t3", "t4"), c("s1", "s2")))
  part <- suppressMessages(hillpart(counts, q = 0, tree = tree,
                                    out = "matrix"))
  # Every branch has descendant abundance, so PD is the whole tree length.
  expect_equal(unname(part["q0", "gamma"]), sum(tree$edge.length))
})

test_that("q = 1 is the continuous limit of nearby q for every type", {
  tree <- ape::read.tree(text = "((t1:1,t2:1):1,(t3:1,t4:2):1.5);")
  dist <- as.matrix(stats::dist(seq_len(4)))
  dimnames(dist) <- list(c("t1", "t2", "t3", "t4"), c("t1", "t2", "t3", "t4"))
  counts <- matrix(c(8, 1, 4, 2, 3, 5, 3, 1), nrow = 4,
                   dimnames = list(c("t1", "t2", "t3", "t4"), c("s1", "s2")))

  for (args in list(list(), list(tree = tree), list(dist = dist))) {
    lim <- do.call(hilldiv, c(list(counts, q = 1, out = "matrix"), args)) |>
      suppressMessages()
    near <- do.call(hilldiv,
                    c(list(counts, q = 1 + 1e-6, out = "matrix"), args)) |>
      suppressMessages()
    expect_equal(unname(lim), unname(near), tolerance = 1e-4)
  }
})

test_that("bundled gut data runs through neutral and phylogenetic types", {
  expect_silent(suppressMessages(hilldiv(gut_counts, q = c(0, 1, 2))))
  expect_silent(suppressMessages(hilldiv(gut_counts, q = 1, tree = gut_tree)))
})

test_that("bundled gut data runs through the functional type", {
  skip_if_not_installed("cluster")
  d <- traits2dist(gut_traits)
  expect_silent(suppressMessages(hilldiv(gut_counts, q = 1, dist = d)))
})
