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

  expect_equal(unname(got[, "q0"]), unname(vegan::specnumber(comm)))
  expect_equal(unname(got[, "q1"]),
               unname(exp(vegan::diversity(comm, "shannon"))))
  expect_equal(unname(got[, "q2"]),
               unname(vegan::diversity(comm, "invsimpson")))
})

test_that("a single taxon has diversity 1 at every order", {
  one <- matrix(c(5, 9), nrow = 1, dimnames = list("t1", c("s1", "s2")))
  got <- suppressMessages(hilldiv(one, q = c(0, 1, 2), out = "matrix"))
  expect_equal(unname(got), matrix(1, 2, 3))
})

test_that("an all-zero sample yields zero diversity", {
  counts <- matrix(c(0, 0, 0, 4, 2, 6), nrow = 3,
                   dimnames = list(c("t1", "t2", "t3"), c("empty", "s2")))
  got <- suppressMessages(hilldiv(counts, q = c(0, 1, 2), out = "matrix"))
  expect_equal(unname(got["empty", ]), c(0, 0, 0))
  expect_true(all(got["s2", ] > 0))
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

  # Output is a bare matrix for one type and a named list of matrices when a
  # tree or dist adds further types; strip dimnames through either shape.
  strip <- function(z) if (is.list(z)) lapply(z, unname) else unname(z)
  for (args in list(list(), list(tree = tree), list(dist = dist))) {
    lim <- do.call(hilldiv, c(list(counts, q = 1, out = "matrix"), args)) |>
      suppressMessages()
    near <- do.call(hilldiv,
                    c(list(counts, q = 1 + 1e-6, out = "matrix"), args)) |>
      suppressMessages()
    expect_equal(strip(lim), strip(near), tolerance = 1e-4)
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

test_that("identical assemblages give beta = 1 (Jost 2007)", {
  # Two identical assemblages share all their diversity, so the multiplicative
  # beta is 1 at every order regardless of the within-sample distribution.
  counts <- matrix(c(5, 3, 2, 5, 3, 2), nrow = 3,
                   dimnames = list(c("t1", "t2", "t3"), c("s1", "s2")))
  part <- suppressMessages(hillpart(counts, q = c(0, 1, 2), out = "matrix"))
  expect_equal(unname(part[, "beta"]), c(1, 1, 1))
})

test_that("completely distinct equal-weight assemblages give beta = N (Jost 2007)", {
  # N assemblages with no shared taxa and matched within-sample evenness have
  # beta = N at every order -- Jost's (2007) result that beta is independent of
  # alpha and ranges in [1, N].
  counts <- matrix(c(5, 5, 0, 0, 0, 0, 5, 5), nrow = 4,
                   dimnames = list(c("t1", "t2", "t3", "t4"), c("s1", "s2")))
  N <- ncol(counts)
  part <- suppressMessages(hillpart(counts, q = c(0, 1, 2), out = "matrix"))
  expect_equal(unname(part[, "beta"]), rep(N, 3))
})

test_that("q0 turnover equals classic Sorensen/Jaccard dissimilarity (Chiu et al. 2014)", {
  # For two presence/absence assemblages the q = 0 Sorensen-type turnover (V)
  # and Jaccard-type turnover (S) reduce to the classic incidence indices
  # (b + c) / (2a + b + c) and (b + c) / (a + b + c).
  pa <- matrix(c(1, 0, 1, 0, 1, 1, 0, 1, 0, 1), nrow = 5, byrow = TRUE,
               dimnames = list(paste0("t", 1:5), c("s1", "s2")))
  a <- sum(pa[, 1] > 0 & pa[, 2] > 0)
  b <- sum(pa[, 1] > 0 & pa[, 2] == 0)
  cc <- sum(pa[, 1] == 0 & pa[, 2] > 0)
  diss <- suppressMessages(hilldiss(pa, q = c(0, 1, 2), out = "matrix"))
  expect_equal(unname(diss["q0", "V"]), (b + cc) / (2 * a + b + cc)) # Sorensen
  expect_equal(unname(diss["q0", "S"]), (b + cc) / (a + b + cc))     # Jaccard
})

test_that("similarity is the complement of dissimilarity (Chiu et al. 2014)", {
  pa <- matrix(c(1, 0, 1, 0, 1, 1, 0, 1, 0, 1), nrow = 5, byrow = TRUE,
               dimnames = list(paste0("t", 1:5), c("s1", "s2")))
  diss <- suppressMessages(hilldiss(pa, q = c(0, 1, 2), out = "matrix"))
  sim <- suppressMessages(hillsim(pa, q = c(0, 1, 2), out = "matrix"))
  expect_equal(unname(sim), unname(1 - diss))
})

test_that("phylogenetic gamma prunes empty branches (partial Faith's PD)", {
  # Taxa absent from every sample must not contribute branch length: gamma PD at
  # q0 is Faith's PD of the occupied subtree, not the whole tree.
  tree <- ape::read.tree(text = "((t1:1,t2:1):1,(t3:1,t4:2):1.5);")
  counts <- matrix(c(3, 1, 2, 4, 0, 0, 0, 0), nrow = 4, byrow = TRUE,
                   dimnames = list(c("t1", "t2", "t3", "t4"), c("s1", "s2")))
  part <- suppressMessages(hillpart(counts, q = 0, tree = tree, out = "matrix"))
  # Occupied branches: t1 (1), t2 (1) and their shared internal edge (1) = 3;
  # the (t3, t4) clade contributes nothing.
  expect_equal(unname(part["q0", "gamma"]), 3)
  expect_lt(unname(part["q0", "gamma"]), sum(tree$edge.length))
})

test_that("functional q2 equals Rao's quadratic entropy form (Chiu & Chao 2014)", {
  # With tau = max(d), the order-2 functional Hill number has the closed form
  # 1 / (1 - Q / tau), where Q = sum_ij d_ij p_i p_j is Rao's quadratic entropy.
  counts <- matrix(c(5, 3, 2), nrow = 3,
                   dimnames = list(c("t1", "t2", "t3"), "s1"))
  dist <- matrix(c(0, 1, 0.5, 1, 0, 0.8, 0.5, 0.8, 0), nrow = 3,
                 dimnames = list(c("t1", "t2", "t3"), c("t1", "t2", "t3")))
  p <- as.vector(counts) / sum(counts)
  rao <- as.numeric(t(p) %*% dist %*% p)
  tau <- max(dist)
  got <- suppressMessages(
    hilldiv(counts, dist = dist, type = "functional", q = 2, out = "matrix"))
  expect_equal(unname(got[, "q2"]), 1 / (1 - rao / tau))
})

test_that("functional diversity with maximally distinct taxa reduces to neutral", {
  # When every pair is maximally distinct (d = tau), the similarity matrix is
  # the identity and functional Hill numbers collapse to the neutral ones.
  counts <- matrix(c(5, 0, 3, 2, 8, 1), nrow = 3,
                   dimnames = list(c("t1", "t2", "t3"), c("s1", "s2")))
  dist <- matrix(1, 3, 3)
  diag(dist) <- 0
  dimnames(dist) <- list(c("t1", "t2", "t3"), c("t1", "t2", "t3"))
  func <- suppressMessages(
    hilldiv(counts, dist = dist, type = "functional", q = c(0, 1, 2),
            out = "matrix"))
  neutral <- suppressMessages(
    hilldiv(counts, type = "neutral", q = c(0, 1, 2), out = "matrix"))
  expect_equal(func, neutral)
})
