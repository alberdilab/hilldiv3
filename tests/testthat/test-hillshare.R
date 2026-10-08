test_that("neutral q-power shares handle absent taxa and additive sets", {
  counts <- matrix(c(3, 1, 0, 0, 0, 0), nrow = 3,
                   dimnames = list(c("A", "B", "C"), c("present", "empty")))
  out <- hillshare(counts, q = c(0, 1, 2),
                   sets = list(A = "A", AB = c("A", "B"), C = "C"))

  at <- function(q, set, sample = "present") {
    out$share[out$q == q & out$set == set & out$sample == sample]
  }
  expect_equal(at(0, "A"), 1 / 2)
  expect_equal(at(1, "A"), 3 / 4)
  expect_equal(at(2, "A"), 9 / 10)
  expect_equal(at(2, "AB"), 1)
  expect_equal(at(0, "C"), 0)
  expect_true(all(is.na(out$share[out$sample == "empty"])))
})

test_that("phylogenetic shares allocate the Hill branch power sum", {
  tree <- ape::read.tree(
    text = "(((A:0.300,B:0.300):0.185,C:0.485):0.009,D:0.494);"
  )
  counts <- c(A = 394, B = 202, C = 354, D = 50)
  out <- hillshare(counts, q = c(0, 1, 3), tree = tree,
                   type = "phylogenetic",
                   sets = list(C = "C", ABC = c("A", "B", "C")))
  at <- function(q, set) out$share[out$q == q & out$set == set]

  # At q = 0, C has its own edge and one third of the ABC stem.
  expect_equal(at(0, "C"), (.485 + .009 / 3) / sum(tree$edge.length))
  # On an ultrametric tree, each tip's q = 1 allocation equals abundance.
  expect_equal(at(1, "C"), .354)
  expect_equal(at(1, "ABC"), .950)
  # Worked example from the documented branch ledger.
  expect_equal(at(3, "C"), .27468656308105, tolerance = 1e-12)
  expect_equal(at(3, "ABC"), .9993083663, tolerance = 1e-9)
})

test_that("phylogenetic shares on a star reduce to neutral shares", {
  tree <- ape::read.tree(text = "(A:1,B:1,C:1);")
  counts <- c(A = 5, B = 2, C = 0)
  out <- hillshare(counts, q = c(0, .5, 1, 2, 3),
                   tree = tree, sets = list(A = "A", B = "B", C = "C"))
  neutral <- out[out$type == "neutral", "share"]
  phylo <- out[out$type == "phylogenetic", "share"]
  expect_equal(phylo, neutral)
})

test_that("functional shares use the summands of the existing Hill power sum", {
  counts <- c(A = 5, B = 3, C = 2)
  d <- matrix(c(0, .5, 1, .5, 0, .5, 1, .5, 0), 3,
              dimnames = list(names(counts), names(counts)))
  out <- hillshare(counts, q = c(0, 1, 2), dist = d,
                   type = "functional")
  # The similarities give ordinariness (0.65, 0.65, 0.35).
  expect_equal(out$share[out$q == 0],
               c(.426829268293, .256097560976, .317073170732),
               tolerance = 1e-11)
  expect_equal(out$share[out$q == 1], c(.5, .3, .2))
  expect_equal(out$share[out$q == 2],
               c(.550847457627, .330508474576, .118644067797),
               tolerance = 1e-11)

  # A set is the sum of the taxon shares.
  grouped <- hillshare(counts, q = 2, dist = d, type = "functional",
                       sets = list(AB = c("A", "B")))
  expect_equal(grouped$share, sum(out$share[out$q == 2 & out$set %in% c("A", "B")]))
})

test_that("functional identity similarity reduces to neutral q-power shares", {
  counts <- c(A = 5, B = 3, C = 2)
  d <- matrix(1, 3, 3, dimnames = list(names(counts), names(counts)))
  diag(d) <- 0
  out <- hillshare(counts, q = c(0, .5, 1, 2), dist = d)
  expect_equal(out$share[out$type == "functional"],
               out$share[out$type == "neutral"])
})

test_that("fully identical functional taxa retain abundance shares", {
  counts <- c(A = 5, B = 3, C = 2)
  d <- matrix(0, 3, 3, dimnames = list(names(counts), names(counts)))
  out <- hillshare(counts, q = c(0, 1, 2), dist = d,
                   type = "functional")
  for (q in c(0, 1, 2)) {
    expect_equal(out$share[out$q == q], c(.5, .3, .2))
  }
  hd <- suppressMessages(hilldiv(counts, q = c(0, 1, 2), dist = d,
                                 type = "functional", out = "matrix"))
  expect_equal(unname(hd[1, ]), c(1, 1, 1))
})

test_that("empty functional samples have zero diversity and undefined share", {
  counts <- matrix(0, 2, 1,
                   dimnames = list(c("A", "B"), "empty"))
  d <- matrix(c(0, 1, 1, 0), 2,
              dimnames = list(c("A", "B"), c("A", "B")))
  shares <- hillshare(counts, q = c(0, 1, 2), dist = d,
                      type = "functional")
  expect_true(all(is.na(shares$share)))
  hd <- suppressMessages(hilldiv(counts, q = c(0, 1, 2), dist = d,
                                 type = "functional", out = "matrix"))
  expect_equal(unname(hd[1, ]), c(0, 0, 0))
})

test_that("hillshare validates sets and offers only qpower allocation", {
  counts <- c(A = 1, B = 2)
  expect_error(hillshare(counts, sets = list(unknown = "C")), "Unknown taxa")
  expect_error(hillshare(counts, allocation = "shapley"), "qpower")
  expect_error(hillshare(counts, q = NA_real_), "finite")
})
