test_that("neutral contributions add to pairwise and collective dissimilarity", {
  counts <- matrix(c(9, 1, 0, 0, 1, 9, 4, 1, 5), nrow = 3,
                   dimnames = list(c("A", "B", "C"), c("s1", "s2", "s3")))
  for (q in c(0, 1, 2)) {
    got <- hillcontrib(counts, q = q)
    pair <- suppressMessages(hillpair(counts, q = q, out = "tibble"))
    for (i in seq_len(nrow(pair))) {
      rows <- got[got$first == pair$first[i] &
                    got$second == pair$second[i] &
                    got$metric == pair$metric[i], ]
      expect_equal(sum(rows$contribution), pair$value[i], tolerance = 1e-12)
      if (pair$value[i] > 0) expect_equal(sum(rows$share), 1)
      expect_equal(rows$rank, seq_len(nrow(rows)))
      expect_true(all(diff(rows$contribution) <= 0))
    }

    collective <- hillcontrib(counts, q = q, by = "collective")
    diss <- suppressMessages(hilldiss(counts, q = q))
    for (i in seq_len(nrow(diss))) {
      rows <- collective[collective$metric == diss$metric[i], ]
      expect_equal(sum(rows$contribution), diss$value[i], tolerance = 1e-12)
    }
  }
})

test_that("unique taxa lead pairwise turnover, and zero turnover has undefined shares", {
  counts <- matrix(c(9, 1, 0, 0, 1, 9), nrow = 3,
                   dimnames = list(c("A", "B", "C"), c("s1", "s2")))
  q0 <- hillcontrib(counts, q = 0, metric = "C")
  expect_equal(q0$set, c("A", "C", "B"))
  expect_equal(q0$contribution[1:2], rep(q0$contribution[1], 2))
  expect_equal(q0$contribution[3], 0)
  for (q in c(1, 2)) {
    got <- hillcontrib(counts, q = q, metric = "C")
    expect_equal(got$set, c("A", "C", "B"))
    expect_equal(got$share, c(.5, .5, 0), tolerance = 1e-12)
  }

  same <- matrix(c(9, 1, 0, 9, 1, 0), nrow = 3,
                 dimnames = list(c("A", "B", "C"), c("s1", "s2")))
  zero <- hillcontrib(same, q = c(0, 1, 2), metric = "C")
  expect_equal(zero$contribution, rep(0, nrow(zero)), tolerance = 1e-12)
  expect_true(all(is.na(zero$share)))
  expect_true(all(is.na(zero$rank)))
})

test_that("named groups and sets aggregate taxon contributions", {
  counts <- matrix(c(9, 1, 0, 0, 1, 9, 4, 1, 5), nrow = 3,
                   dimnames = list(c("A", "B", "C"), c("s1", "s2", "s3")))
  groups <- list(early = c("s1", "s2"), late = c("s2", "s3"))
  full <- hillcontrib(counts, q = 1, metric = "C", by = "collective",
                      groups = groups)
  sets <- hillcontrib(counts, q = 1, metric = "C", by = "collective",
                      groups = groups, sets = list(AC = c("A", "C"), B = "B"))
  for (g in names(groups)) {
    expect_equal(sets$contribution[sets$group == g & sets$set == "AC"],
                 sum(full$contribution[full$group == g & full$set %in% c("A", "C")]))
  }
  expect_true(all(is.na(sets$first)))
  expect_true(all(is.na(sets$second)))
})

test_that("phylogenetic and functional contributions match their beta engines", {
  counts <- matrix(c(9, 1, 0, 0, 1, 9, 4, 1, 5), nrow = 3,
                   dimnames = list(c("A", "B", "C"), c("s1", "s2", "s3")))
  tree <- ape::read.tree(text = "((A:1,B:1):1,C:2);")
  d <- matrix(c(0, .5, 1, .5, 0, .5, 1, .5, 0), 3,
              dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  for (kind in c("phylogenetic", "functional")) {
    extra <- if (kind == "phylogenetic") list(tree = tree) else list(dist = d)
    got <- do.call(hillcontrib, c(list(data = counts, q = c(0, 1, 2),
                                      metric = "C", by = "collective"), extra))
    ref <- suppressMessages(do.call(hilldiss,
                                    c(list(data = counts, q = c(0, 1, 2),
                                           metric = "C"), extra)))
    for (q in c(0, 1, 2)) {
      expect_equal(sum(got$contribution[got$q == q]),
                   ref$value[ref$q == q], tolerance = 1e-12)
    }
    expect_true(all(got$contribution >= 0))
  }
})

test_that("phylogenetic star and functional identity agree with neutral ranking", {
  counts <- matrix(c(9, 1, 0, 0, 1, 9), nrow = 3,
                   dimnames = list(c("A", "B", "C"), c("s1", "s2")))
  tree <- ape::read.tree(text = "(A:1,B:1,C:1);")
  d <- matrix(1, 3, 3, dimnames = list(rownames(counts), rownames(counts)))
  diag(d) <- 0
  base <- hillcontrib(counts, q = c(0, 1, 2), metric = "C")
  phy <- hillcontrib(counts, q = c(0, 1, 2), metric = "C", tree = tree)
  fun <- hillcontrib(counts, q = c(0, 1, 2), metric = "C", dist = d)
  expect_equal(phy$set, base$set)
  expect_equal(phy$contribution, base$contribution, tolerance = 1e-12)
  expect_equal(fun$set, base$set)
  expect_equal(fun$contribution, base$contribution, tolerance = 1e-12)
})

test_that("pairwise attribution follows phylogenetic and functional engines", {
  counts <- matrix(c(8, 2, 0, 1, 4, 5, 3, 1, 2), nrow = 3,
                   dimnames = list(c("A", "B", "C"), c("s1", "s2", "s3")))
  tree <- ape::read.tree(text = "((A:1,B:2):.5,C:1);")
  d <- matrix(c(0, .5, 1, .5, 0, .5, 1, .5, 0), 3,
              dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  for (extra in list(list(tree = tree), list(dist = d, tau = .75))) {
    got <- do.call(hillcontrib,
                   c(list(data = counts, q = c(0, .5, 1, 3), metric = "U"), extra))
    ref <- suppressMessages(do.call(hillpair,
                                    c(list(data = counts, q = c(0, .5, 1, 3),
                                           metric = "U", out = "tibble"), extra)))
    for (i in seq_len(nrow(ref))) {
      rows <- got[got$first == ref$first[i] &
                    got$second == ref$second[i] & got$q == ref$q[i], ]
      expect_equal(sum(rows$contribution), ref$value[i], tolerance = 1e-12)
    }
  }
})

test_that("functional all-zero distances give zero turnover for equal totals", {
  counts <- matrix(c(5, 0, 0, 5), 2,
                   dimnames = list(c("A", "B"), c("s1", "s2")))
  d <- matrix(0, 2, 2, dimnames = list(rownames(counts), rownames(counts)))
  got <- hillcontrib(counts, q = c(0, 1, 2), metric = "C", dist = d)
  expect_equal(got$contribution, rep(0, nrow(got)))
  expect_true(all(is.na(got$share)))
})

test_that("input and group validation is explicit", {
  counts <- matrix(c(1, 0, 0, 1), 2,
                   dimnames = list(c("A", "B"), c("s1", "s2")))
  expect_error(hillcontrib(counts, q = NA_real_), "finite")
  expect_error(hillcontrib(counts, groups = list(g = c("s1", "s2"))), "collective")
  expect_error(hillcontrib(counts, by = "collective",
                           groups = list(g = c("s1", "bad"))), "known sample")
  expect_error(hillcontrib(counts, sets = list(bad = "C")), "Unknown taxa")
  counts[, 2] <- 0
  expect_error(hillcontrib(counts), "positive total")
})
