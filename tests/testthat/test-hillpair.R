counts <- matrix(
  c(10, 0, 5,
    2, 8, 1,
    3, 4, 0,
    6, 2, 7),
  nrow = 3,
  dimnames = list(c("t1", "t2", "t3"), c("s1", "s2", "s3", "s4"))
)
tree <- ape::read.tree(text = "((t1:1,t2:1):1,t3:2);")
dist <- matrix(c(0, 0.4, 0.9, 0.4, 0, 0.6, 0.9, 0.6, 0), nrow = 3,
               dimnames = list(rownames(counts), rownames(counts)))
q_set <- c(0, 1, 2)

test_that("pairwise dissimilarities match hilldiss on each pair", {
  dp <- suppressMessages(hillpair(counts, q = q_set, metric = c("C", "U")))
  pairs <- utils::combn(ncol(counts), 2, simplify = FALSE)
  for (idx in pairs) {
    ref <- suppressMessages(
      hilldiss(counts[, idx], q = q_set, metric = c("C", "U")))
    s <- colnames(counts)[idx]
    for (qi in seq_along(q_set)) {
      for (m in c("C", "U")) {
        got <- as.matrix(dp[[paste0("q", q_set[qi], m)]])[s[1], s[2]]
        expect_equal(got, unname(ref[qi, m]))
      }
    }
  }
})

test_that("a single order/metric returns one dist object", {
  one <- suppressMessages(hillpair(counts, q = 1, metric = "C"))
  expect_s3_class(one, "dist")
  expect_equal(attr(one, "Size"), ncol(counts))
  expect_equal(attr(one, "Labels"), colnames(counts))
})

test_that("tibble output is long and consistent with the dist output", {
  tb <- suppressMessages(
    hillpair(counts, q = q_set, metric = c("S", "V"), out = "tibble"))
  np <- choose(ncol(counts), 2)
  expect_equal(nrow(tb), np * length(q_set) * 2)
  expect_named(tb, c("first", "second", "q", "metric", "value"))

  dp <- suppressMessages(hillpair(counts, q = q_set, metric = c("S", "V")))
  row <- tb[tb$first == "s1" & tb$second == "s3" & tb$q == 2 &
              tb$metric == "V", "value"]
  expect_equal(row, as.matrix(dp$q2V)["s1", "s3"])
})

test_that("pairwise works for phylogenetic and functional types", {
  dphylo <- suppressMessages(hillpair(counts, q = 1, metric = "C", tree = tree))
  dfunc <- suppressMessages(hillpair(counts, q = 1, metric = "C", dist = dist))
  expect_s3_class(dphylo, "dist")
  expect_s3_class(dfunc, "dist")

  ref <- suppressMessages(
    hilldiss(counts[, c(1, 2)], q = 1, metric = c("C", "U"), tree = tree))
  expect_equal(as.matrix(dphylo)["s1", "s2"], unname(ref[1, "C"]))
})

test_that("identical samples have zero pairwise dissimilarity", {
  twin <- cbind(a = counts[, 1], b = counts[, 1])
  rownames(twin) <- rownames(counts)
  dp <- suppressMessages(hillpair(twin, q = q_set))
  expect_true(all(vapply(dp, function(d) all(abs(d) < 1e-9), logical(1))))
})

test_that("pairwise needs at least two samples", {
  expect_error(suppressMessages(
    hillpair(counts[, 1, drop = FALSE], q = 0)))
})
