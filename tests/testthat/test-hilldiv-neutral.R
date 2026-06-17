make_counts <- function() {
  matrix(c(10, 0, 5, 2, 8, 1), nrow = 3,
         dimnames = list(c("t1", "t2", "t3"), c("s1", "s2")))
}

test_that("q = 0 returns observed richness", {
  counts <- make_counts()
  out <- suppressMessages(hilldiv(counts, q = 0))
  expect_equal(unname(out["q0", ]), c(2, 3)) # s1 has 2 taxa, s2 has 3
})

test_that("q = 2 returns inverse Simpson", {
  counts <- make_counts()
  out <- suppressMessages(hilldiv(counts, q = 2))
  p <- tss(counts)
  expected <- 1 / colSums(p^2)
  expect_equal(unname(out["q2", ]), unname(expected))
})

test_that("q = 1 returns exp(Shannon)", {
  counts <- make_counts()
  out <- suppressMessages(hilldiv(counts, q = 1))
  p <- tss(counts)
  shannon <- apply(p, 2, function(x) -sum(x[x > 0] * log(x[x > 0])))
  expect_equal(unname(out["q1", ]), unname(exp(shannon)))
})

test_that("Hill numbers are non-increasing in q", {
  counts <- make_counts()
  out <- suppressMessages(hilldiv(counts, q = c(0, 1, 2)))
  expect_true(all(out["q0", ] >= out["q1", ]))
  expect_true(all(out["q1", ] >= out["q2", ]))
})
