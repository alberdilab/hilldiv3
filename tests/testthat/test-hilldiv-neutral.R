make_counts <- function() {
  matrix(c(10, 0, 5, 2, 8, 1), nrow = 3,
         dimnames = list(c("t1", "t2", "t3"), c("s1", "s2")))
}

test_that("q = 0 returns observed richness", {
  counts <- make_counts()
  out <- suppressMessages(hilldiv(counts, q = 0, out = "matrix"))
  expect_equal(unname(out[, "q0"]), c(2, 3)) # s1 has 2 taxa, s2 has 3
})

test_that("q = 2 returns inverse Simpson", {
  counts <- make_counts()
  out <- suppressMessages(hilldiv(counts, q = 2, out = "matrix"))
  p <- tss(counts)
  expected <- 1 / colSums(p^2)
  expect_equal(unname(out[, "q2"]), unname(expected))
})

test_that("q = 1 returns exp(Shannon)", {
  counts <- make_counts()
  out <- suppressMessages(hilldiv(counts, q = 1, out = "matrix"))
  p <- tss(counts)
  shannon <- apply(p, 2, function(x) -sum(x[x > 0] * log(x[x > 0])))
  expect_equal(unname(out[, "q1"]), unname(exp(shannon)))
})

test_that("Hill numbers are non-increasing in q", {
  counts <- make_counts()
  out <- suppressMessages(hilldiv(counts, q = c(0, 1, 2), out = "matrix"))
  expect_true(all(out[, "q0"] >= out[, "q1"]))
  expect_true(all(out[, "q1"] >= out[, "q2"]))
})

test_that("a single sample (named vector) is handled", {
  out <- suppressMessages(hilldiv(c(t1 = 10, t2 = 2, t3 = 3), q = c(0, 1, 2),
                                  out = "matrix"))
  expect_equal(dim(out), c(1L, 3L))
  expect_equal(colnames(out), c("q0", "q1", "q2"))
  expect_equal(unname(out[, "q0"]), 3)
})

test_that("a single sample matches the same data as a one-column matrix", {
  vec <- c(t1 = 10, t2 = 2, t3 = 3)
  mat <- matrix(vec, ncol = 1, dimnames = list(names(vec), "s1"))
  out_vec <- suppressMessages(hilldiv(vec, out = "matrix"))
  out_mat <- suppressMessages(hilldiv(mat, out = "matrix"))
  expect_equal(unname(out_vec), unname(out_mat))
})

test_that("a leading character column is used as taxa names", {
  counts <- make_counts()
  # Column name is arbitrary; the first column is what matters.
  df <- tibble::as_tibble(counts, rownames = "anything")
  expect_message(
    out <- hilldiv(df, q = c(0, 1, 2), out = "matrix"),
    "first column"
  )
  expect_equal(out, suppressMessages(hilldiv(counts, q = c(0, 1, 2),
                                             out = "matrix")))
})

test_that("multiple non-numeric columns ask the user to fix the input", {
  df <- tibble::tibble(genome = c("a", "b"), label = c("x", "y"),
                       s1 = c(1, 2))
  expect_error(hilldiv(df), "Non-numeric columns found")
})

test_that("a non-numeric column that is not first is rejected", {
  df <- tibble::tibble(s1 = c(1, 2), label = c("x", "y"))
  expect_error(hilldiv(df), "Non-numeric column")
})
