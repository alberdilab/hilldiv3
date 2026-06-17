test_that("tss normalises a vector to sum 1", {
  v <- c(a = 1, b = 3)
  expect_equal(tss(v), c(a = 0.25, b = 0.75))
})

test_that("tss handles all-zero input without NaN", {
  expect_equal(tss(c(0, 0, 0)), c(0, 0, 0))
})

test_that("tss normalises matrix columns to sum 1", {
  m <- matrix(c(1, 1, 0, 2, 0, 2), nrow = 3)
  out <- tss(m)
  expect_equal(unname(colSums(out)), c(1, 1))
})

test_that("tss maps all-zero columns to 0 (not NaN)", {
  m <- matrix(c(0, 0, 1, 1), nrow = 2)
  out <- tss(m)
  expect_false(any(is.nan(out)))
  expect_equal(unname(out[, 1]), c(0, 0))
})
