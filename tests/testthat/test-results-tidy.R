counts <- matrix(
  c(10, 0, 5, 2,
    8, 1, 3, 4,
    0, 6, 2, 7),
  nrow = 4,
  dimnames = list(c("t1", "t2", "t3", "t4"), c("s1", "s2", "s3"))
)
tree <- ape::read.tree(text = "((t1:1,t2:1):1,(t3:1,t4:2):1.5);")

test_that("hilldiv returns a tidy hill_diversity by default", {
  out <- suppressMessages(hilldiv(counts, q = c(0, 1, 2)))
  expect_s3_class(out, "hill_diversity")
  expect_s3_class(out, "hill_result")
  expect_named(out, c("q", "sample", "value"))
  expect_equal(nrow(out), 3L * ncol(counts))
  expect_equal(attr(out, "hill_type"), "neutral")

  # tidy and matrix agree cell by cell.
  mat <- suppressMessages(hilldiv(counts, q = c(0, 1, 2), out = "matrix"))
  v <- out$value[out$q == 1 & out$sample == "s2"]
  expect_equal(v, unname(mat["s2", "q1"]))
})

test_that("hillpart tidy output uses a component column", {
  out <- suppressMessages(hillpart(counts, q = c(0, 1)))
  expect_s3_class(out, "hill_partition")
  expect_named(out, c("q", "component", "value"))
  expect_setequal(unique(out$component), c("alpha", "gamma", "beta"))
})

test_that("hilldiss/hillsim tidy output uses a metric column and complements", {
  d <- suppressMessages(hilldiss(counts, q = c(0, 1, 2)))
  s <- suppressMessages(hillsim(counts, q = c(0, 1, 2)))
  expect_s3_class(d, "hill_dissimilarity")
  expect_s3_class(s, "hill_similarity")
  expect_named(d, c("q", "metric", "value"))
  # similarity = 1 - dissimilarity row for row (same ordering).
  expect_equal(d$value + s$value, rep(1, nrow(d)))
})

test_that("hilleven tidy output is bounded in (0, 1]", {
  out <- suppressMessages(hilleven(counts, q = c(1, 2)))
  expect_s3_class(out, "hill_evenness")
  expect_named(out, c("q", "sample", "value"))
  expect_true(all(out$value > 0 & out$value <= 1 + 1e-9))
})

test_that("hillred tidy output is one row per q", {
  out <- suppressWarnings(
    suppressMessages(hillred(counts, q = c(0, 1, 2), tree = tree)))
  expect_s3_class(out, "hill_redundancy")
  expect_named(out, c("q", "redundancy", "a", "b", "c"))
  expect_equal(nrow(out), 3L)
})

test_that("print method returns its input invisibly", {
  out <- suppressMessages(hilldiv(counts, q = 1))
  expect_output(print(out), "hilldiv3 result")
  expect_invisible(print(out))
})

test_that("plot methods draw without error for each result type", {
  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp)
  on.exit({ grDevices::dev.off(); unlink(tmp) })
  expect_invisible(plot(suppressMessages(hilldiv(counts, q = c(0, 1, 2)))))
  expect_invisible(plot(suppressMessages(hillpart(counts, q = c(0, 1, 2)))))
  expect_invisible(plot(suppressMessages(hilldiss(counts, q = c(0, 1, 2)))))
  expect_invisible(plot(suppressMessages(hilleven(counts, q = c(1, 2)))))
})

test_that("autoplot returns a ggplot when ggplot2 is available", {
  skip_if_not_installed("ggplot2")
  out <- suppressMessages(hilldiv(counts, q = c(0, 1, 2)))
  p <- ggplot2::autoplot(out)
  expect_s3_class(p, "ggplot")
})
