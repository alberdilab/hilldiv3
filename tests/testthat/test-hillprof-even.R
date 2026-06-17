counts <- matrix(
  c(10, 0, 5, 2,
    8, 1, 3, 4,
    0, 6, 2, 7),
  nrow = 4,
  dimnames = list(c("t1", "t2", "t3", "t4"), c("s1", "s2", "s3"))
)
tree <- ape::read.tree(text = "((t1:1,t2:1):1,(t3:1,t4:2):1.5);")

test_that("hillprof matrix equals hilldiv across orders", {
  q <- c(0, 0.5, 1, 2, 3)
  prof <- suppressMessages(hillprof(counts, q = q, out = "matrix"))
  ref <- suppressMessages(hilldiv(counts, q = q, out = "matrix"))
  expect_equal(prof, ref)
})

test_that("hillprof tibble is long-format and consistent with the matrix", {
  q <- c(0, 1, 2)
  prof <- suppressMessages(hillprof(counts, q = q))
  expect_s3_class(prof, "hill_profile")
  expect_named(prof, c("q", "sample", "value"))
  expect_equal(nrow(prof), length(q) * ncol(counts))

  mat <- suppressMessages(hillprof(counts, q = q, out = "matrix"))
  val <- prof[prof$sample == "s2" & prof$q == 1, "value"]
  expect_equal(val, mat["q1", "s2"])
})

test_that("hill_profile has a working plot method", {
  prof <- suppressMessages(hillprof(counts, q = seq(0, 3, 0.5)))
  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp)
  on.exit({ grDevices::dev.off(); unlink(tmp) })
  expect_invisible(plot(prof))
})

test_that("hilleven returns qD / 0D in (0, 1]", {
  q <- c(1, 2)
  ev <- suppressMessages(hilleven(counts, q = q, out = "matrix"))
  qd <- suppressMessages(hilldiv(counts, q = q, out = "matrix"))
  q0 <- suppressMessages(hilldiv(counts, q = 0, out = "matrix"))
  expect_equal(ev, sweep(qd, 2, q0[1, ], "/"))
  expect_true(all(ev > 0 & ev <= 1 + 1e-9))
})

test_that("a perfectly even community has evenness 1", {
  even <- matrix(rep(5, 9), nrow = 3,
                 dimnames = list(c("t1", "t2", "t3"), c("s1", "s2", "s3")))
  ev <- suppressMessages(hilleven(even, q = c(1, 2), out = "matrix"))
  expect_equal(unname(ev), matrix(1, 2, 3))
})

test_that("profile and evenness route through phylogenetic type", {
  prof <- suppressMessages(hillprof(counts, q = c(0, 1), tree = tree,
                                    out = "matrix"))
  ref <- suppressMessages(hilldiv(counts, q = c(0, 1), tree = tree,
                                  out = "matrix"))
  expect_equal(prof, ref)

  ev <- suppressMessages(hilleven(counts, q = 2, tree = tree, out = "matrix"))
  expect_equal(dim(ev), c(1, ncol(counts)))
})
