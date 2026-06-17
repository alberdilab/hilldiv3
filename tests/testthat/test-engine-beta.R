test_that("similarity is the complement of dissimilarity", {
  for (q in c(0, 1, 2)) {
    d <- hilldiv3:::beta_to_dissim(beta = 1.5, N = 4, q = q)
    s <- hilldiv3:::beta_to_sim(beta = 1.5, N = 4, q = q)
    expect_equal(unname(d + s), rep(1, 4))
  }
})

test_that("identical samples (beta = 1) give zero dissimilarity", {
  for (q in c(0, 1, 2)) {
    d <- hilldiv3:::beta_to_dissim(beta = 1, N = 5, q = q)
    expect_equal(unname(d), rep(0, 4))
  }
})

test_that("maximally distinct samples (beta = N) give dissimilarity 1", {
  N <- 5
  for (q in c(0, 1, 2)) {
    d <- hilldiv3:::beta_to_dissim(beta = N, N = N, q = q)
    expect_equal(unname(d), rep(1, 4))
  }
})

test_that("q = 1 overlap limit matches the limit of nearby q", {
  near <- hilldiv3:::beta_to_dissim(beta = 1.5, N = 4, q = 1 + 1e-6)
  lim <- hilldiv3:::beta_to_dissim(beta = 1.5, N = 4, q = 1)
  expect_equal(lim[["C"]], near[["C"]], tolerance = 1e-4)
  expect_equal(lim[["U"]], near[["U"]], tolerance = 1e-4)
})
