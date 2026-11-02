test_that("hierarchical overlap matches pooled-unit overlap at the top", {
  set.seed(42)
  counts <- matrix(rpois(80, 4), nrow = 10,
                   dimnames = list(paste0("t", 1:10), paste0("s", 1:8)))
  md <- data.frame(region = rep(c("N", "S"), each = 4),
                   site = rep(letters[1:4], each = 2),
                   row.names = colnames(counts))
  p <- tss(counts)
  pooled <- cbind(N = rowMeans(p[, 1:4, drop = FALSE]),
                  S = rowMeans(p[, 5:8, drop = FALSE]))
  ref <- suppressMessages(hilldiss(pooled, q = c(1, 2)))
  got <- suppressMessages(hilldiss(counts, q = c(1, 2),
                                    hierarchy = ~ region / site,
                                    metadata = md))
  top <- got[got$scale == "total", ]
  for (qv in c(1, 2)) {
    for (m in c("S", "C", "U", "V")) {
      expect_equal(top$value[top$q == qv & top$metric == m],
                   ref$value[ref$q == qv & ref$metric == m],
                   tolerance = 1e-10)
    }
  }
  sim <- suppressMessages(hillsim(counts, q = c(1, 2),
                                   hierarchy = ~ region / site,
                                   metadata = md))
  expect_s3_class(got, "hill_dissimilarity")
  expect_s3_class(sim, "hill_similarity")
  expect_equal(sim$value, 1 - got$value, tolerance = 1e-12)
  expect_named(got, c("q", "scale", "metric", "value"))

  mat <- suppressMessages(hilldiss(counts, q = c(1, 2), metric = "C",
                                    hierarchy = ~ region / site,
                                    metadata = md, out = "matrix"))
  expect_equal(dim(mat), c(6L, 1L))
  expect_equal(rownames(mat), paste0("q", rep(c(1, 2), each = 3), "_",
                                       rep(c("site", "region", "total"), 2)))
})

test_that("hierarchical overlap supports tau grids and plots", {
  counts <- gut_counts
  d <- traits2dist(gut_traits)
  md <- data.frame(region = rep(c("N", "S"), each = 6),
                   row.names = colnames(counts))
  tau_grid <- c(0.25, 0.5)
  args <- list(data = counts, dist = d, hierarchy = ~ region,
               metadata = md, q = c(0, 1, 2), weights = seq_len(ncol(counts)))
  got <- suppressMessages(do.call(hilldiss, c(args, list(tau = tau_grid))))
  sim <- suppressMessages(do.call(hillsim, c(args, list(tau = tau_grid))))
  expect_named(got, c("q", "scale", "tau", "metric", "value"))
  expect_equal(sim$value, 1 - got$value, tolerance = 1e-12)
  expect_true(all(got$value >= -1e-10 & got$value <= 1 + 1e-10))

  matrices <- suppressMessages(do.call(hilldiss, c(
    args, list(tau = tau_grid, out = "matrix")
  )))
  expect_named(matrices, paste0("tau_", tau_grid))
  for (i in seq_along(tau_grid)) {
    single <- suppressMessages(do.call(hilldiss, c(
      args, list(tau = tau_grid[i])
    )))
    expect_false("tau" %in% names(single))
    expect_equal(got$value[got$tau == tau_grid[i]], single$value,
                 tolerance = 1e-12)
    one_matrix <- suppressMessages(do.call(hilldiss, c(
      args, list(tau = tau_grid[i], out = "matrix")
    )))
    expect_identical(matrices[[i]], one_matrix)
  }

  grDevices::pdf(file = "/dev/null")
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(plot(got))
  expect_silent(plot(sim))
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    expect_s3_class(ggplot2::autoplot(got), "ggplot")
  }
})

test_that("hierarchical overlap validates weights, depth and tau", {
  counts <- gut_counts
  md <- data.frame(region = rep(c("N", "S"), each = 6),
                   row.names = colnames(counts))
  expect_error(suppressMessages(hilldiss(counts, weights = "equal")),
               "hierarchy")
  expect_error(suppressMessages(hillsim(counts, depth = 0.5)), "hierarchy")
  expect_error(suppressMessages(hilldiss(counts, hierarchy = ~ region,
                                         metadata = md, depth = 0.5)),
               "phylogenetic")
  expect_error(suppressMessages(hilldiss(counts, hierarchy = ~ region,
                                         metadata = md, tau = c(0.2, 0.5))),
               "functional")

  custom <- suppressMessages(hilldiss(counts, tree = gut_tree,
                                       hierarchy = ~ region, metadata = md,
                                       depth = 0.5,
                                       weights = seq_len(ncol(counts))))
  full <- suppressMessages(hilldiss(counts, tree = gut_tree,
                                     hierarchy = ~ region, metadata = md,
                                     weights = seq_len(ncol(counts))))
  expect_false(isTRUE(all.equal(custom$value, full$value)))
  expect_true(all(custom$value >= -1e-10 & custom$value <= 1 + 1e-10))
})
