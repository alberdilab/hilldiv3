q_set <- c(0, 0.5, 1, 2)

# 8 samples nested as: region {N, S}, each region two sites, each site two
# samples (a balanced 2-2-2 design), plus a deterministic count table.
md <- data.frame(
  region = rep(c("N", "S"), each = 4),
  site   = rep(c("a", "b", "c", "d"), each = 2),
  row.names = paste0("s", 1:8)
)
set.seed(42)
tab <- matrix(rpois(10 * 8, lambda = 4), nrow = 10,
              dimnames = list(paste0("t", 1:10), paste0("s", 1:8)))

hier <- function(...) suppressMessages(hillpart(tab, metadata = md, ...))

test_that("the nested chain telescopes: gamma == alpha * prod(beta)", {
  df <- hier(q = q_set, hierarchy = ~ region / site)
  for (qv in q_set) {
    rows <- df[df$q == qv, ]
    alpha <- rows$diversity[rows$scale == "sample"]
    gamma <- rows$diversity[rows$scale == "total"]
    betas <- rows$beta[!is.na(rows$beta)]
    expect_equal(prod(betas), gamma / alpha, tolerance = 1e-10)
    expect_equal(alpha * prod(betas), gamma, tolerance = 1e-10)
  }
})

test_that("finest alpha and pooled gamma match the single-level engine", {
  std <- suppressMessages(hillpart(tab, q = q_set, out = "matrix"))
  df <- hier(q = q_set, hierarchy = ~ region / site)
  for (i in seq_along(q_set)) {
    rows <- df[df$q == q_set[i], ]
    expect_equal(rows$diversity[rows$scale == "sample"],
                 unname(std[i, "alpha"]), tolerance = 1e-10)
    expect_equal(rows$diversity[rows$scale == "total"],
                 unname(std[i, "gamma"]), tolerance = 1e-10)
    # Overall beta (product) equals the single-level Whittaker beta.
    betas <- rows$beta[!is.na(rows$beta)]
    expect_equal(prod(betas), unname(std[i, "beta"]), tolerance = 1e-10)
  }
})

test_that("every beta is >= 1 and diversity is non-decreasing across scales", {
  df <- hier(q = q_set, hierarchy = ~ region / site)
  expect_true(all(df$beta >= 1 - 1e-9, na.rm = TRUE))
  for (qv in q_set) {
    rows <- df[df$q == qv, ]
    rows <- rows[order(rows$scale), ]
    expect_true(all(diff(rows$diversity) >= -1e-9))
  }
})

test_that("each beta is bounded by the number of child units it pools", {
  # sample(8) -> site(4) -> region(2) -> total(1): the per-step turnover cannot
  # exceed the number of finer units merged into each coarser unit.
  df <- hier(q = 0, hierarchy = ~ region / site)
  rows <- df[order(df$scale), ]
  # n_units finest->coarsest: 8, 4, 2, 1.
  expect_equal(rows$n_units, c(8, 4, 2, 1))
  # beta_site <= 8/4=2 per site, beta_region <= 4/2, beta_total <= 2/1.
  expect_true(rows$beta[rows$scale == "site"]   <= 2 + 1e-9)
  expect_true(rows$beta[rows$scale == "region"] <= 2 + 1e-9)
  expect_true(rows$beta[rows$scale == "total"]  <= 2 + 1e-9)
})

test_that("identical samples give no turnover (all betas = 1)", {
  one <- tab[, 1]
  twin <- matrix(rep(one, 8), nrow = nrow(tab),
                 dimnames = list(rownames(tab), paste0("s", 1:8)))
  df <- suppressMessages(
    hillpart(twin, metadata = md, q = q_set, hierarchy = ~ region / site)
  )
  expect_equal(df$beta[!is.na(df$beta)],
               rep(1, sum(!is.na(df$beta))), tolerance = 1e-9)
})

test_that("a single grouping level works and still telescopes", {
  df <- hier(q = q_set, hierarchy = ~ region)
  expect_setequal(levels(df$scale), c("sample", "region", "total"))
  for (qv in q_set) {
    rows <- df[df$q == qv, ]
    alpha <- rows$diversity[rows$scale == "sample"]
    gamma <- rows$diversity[rows$scale == "total"]
    expect_equal(alpha * prod(rows$beta, na.rm = TRUE), gamma, tolerance = 1e-10)
  }
})

test_that("deep hierarchies with arbitrarily many levels telescope", {
  md4 <- data.frame(
    a = rep(c("A", "B"), each = 8),
    b = rep(c("p", "q", "r", "s"), each = 4),
    c = rep(letters[1:8], each = 2),
    row.names = paste0("x", 1:16)
  )
  set.seed(7)
  tab4 <- matrix(rpois(15 * 16, 5), nrow = 15,
                 dimnames = list(paste0("t", 1:15), paste0("x", 1:16)))
  df <- suppressMessages(
    hillpart(tab4, metadata = md4, q = q_set, hierarchy = ~ a / b / c)
  )
  expect_equal(nlevels(df$scale), 5L)  # sample < c < b < a < total
  for (qv in q_set) {
    rows <- df[df$q == qv, ]
    alpha <- rows$diversity[rows$scale == "sample"]
    gamma <- rows$diversity[rows$scale == "total"]
    expect_equal(alpha * prod(rows$beta, na.rm = TRUE), gamma, tolerance = 1e-9)
  }
})

test_that("nesting keeps reused child labels distinct across parents", {
  # Both regions reuse site labels "a"/"b"; they must not be pooled together.
  md2 <- data.frame(
    region = rep(c("N", "S"), each = 4),
    site   = rep(c("a", "b"), times = 4),
    row.names = paste0("s", 1:8)
  )
  df <- suppressMessages(
    hillpart(tab, metadata = md2, q = 0, hierarchy = ~ region / site)
  )
  expect_equal(df$n_units[df$scale == "site"], 4)  # 2 regions x 2 sites
})

test_that("hierarchy resolves variables from the calling environment", {
  region <- md$region
  site <- md$site
  df_env <- suppressMessages(
    hillpart(tab, q = 0, hierarchy = ~ region / site)
  )
  df_md <- hier(q = 0, hierarchy = ~ region / site)
  expect_equal(df_env$diversity, df_md$diversity, tolerance = 1e-12)
})

test_that("matrix output has alpha, per-level betas and gamma", {
  m <- suppressMessages(
    hillpart(tab, metadata = md, q = q_set, hierarchy = ~ region / site,
             out = "matrix")
  )
  expect_equal(colnames(m),
               c("alpha", "beta_site", "beta_region", "beta_total", "gamma"))
  expect_equal(rownames(m), paste0("q", q_set))
  # Row-wise telescoping in the wide form: alpha * prod(betas) == gamma.
  beta_cols <- c("beta_site", "beta_region", "beta_total")
  for (i in seq_along(q_set)) {
    expect_equal(m[i, "alpha"] * prod(m[i, beta_cols]),
                 m[i, "gamma"], tolerance = 1e-10)
  }
})

test_that("type is auto-detected from inputs for hierarchical partitioning", {
  tree <- ape::rtree(nrow(tab), tip.label = rownames(tab))
  df <- suppressMessages(
    hillpart(tab, metadata = md, tree = tree, hierarchy = ~ region / site)
  )
  expect_identical(attr(df, "hill_type"), "phylogenetic")
  expect_s3_class(df, "hill_hierarchy")
})

test_that("a length-mismatched hierarchy variable is rejected", {
  bad <- data.frame(region = rep("N", 3), row.names = paste0("z", 1:3))
  expect_error(
    suppressMessages(hillpart(tab, metadata = bad, hierarchy = ~ region)),
    "length"
  )
})

# ---- Phylogenetic ----------------------------------------------------------

test_that("phylogenetic chain telescopes and matches the single-level engine", {
  set.seed(11)
  tree <- ape::compute.brlen(ape::rtree(nrow(tab), tip.label = rownames(tab)),
                             method = "Grafen")  # ultrametric
  std <- suppressMessages(hillpart(tab, q = q_set, tree = tree, out = "matrix"))
  df <- suppressMessages(
    hillpart(tab, metadata = md, q = q_set, tree = tree,
             hierarchy = ~ region / site)
  )
  for (i in seq_along(q_set)) {
    rows <- df[df$q == q_set[i], ]
    alpha <- rows$diversity[rows$scale == "sample"]
    gamma <- rows$diversity[rows$scale == "total"]
    betas <- rows$beta[!is.na(rows$beta)]
    expect_equal(alpha * prod(betas), gamma, tolerance = 1e-9)
    # Endpoints reproduce the single-level phylogenetic engine (ultrametric).
    expect_equal(alpha, unname(std[i, "alpha"]), tolerance = 1e-9)
    expect_equal(gamma, unname(std[i, "gamma"]), tolerance = 1e-9)
    expect_equal(prod(betas), unname(std[i, "beta"]), tolerance = 1e-9)
    expect_true(all(betas >= 1 - 1e-9))
  }
})

test_that("phylogenetic hierarchy telescopes on a non-ultrametric tree", {
  set.seed(12)
  tree <- ape::rtree(nrow(tab), tip.label = rownames(tab))  # not ultrametric
  df <- suppressMessages(
    hillpart(tab, metadata = md, q = q_set, tree = tree,
             hierarchy = ~ region / site)
  )
  for (qv in q_set) {
    rows <- df[df$q == qv, ]
    alpha <- rows$diversity[rows$scale == "sample"]
    gamma <- rows$diversity[rows$scale == "total"]
    expect_equal(alpha * prod(rows$beta, na.rm = TRUE), gamma, tolerance = 1e-9)
    expect_true(all(rows$beta >= 1 - 1e-9, na.rm = TRUE))
  }
})

# ---- Functional ------------------------------------------------------------

test_that("functional chain telescopes and matches the single-level engine", {
  set.seed(13)
  d <- as.matrix(dist(matrix(runif(nrow(tab) * 3), nrow = nrow(tab))))
  dimnames(d) <- list(rownames(tab), rownames(tab))
  std <- suppressMessages(hillpart(tab, q = q_set, dist = d, out = "matrix"))
  df <- suppressMessages(
    hillpart(tab, metadata = md, q = q_set, dist = d,
             hierarchy = ~ region / site)
  )
  for (i in seq_along(q_set)) {
    rows <- df[df$q == q_set[i], ]
    alpha <- rows$diversity[rows$scale == "sample"]
    gamma <- rows$diversity[rows$scale == "total"]
    betas <- rows$beta[!is.na(rows$beta)]
    expect_equal(alpha * prod(betas), gamma, tolerance = 1e-9)
    expect_equal(alpha, unname(std[i, "alpha"]), tolerance = 1e-9)
    expect_equal(gamma, unname(std[i, "gamma"]), tolerance = 1e-9)
    expect_equal(prod(betas), unname(std[i, "beta"]), tolerance = 1e-9)
    expect_true(all(betas >= 1 - 1e-9))
  }
})

test_that("functional respects a custom tau and stays consistent across scales", {
  set.seed(14)
  d <- as.matrix(dist(matrix(runif(nrow(tab) * 3), nrow = nrow(tab))))
  dimnames(d) <- list(rownames(tab), rownames(tab))
  df <- suppressMessages(
    hillpart(tab, metadata = md, q = q_set, dist = d, tau = 0.5,
             hierarchy = ~ region / site)
  )
  for (qv in q_set) {
    rows <- df[df$q == qv, ]
    alpha <- rows$diversity[rows$scale == "sample"]
    gamma <- rows$diversity[rows$scale == "total"]
    expect_equal(alpha * prod(rows$beta, na.rm = TRUE), gamma, tolerance = 1e-9)
  }
})
