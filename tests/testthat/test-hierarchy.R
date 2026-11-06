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

test_that("equal weights preserve the default hierarchy output", {
  default <- hier(q = q_set, hierarchy = ~ region / site)
  expect_identical(hier(q = q_set, hierarchy = ~ region / site,
                        weights = "equal"), default)
  expect_identical(hier(q = q_set, hierarchy = ~ region / site,
                        weights = rep(3, ncol(tab))), default)
  expect_identical(hier(q = q_set, hierarchy = ~ region / site,
                        weights = setNames(rep(3, ncol(tab)), colnames(tab))),
                   default)
})

test_that("unequal weights telescope and give valid betas for all types", {
  example_counts <- gut_counts
  example_md <- data.frame(
    region = rep(c("control", "treatment"), each = 6),
    site = rep(letters[1:6], each = 2),
    row.names = colnames(example_counts)
  )
  sample_weights <- seq_len(ncol(example_counts))
  args <- list(data = example_counts, hierarchy = ~ region / site,
               metadata = example_md, q = c(0, 1, 2),
               weights = sample_weights)
  inputs <- list(
    neutral = list(),
    phylogenetic = list(tree = gut_tree),
    functional = list(dist = traits2dist(gut_traits))
  )
  for (diversity_type in names(inputs)) {
    input <- inputs[[diversity_type]]
    prep <- do.call(hier_part_prep, c(
      list(counts = example_counts, type = diversity_type,
           sample_weights = sample_weights / sum(sample_weights)), input
    ))
    expect_equal(sum(prep$measure * rowSums(prep$contrib)) / prep$C,
                 1, tolerance = 1e-12)
    df <- suppressMessages(do.call(hillpart, c(args, input)))
    for (qv in c(0, 1, 2)) {
      rows <- df[df$q == qv, ]
      alpha <- rows$diversity[rows$scale == "sample"]
      gamma <- rows$diversity[rows$scale == "total"]
      expect_equal(alpha * prod(rows$beta, na.rm = TRUE), gamma,
                   tolerance = 1e-10)
      expect_true(all(rows$beta >= 1 - 1e-10, na.rm = TRUE))
    }
  }
})

test_that("balanced equal-weight beta values match the previous engine", {
  set.seed(11)
  tree <- ape::compute.brlen(ape::rtree(nrow(tab), tip.label = rownames(tab)),
                             method = "Grafen")
  set.seed(13)
  d <- as.matrix(dist(matrix(runif(nrow(tab) * 3), nrow = nrow(tab))))
  dimnames(d) <- list(rownames(tab), rownames(tab))
  inputs <- list(
    neutral = list(),
    phylogenetic = list(tree = tree),
    functional = list(dist = d)
  )
  # q = 0, 1, 2 beta_site from the prior joint-distribution engine.
  reference <- list(
    neutral = c(1.05263157894737, 1.06358818218746, 1.08895431881670),
    phylogenetic = c(1.01657458563536, 1.02731098237962,
                     1.02226525821352),
    functional = c(1, 1.02182951696134, 1.03726656329460)
  )
  for (diversity_type in names(inputs)) {
    result <- suppressMessages(do.call(hillpart, c(
      list(data = tab, metadata = md, q = c(0, 1, 2),
           hierarchy = ~ region / site, out = "matrix"),
      inputs[[diversity_type]]
    )))
    expect_equal(unname(result[, "beta_site"]), reference[[diversity_type]],
                 tolerance = 1e-12)
  }
})

test_that("linear unit weighting keeps beta >= 1 in unbalanced designs", {
  set.seed(2026)
  counts <- gut_counts
  tree <- gut_tree
  d <- traits2dist(gut_traits)
  orders <- c(0, 0.5, 0.999, 1, 2, 3)
  types <- c("neutral", "phylogenetic", "functional")
  min_beta <- setNames(rep(Inf, length(types)), types)
  min_ceiling_gap <- setNames(rep(Inf, length(types)), types)
  min_overlap <- setNames(rep(Inf, length(types)), types)
  max_overlap <- setNames(rep(-Inf, length(types)), types)
  max_gap <- setNames(rep(0, length(types)), types)

  for (design in seq_len(200)) {
    sample_weights <- rexp(ncol(counts))
    sample_weights <- sample_weights / sum(sample_weights)
    second_size <- sample(2:9, 1)
    site <- sample(rep(1:3, c(1, second_size, 11 - second_size)))
    region <- ifelse(site == 1, 1, 2)
    for (diversity_type in types) {
      prep <- hier_part_prep(counts, type = diversity_type, tree = tree,
                             dist = d, sample_weights = sample_weights)
      df <- hier_partition(prep, list(site, region), orders,
                           c("site", "region"), with_ceiling = TRUE)
      min_beta[diversity_type] <- min(min_beta[diversity_type], df$beta,
                                      na.rm = TRUE)
      transitions <- df[!is.na(df$beta), ]
      min_ceiling_gap[diversity_type] <- min(
        min_ceiling_gap[diversity_type],
        transitions$n_eff - transitions$beta
      )
      overlap <- unlist(lapply(seq_len(nrow(transitions)), function(i) {
        if (transitions$n_eff[i] <= 1 + 1e-12) return(rep(0, 4))
        beta_to_dissim(transitions$beta[i], transitions$n_eff[i],
                       transitions$q[i])
      }))
      min_overlap[diversity_type] <- min(min_overlap[diversity_type], overlap)
      max_overlap[diversity_type] <- max(max_overlap[diversity_type], overlap)
      for (qv in orders) {
        rows <- df[df$q == qv, ]
        alpha <- rows$diversity[rows$scale == "sample"]
        gamma <- rows$diversity[rows$scale == "total"]
        max_gap[diversity_type] <- max(max_gap[diversity_type],
                                       abs(alpha * prod(rows$beta, na.rm = TRUE) -
                                             gamma))
      }
    }
  }
  for (diversity_type in types) {
    expect_gte(min_beta[diversity_type], 1 - 1e-10)
    expect_gte(min_ceiling_gap[diversity_type], -1e-10)
    expect_gte(min_overlap[diversity_type], -1e-10)
    expect_lte(max_overlap[diversity_type], 1 + 1e-10)
    expect_lt(max_gap[diversity_type], 1e-10)
  }
})

test_that("balanced hierarchy ceilings equal child counts", {
  prep <- hier_part_prep(tab)
  spec <- .parse_hierarchy(~ region / site, md, colnames(tab))
  df <- hier_partition(prep, spec$groupings, c(0, 0.5, 1, 2, 3),
                       spec$level_names, with_ceiling = TRUE)
  expect_equal(df$n_eff[!is.na(df$n_eff)], rep(2, 15),
               tolerance = 1e-10)
})

test_that("the q = 1 ceiling includes child measure mass", {
  tree <- ape::read.tree(text = "(t1:3,t2:1);")
  counts <- matrix(c(1, 0, 0, 1), nrow = 2,
                   dimnames = list(c("t1", "t2"), c("s1", "s2")))
  sample_weights <- c(0.3, 0.7)
  prep <- hier_part_prep(counts, type = "phylogenetic", tree = tree,
                         sample_weights = sample_weights)
  df <- hier_partition(prep, list(c(1, 2)), c(0.999, 1), "unit",
                       with_ceiling = TRUE)
  top <- df[df$scale == "total", ]
  expect_equal(top$beta, top$n_eff, tolerance = 1e-10)
  expect_equal(top$n_eff[top$q == 1], 2.30085315338, tolerance = 1e-10)
  expect_lt(abs(diff(top$n_eff)), 1e-3)

  md <- data.frame(unit = c("a", "b"), row.names = colnames(counts))
  overlap <- suppressMessages(hilldiss(counts, tree = tree,
                                        hierarchy = ~ unit, metadata = md,
                                        weights = sample_weights,
                                        q = c(0.999, 1)))
  expect_true(all(overlap$value >= -1e-10 & overlap$value <= 1 + 1e-10))
})

test_that("named and abundance sample weights are aligned and normalised", {
  w <- seq_len(ncol(tab))
  positional <- hier(q = q_set, hierarchy = ~ region / site, weights = w)
  named <- setNames(w, colnames(tab))
  expect_identical(hier(q = q_set, hierarchy = ~ region / site,
                        weights = named[rev(seq_along(named))]), positional)
  expect_identical(hier(q = q_set, hierarchy = ~ region / site,
                        weights = 10 * w), positional)
  expect_equal(hier(q = q_set, hierarchy = ~ region / site,
                    weights = "abundance")$diversity,
               hier(q = q_set, hierarchy = ~ region / site,
                    weights = colSums(tab))$diversity)
})

test_that("sample weights are validated for hierarchical partitioning", {
  expect_error(suppressMessages(hillpart(tab, weights = "equal")), "hierarchy")
  expect_error(suppressMessages(hillpart(tab, weights = 1:8)), "hierarchy")
  expect_error(hier(hierarchy = ~ region, weights = c(1, 2)), "one finite")
  expect_error(hier(hierarchy = ~ region, weights = rep(0, 8)), "positive")
  expect_error(hier(hierarchy = ~ region,
                    weights = setNames(1:8, paste0("x", 1:8))), "match")
})

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

test_that("phylogenetic reference depth truncates one shared tree", {
  counts <- gut_counts
  tree <- gut_tree
  md <- data.frame(region = rep(c("control", "treatment"), each = 6),
                   site = rep(letters[1:6], each = 2),
                   row.names = colnames(counts))
  height <- max(ape::node.depth.edgelength(tree)[seq_along(tree$tip.label)])
  args <- list(data = counts, tree = tree, metadata = md,
               hierarchy = ~ region / site, q = c(0, 1, 2))
  default <- suppressMessages(do.call(hillpart, args))
  full_depth <- suppressMessages(do.call(hillpart,
                                         c(args, list(depth = height))))
  expect_identical(full_depth, default)

  depth <- height / 2
  prep <- hier_part_prep(counts, type = "phylogenetic", tree = tree,
                         depth = depth)
  expect_equal(prep$C, ncol(counts) * depth, tolerance = 1e-10)
  expect_true(any(prep$measure == 0))
  custom <- suppressMessages(do.call(hillpart, c(
    args, list(depth = depth, weights = seq_len(ncol(counts)))
  )))
  expect_false(isTRUE(all.equal(custom$diversity, default$diversity)))
  for (qv in c(0, 1, 2)) {
    rows <- custom[custom$q == qv, ]
    alpha <- rows$diversity[rows$scale == "sample"]
    gamma <- rows$diversity[rows$scale == "total"]
    expect_equal(alpha * prod(rows$beta, na.rm = TRUE), gamma,
                 tolerance = 1e-10)
    expect_true(all(rows$beta >= 1 - 1e-10, na.rm = TRUE))
  }
})

test_that("phylogenetic reference depth clips partial branches", {
  tree <- ape::read.tree(text = "((t1:1,t2:1):1,t3:2);")
  counts <- matrix(1, nrow = 3, ncol = 2,
                   dimnames = list(tree$tip.label, c("s1", "s2")))
  prep <- hier_part_prep(counts, type = "phylogenetic", tree = tree,
                         depth = 1.5)
  child <- tree$edge[, 2]
  expected <- ifelse(child > length(tree$tip.label), 0.5,
                     ifelse(child == match("t3", tree$tip.label), 1.5, 1))
  expect_equal(prep$measure, expected)
  expect_equal(prep$C, 3)
})

test_that("phylogenetic reference depth is validated", {
  counts <- gut_counts
  tree <- gut_tree
  md <- data.frame(region = rep(c("control", "treatment"), each = 6),
                   row.names = colnames(counts))
  height <- max(ape::node.depth.edgelength(tree)[seq_along(tree$tip.label)])
  args <- list(data = counts, tree = tree, metadata = md,
               hierarchy = ~ region)
  for (bad in list(0, -1, NA_real_, Inf, c(0.5, 0.8))) {
    expect_error(suppressMessages(do.call(hillpart,
                                          c(args, list(depth = bad)))),
                 "depth")
  }
  expect_error(suppressMessages(do.call(hillpart,
                                        c(args, list(depth = height + 0.1)))),
               "tree height")
  nonultra <- tree
  tip_edge <- which(nonultra$edge[, 2] == 1)
  nonultra$edge.length[tip_edge] <- nonultra$edge.length[tip_edge] / 2
  args$tree <- nonultra
  expect_error(suppressMessages(do.call(hillpart,
                                        c(args, list(depth = height / 2)))),
               "ultrametric")
  expect_error(suppressMessages(hillpart(counts, tree = tree,
                                         depth = height / 2)), "hierarchy")
  expect_error(suppressMessages(hillpart(counts, metadata = md,
                                         hierarchy = ~ region, depth = 0.5)),
               "phylogenetic")
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

test_that("a functional tau grid matches separate hierarchy calls", {
  counts <- gut_counts
  d <- traits2dist(gut_traits)
  md <- data.frame(region = rep(c("control", "treatment"), each = 6),
                   site = rep(letters[1:6], each = 2),
                   row.names = colnames(counts))
  tau_grid <- c(0.25, 0.5, 1)
  args <- list(data = counts, dist = d, metadata = md,
               hierarchy = ~ region / site, q = c(0, 1, 2),
               weights = seq_len(ncol(counts)))
  grid <- suppressMessages(do.call(hillpart, c(args, list(tau = tau_grid))))
  expect_s3_class(grid, "hill_hierarchy")
  expect_true("tau" %in% names(grid))
  expect_equal(unique(grid$tau), tau_grid)
  wide <- suppressMessages(do.call(hillpart, c(
    args, list(tau = tau_grid, out = "matrix")
  )))
  expect_named(wide, paste0("tau_", tau_grid))
  for (i in seq_along(tau_grid)) {
    single <- suppressMessages(do.call(hillpart, c(
      args, list(tau = tau_grid[i])
    )))
    expect_false("tau" %in% names(single))
    rows <- grid[grid$tau == tau_grid[i], names(single)]
    rownames(rows) <- NULL
    expect_equal(rows$diversity, single$diversity, tolerance = 1e-12)
    expect_equal(rows$beta, single$beta, tolerance = 1e-12)
    one_wide <- suppressMessages(do.call(hillpart, c(
      args, list(tau = tau_grid[i], out = "matrix")
    )))
    expect_identical(wide[[i]], one_wide)
  }
  pdf_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file = pdf_file)
  on.exit({
    grDevices::dev.off()
    unlink(pdf_file)
  }, add = TRUE)
  expect_silent(plot(grid))
})

test_that("tau grids require functional hierarchy partitioning", {
  expect_error(hier(hierarchy = ~ region, tau = c(0.2, 0.5)),
               "functional")
  expect_error(suppressMessages(hillpart(gut_counts, tree = gut_tree,
                                         hierarchy = ~ region,
                                         metadata = data.frame(
                                           region = rep(1:2, each = 6)),
                                         tau = c(0.2, 0.5))), "functional")
  d <- traits2dist(gut_traits)
  md <- data.frame(region = rep(1:2, each = 6),
                   row.names = colnames(gut_counts))
  expect_error(suppressMessages(hillpart(gut_counts, dist = d,
                                         hierarchy = ~ region, metadata = md,
                                         tau = c(0.2, -0.5))), "tau")
})
