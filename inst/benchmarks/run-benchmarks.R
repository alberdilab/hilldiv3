#!/usr/bin/env Rscript

# Reproducible performance comparisons for hilldiv3 and related packages.
# Run from the package root with:
# Rscript inst/benchmarks/run-benchmarks.R

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || is.na(x)) y else x

args <- commandArgs(trailingOnly = TRUE)
repo_root_arg <- if (length(args) >= 1L) args[[1]] else getwd()
repo_root <- normalizePath(repo_root_arg, mustWork = TRUE)
out_dir <- normalizePath(
  Sys.getenv("BENCH_OUTPUT_DIR", file.path(repo_root, "inst/benchmarks/results")),
  mustWork = FALSE
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

iterations <- as.integer(Sys.getenv("BENCH_ITERATIONS", "5"))
n_taxa <- as.integer(Sys.getenv("BENCH_N_TAXA", "200"))
n_samples <- as.integer(Sys.getenv("BENCH_N_SAMPLES", "40"))
seed <- as.integer(Sys.getenv("BENCH_SEED", "20260617"))
q_all <- c(0, 1, 2)
q_one <- 1

load_hilldiv3 <- function(root) {
  if (basename(root) == "hilldiv3" &&
      file.exists(file.path(root, "DESCRIPTION")) &&
      requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(root, quiet = TRUE)
  } else {
    library(hilldiv3)
  }
  asNamespace("hilldiv3")
}

hd3 <- load_hilldiv3(repo_root)

make_data <- function(n_taxa, n_samples, seed) {
  set.seed(seed)
  taxa <- sprintf("sp%04d", seq_len(n_taxa))
  samples <- sprintf("sample%03d", seq_len(n_samples))
  base_abundance <- rgamma(n_taxa, shape = 0.45, rate = 0.08)
  counts <- vapply(seq_len(n_samples), function(i) {
    mu <- base_abundance * runif(n_taxa, 0.25, 2.5)
    stats::rpois(n_taxa, lambda = mu)
  }, integer(n_taxa))
  counts[rowSums(counts) == 0L, sample.int(n_samples, 1L)] <- 1L
  dimnames(counts) <- list(taxa, samples)

  tree <- ape::rcoal(n_taxa)
  tree$tip.label <- taxa

  traits <- data.frame(
    trait1 = stats::rnorm(n_taxa),
    trait2 = stats::runif(n_taxa),
    trait3 = stats::rpois(n_taxa, lambda = 4),
    row.names = taxa
  )
  distance <- stats::dist(traits)
  z <- exp(-as.matrix(distance) / max(distance))
  diag(z) <- 1

  metadata <- data.frame(
    region = rep(c("north", "south"), length.out = n_samples),
    site = rep(sprintf("site%02d", seq_len(max(2, ceiling(n_samples / 4)))),
               each = 4, length.out = n_samples),
    row.names = samples
  )

  list(
    counts = counts,
    comm = t(counts),
    tree = tree,
    traits = traits,
    dist = distance,
    z = z,
    metadata = metadata
  )
}

dat <- make_data(n_taxa, n_samples, seed)

task_table <- data.frame(
  task = c(
    "alpha_neutral",
    "alpha_phylogenetic",
    "alpha_functional",
    "partition_neutral",
    "partition_phylogenetic",
    "partition_functional",
    "overall_dissimilarity_neutral",
    "pairwise_neutral",
    "pairwise_phylogenetic",
    "pairwise_functional",
    "hierarchical_partition_neutral"
  ),
  operation = c(
    "Alpha Hill numbers",
    "Alpha phylogenetic Hill numbers",
    "Alpha functional Hill numbers",
    "Alpha/gamma/beta partitioning",
    "Phylogenetic alpha/gamma/beta partitioning",
    "Functional alpha/gamma/beta partitioning",
    "Overall Hill-number dissimilarity or similarity",
    "Pairwise Hill-number dissimilarity",
    "Pairwise phylogenetic Hill-number dissimilarity",
    "Pairwise functional Hill-number dissimilarity",
    "Nested hierarchical Hill-number partitioning"
  ),
  facet = c(
    "neutral", "phylogenetic", "functional",
    "neutral", "phylogenetic", "functional",
    "neutral", "neutral", "phylogenetic", "functional",
    "neutral"
  ),
  stringsAsFactors = FALSE
)

pkg_installed <- function(pkg) requireNamespace(pkg, quietly = TRUE)

call_if <- function(pkg, support, note, fun) {
  if (!pkg_installed(pkg)) {
    return(list(support = "unavailable", note = paste0("Package ", pkg, " is not installed."), fun = NULL))
  }
  list(support = support, note = note, fun = fun)
}

hd2 <- function(support = "yes", note = "Direct hilldiv2 call.") {
  call_if("hilldiv2", support, note, NULL)
}

adapters <- list(
  hilldiv3 = list(
    alpha_neutral = list(support = "yes", note = "Direct hilldiv3 call.",
                         fun = function() hd3$hilldiv(dat$counts, q = q_all, out = "matrix")),
    alpha_phylogenetic = list(support = "yes", note = "Direct hilldiv3 call.",
                              fun = function() hd3$hilldiv(dat$counts, q = q_all, tree = dat$tree, out = "matrix")),
    alpha_functional = list(support = "yes", note = "Direct hilldiv3 call.",
                            fun = function() hd3$hilldiv(dat$counts, q = q_all, dist = dat$dist, out = "matrix")),
    partition_neutral = list(support = "yes", note = "Direct hilldiv3 call.",
                             fun = function() hd3$hillpart(dat$counts, q = q_all, out = "matrix")),
    partition_phylogenetic = list(support = "yes", note = "Direct hilldiv3 call.",
                                  fun = function() hd3$hillpart(dat$counts, q = q_all, tree = dat$tree, out = "matrix")),
    partition_functional = list(support = "yes", note = "Direct hilldiv3 call.",
                                fun = function() hd3$hillpart(dat$counts, q = q_all, dist = dat$dist, out = "matrix")),
    overall_dissimilarity_neutral = list(support = "yes", note = "Direct hilldiv3 call.",
                                         fun = function() hd3$hilldiss(dat$counts, q = q_one, metric = "C", out = "matrix")),
    pairwise_neutral = list(support = "yes", note = "Direct hilldiv3 call.",
                            fun = function() hd3$hillpair(dat$counts, q = q_one, metric = "C")),
    pairwise_phylogenetic = list(support = "yes", note = "Direct hilldiv3 call.",
                                 fun = function() hd3$hillpair(dat$counts, q = q_one, metric = "C", tree = dat$tree)),
    pairwise_functional = list(support = "yes", note = "Direct hilldiv3 call.",
                               fun = function() hd3$hillpair(dat$counts, q = q_one, metric = "C", dist = dat$dist)),
    hierarchical_partition_neutral = list(support = "yes", note = "Direct hilldiv3 call.",
                                          fun = function() hd3$hillpart(dat$counts, q = q_one,
                                                                       hierarchy = ~ region / site,
                                                                       metadata = dat$metadata))
  ),
  hilldiv2 = list(
    alpha_neutral = call_if("hilldiv2", "yes", "Direct hilldiv2 call.",
                            function() hilldiv2::hilldiv(data = dat$counts, q = q_all)),
    alpha_phylogenetic = call_if("hilldiv2", "yes", "Direct hilldiv2 call.",
                                 function() hilldiv2::hilldiv(data = dat$counts, q = q_all, tree = dat$tree)),
    alpha_functional = call_if("hilldiv2", "yes", "Direct hilldiv2 call.",
                               function() hilldiv2::hilldiv(data = dat$counts, q = q_all, dist = as.matrix(dat$dist))),
    partition_neutral = call_if("hilldiv2", "yes", "Direct hilldiv2 call.",
                                function() hilldiv2::hillpart(data = dat$counts, q = q_all)),
    partition_phylogenetic = call_if("hilldiv2", "yes", "Direct hilldiv2 call.",
                                     function() hilldiv2::hillpart(data = dat$counts, q = q_all, tree = dat$tree)),
    partition_functional = call_if("hilldiv2", "yes", "Direct hilldiv2 call.",
                                   function() hilldiv2::hillpart(data = dat$counts, q = q_all, dist = as.matrix(dat$dist))),
    overall_dissimilarity_neutral = call_if("hilldiv2", "yes", "Direct hilldiv2 call.",
                                            function() hilldiv2::hilldiss(data = dat$counts, q = q_one, metric = "C")),
    pairwise_neutral = call_if("hilldiv2", "yes", "Direct hilldiv2 call.",
                               function() hilldiv2::hillpair(data = dat$counts, q = q_one, metric = "C")),
    pairwise_phylogenetic = call_if("hilldiv2", "yes", "Direct hilldiv2 call.",
                                    function() hilldiv2::hillpair(data = dat$counts, q = q_one, metric = "C", tree = dat$tree)),
    pairwise_functional = call_if("hilldiv2", "yes", "Direct hilldiv2 call.",
                                  function() hilldiv2::hillpair(data = dat$counts, q = q_one, metric = "C", dist = as.matrix(dat$dist))),
    hierarchical_partition_neutral = list(support = "no", note = "No nested hierarchy argument in hilldiv2.", fun = NULL)
  ),
  hillR = list(
    alpha_neutral = call_if("hillR", "yes", "Taxonomic Hill numbers via hill_taxa().",
                            function() lapply(q_all, function(q) hillR::hill_taxa(dat$comm, q = q))),
    alpha_phylogenetic = call_if("hillR", "yes", "Phylogenetic Hill numbers via hill_phylo().",
                                 function() lapply(q_all, function(q) hillR::hill_phylo(dat$comm, dat$tree, q = q, show_warning = FALSE))),
    alpha_functional = call_if("hillR", "yes", "Functional Hill numbers via hill_func().",
                               function() lapply(q_all, function(q) hillR::hill_func(dat$comm, dat$traits, q = q, check_data = FALSE))),
    partition_neutral = call_if("hillR", "yes", "Partitioning via hill_taxa_parti().",
                                function() lapply(q_all, function(q) hillR::hill_taxa_parti(dat$comm, q = q, show_warning = FALSE))),
    partition_phylogenetic = call_if("hillR", "yes", "Partitioning via hill_phylo_parti().",
                                     function() lapply(q_all, function(q) hillR::hill_phylo_parti(dat$comm, dat$tree, q = q, show_warning = FALSE))),
    partition_functional = call_if("hillR", "yes", "Partitioning via hill_func_parti().",
                                   function() lapply(q_all, function(q) hillR::hill_func_parti(dat$comm, dat$traits, q = q, show_warning = FALSE))),
    overall_dissimilarity_neutral = call_if("hillR", "partial", "Returns beta and local/regional similarity, not hilldiv3's four S/C/U/V dissimilarities.",
                                            function() hillR::hill_taxa_parti(dat$comm, q = q_one, show_warning = FALSE)),
    pairwise_neutral = call_if("hillR", "yes", "Pairwise partitioning via hill_taxa_parti_pairwise().",
                               function() hillR::hill_taxa_parti_pairwise(dat$comm, q = q_one, .progress = FALSE, show_warning = FALSE)),
    pairwise_phylogenetic = call_if("hillR", "yes", "Pairwise partitioning via hill_phylo_parti_pairwise().",
                                    function() hillR::hill_phylo_parti_pairwise(dat$comm, dat$tree, q = q_one, .progress = FALSE, show_warning = FALSE)),
    pairwise_functional = call_if("hillR", "yes", "Pairwise partitioning via hill_func_parti_pairwise().",
                                  function() hillR::hill_func_parti_pairwise(dat$comm, dat$traits, q = q_one, .progress = FALSE, show_warning = FALSE)),
    hierarchical_partition_neutral = list(support = "no", note = "No nested hierarchy argument in hillR.", fun = NULL)
  ),
  entropart = list(
    alpha_neutral = call_if("entropart", "yes", "Neutral alpha diversity via AlphaDiversity().",
                            function() {
                              mc <- entropart::MetaCommunity(as.data.frame(dat$counts))
                              lapply(q_all, function(q) entropart::AlphaDiversity(mc, q = q, Correction = "None"))
                            }),
    alpha_phylogenetic = call_if("entropart", "partial", "Phylogenetic alpha requires an ultrametric tree and returns entropart's normalized phylodiversity.",
                                 function() {
                                   mc <- entropart::MetaCommunity(as.data.frame(dat$counts))
                                   lapply(q_all, function(q) entropart::AlphaDiversity(mc, q = q, Correction = "None", Tree = dat$tree))
                                 }),
    alpha_functional = call_if("entropart", "partial", "Similarity-based alpha uses Z, not hilldiv3's distance-threshold functional definition.",
                               function() {
                                 mc <- entropart::MetaCommunity(as.data.frame(dat$counts))
                                 lapply(q_all, function(q) entropart::AlphaDiversity(mc, q = q, Correction = "None", Z = dat$z))
                               }),
    partition_neutral = call_if("entropart", "yes", "Metacommunity partitioning via DivPart().",
                                function() {
                                  mc <- entropart::MetaCommunity(as.data.frame(dat$counts))
                                  lapply(q_all, function(q) entropart::DivPart(q = q, MC = mc, Biased = TRUE, Correction = "None"))
                                }),
    partition_phylogenetic = call_if("entropart", "partial", "Phylogenetic partitioning requires an ultrametric tree and uses entropart's normalization.",
                                     function() {
                                       mc <- entropart::MetaCommunity(as.data.frame(dat$counts))
                                       lapply(q_all, function(q) entropart::DivPart(q = q, MC = mc, Biased = TRUE, Correction = "None", Tree = dat$tree))
                                     }),
    partition_functional = call_if("entropart", "partial", "Similarity-based partitioning uses Z, not hilldiv3's distance-threshold functional definition.",
                                   function() {
                                     mc <- entropart::MetaCommunity(as.data.frame(dat$counts))
                                     lapply(q_all, function(q) entropart::DivPart(q = q, MC = mc, Biased = TRUE, Correction = "None", Z = dat$z))
                                   }),
    overall_dissimilarity_neutral = call_if("entropart", "partial", "Returns beta diversity, not hilldiv3's four S/C/U/V dissimilarities.",
                                            function() {
                                              mc <- entropart::MetaCommunity(as.data.frame(dat$counts))
                                              entropart::DivPart(q = q_one, MC = mc, Biased = TRUE, Correction = "None")
                                            }),
    pairwise_neutral = list(support = "no", note = "No direct all-pairs Hill dissimilarity API.", fun = NULL),
    pairwise_phylogenetic = list(support = "no", note = "No direct all-pairs Hill dissimilarity API.", fun = NULL),
    pairwise_functional = list(support = "no", note = "No direct all-pairs Hill dissimilarity API.", fun = NULL),
    hierarchical_partition_neutral = call_if("entropart", "partial", "MergeMC supports hierarchical metacommunities, but not the same formula API/output.",
                                             function() {
                                               groups <- split(colnames(dat$counts), dat$metadata$region)
                                               mcs <- lapply(groups, function(cols) entropart::MetaCommunity(as.data.frame(dat$counts[, cols, drop = FALSE])))
                                               merged <- entropart::MergeMC(mcs)
                                               entropart::DivPart(q = q_one, MC = merged, Biased = TRUE, Correction = "None")
                                             })
  ),
  vegan = list(
    alpha_neutral = call_if("vegan", "yes", "Neutral Hill numbers via renyi(..., hill = TRUE).",
                            function() vegan::renyi(dat$comm, scales = q_all, hill = TRUE)),
    alpha_phylogenetic = list(support = "no", note = "vegan does not compute phylogenetic Hill numbers.", fun = NULL),
    alpha_functional = list(support = "no", note = "vegan does not compute functional Hill numbers.", fun = NULL),
    partition_neutral = list(support = "no", note = "No Hill alpha/gamma/beta partitioning API.", fun = NULL),
    partition_phylogenetic = list(support = "no", note = "No phylogenetic Hill partitioning API.", fun = NULL),
    partition_functional = list(support = "no", note = "No functional Hill partitioning API.", fun = NULL),
    overall_dissimilarity_neutral = list(support = "no", note = "vegdist computes pairwise ecological distances, not overall Hill S/C/U/V metrics.", fun = NULL),
    pairwise_neutral = call_if("vegan", "partial", "Pairwise Bray-Curtis distance; related, but not Hill-number dissimilarity.",
                               function() vegan::vegdist(dat$comm, method = "bray")),
    pairwise_phylogenetic = list(support = "no", note = "No phylogenetic Hill pairwise dissimilarity API.", fun = NULL),
    pairwise_functional = list(support = "no", note = "No functional Hill pairwise dissimilarity API.", fun = NULL),
    hierarchical_partition_neutral = list(support = "no", note = "No nested Hill partitioning API.", fun = NULL)
  ),
  BAT = list(
    alpha_neutral = call_if("BAT", "yes", "Neutral Hill numbers via hill().",
                            function() lapply(q_all, function(q) BAT::hill(dat$comm, q = q))),
    alpha_phylogenetic = call_if("BAT", "partial", "BAT::alpha() computes Faith-style PD richness, not q = 0, 1, 2 phylogenetic Hill numbers.",
                                 function() BAT::alpha(dat$comm, dat$tree)),
    alpha_functional = call_if("BAT", "partial", "BAT::alpha() computes tree-based functional richness, not q = 0, 1, 2 functional Hill numbers.",
                               function() BAT::alpha(dat$comm, dat$traits)),
    partition_neutral = list(support = "no", note = "No Hill alpha/gamma/beta partitioning API.", fun = NULL),
    partition_phylogenetic = list(support = "no", note = "No phylogenetic Hill alpha/gamma/beta partitioning API.", fun = NULL),
    partition_functional = list(support = "no", note = "No functional Hill alpha/gamma/beta partitioning API.", fun = NULL),
    overall_dissimilarity_neutral = list(support = "no", note = "BAT beta functions are pairwise replacement/richness components, not overall Hill S/C/U/V metrics.", fun = NULL),
    pairwise_neutral = call_if("BAT", "partial", "Pairwise Jaccard/Sorensen beta components, not Hill-number dissimilarity.",
                               function() BAT::beta(dat$comm, func = "jaccard", abund = TRUE)),
    pairwise_phylogenetic = call_if("BAT", "partial", "Pairwise PD beta components, not Hill-number dissimilarity.",
                                    function() BAT::beta(dat$comm, dat$tree, func = "jaccard", abund = TRUE)),
    pairwise_functional = call_if("BAT", "partial", "Pairwise FD beta components, not Hill-number dissimilarity.",
                                  function() BAT::beta(dat$comm, dat$traits, func = "jaccard", abund = TRUE)),
    hierarchical_partition_neutral = list(support = "no", note = "No nested Hill partitioning API.", fun = NULL)
  )
)

quiet_run <- function(fun) {
  tmp <- tempfile()
  con <- file(tmp, open = "wt")
  sink(con, type = "output")
  on.exit({
    sink(type = "output")
    close(con)
    unlink(tmp)
  }, add = TRUE)
  suppressMessages(suppressWarnings(fun()))
}

object_bytes <- function(x) as.numeric(utils::object.size(x))

measure_base <- function(fun, iterations) {
  elapsed <- numeric(iterations)
  result <- NULL
  for (i in seq_len(iterations)) {
    gc()
    start <- proc.time()[["elapsed"]]
    result <- quiet_run(fun)
    elapsed[[i]] <- proc.time()[["elapsed"]] - start
  }
  list(
    median_seconds = stats::median(elapsed),
    min_seconds = min(elapsed),
    max_seconds = max(elapsed),
    memory_bytes = NA_real_,
    result_bytes = object_bytes(result),
    backend = "system.time"
  )
}

measure_bench <- function(fun, iterations) {
  mark <- tryCatch(
    bench::mark(quiet_run(fun), iterations = iterations, check = FALSE),
    error = function(e) e
  )
  memory_profiled <- TRUE
  if (inherits(mark, "error") &&
      grepl("Memory profiling failed", conditionMessage(mark), fixed = TRUE)) {
    mark <- bench::mark(
      quiet_run(fun),
      iterations = iterations,
      check = FALSE,
      memory = FALSE
    )
    memory_profiled <- FALSE
  }
  if (inherits(mark, "error")) {
    stop(mark)
  }
  result <- quiet_run(fun)
  list(
    median_seconds = as.numeric(mark$median[[1]], units = "secs"),
    min_seconds = as.numeric(min(mark$time[[1]]), units = "secs"),
    max_seconds = as.numeric(max(mark$time[[1]]), units = "secs"),
    memory_bytes = if (memory_profiled) {
      as.numeric(mark$mem_alloc[[1]], units = "bytes")
    } else {
      NA_real_
    },
    result_bytes = object_bytes(result),
    backend = "bench"
  )
}

measure <- function(fun, iterations) {
  if (pkg_installed("bench")) measure_bench(fun, iterations) else measure_base(fun, iterations)
}

record_one <- function(package, task, adapter) {
  task_meta <- task_table[task_table$task == task, , drop = FALSE]
  base <- data.frame(
    package = package,
    package_version = if (pkg_installed(package)) as.character(utils::packageVersion(package)) else NA_character_,
    task = task,
    operation = task_meta$operation,
    facet = task_meta$facet,
    support = adapter$support,
    median_seconds = NA_real_,
    min_seconds = NA_real_,
    max_seconds = NA_real_,
    memory_bytes = NA_real_,
    result_bytes = NA_real_,
    backend = NA_character_,
    iterations = iterations,
    n_taxa = n_taxa,
    n_samples = n_samples,
    note = adapter$note,
    error = NA_character_,
    stringsAsFactors = FALSE
  )

  if (is.null(adapter$fun) || !(adapter$support %in% c("yes", "partial"))) {
    return(base)
  }

  measured <- tryCatch(measure(adapter$fun, iterations), error = function(e) e)
  if (inherits(measured, "error")) {
    base$support <- "failed"
    base$error <- conditionMessage(measured)
    return(base)
  }

  base$median_seconds <- measured$median_seconds
  base$min_seconds <- measured$min_seconds
  base$max_seconds <- measured$max_seconds
  base$memory_bytes <- measured$memory_bytes
  base$result_bytes <- measured$result_bytes
  base$backend <- measured$backend
  base
}

results <- do.call(rbind, unlist(lapply(names(adapters), function(pkg) {
  lapply(task_table$task, function(task) {
    record_one(pkg, task, adapters[[pkg]][[task]])
  })
}), recursive = FALSE))

results$relative_to_hilldiv3 <- NA_real_
for (task in unique(results$task)) {
  baseline <- results$median_seconds[results$package == "hilldiv3" & results$task == task]
  if (length(baseline) == 1L && !is.na(baseline) && baseline > 0) {
    idx <- results$task == task & !is.na(results$median_seconds)
    results$relative_to_hilldiv3[idx] <- results$median_seconds[idx] / baseline
  }
}

session <- utils::capture.output(utils::sessionInfo())

write.csv(results, file.path(out_dir, "performance.csv"), row.names = FALSE)
writeLines(session, file.path(out_dir, "session-info.txt"))

message("Wrote benchmark results to: ", file.path(out_dir, "performance.csv"))
