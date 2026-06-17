# Generates the bundled example data: `gut_counts`, `gut_tree`, `gut_traits`.
# These are SIMULATED data, designed to exercise every hilldiv3 path (neutral,
# phylogenetic and functional) on a small, realistic-looking gut microbiome
# MAG table. Run with: source("data-raw/make-data.R")
#
# Re-running regenerates data/*.rda; the fixed seed keeps them reproducible.

set.seed(2024)

n_mags <- 24L
n_samples <- 12L
mags <- sprintf("mag%02d", seq_len(n_mags))
# Two host groups (e.g. diet treatments), 6 samples each.
samples <- c(sprintf("ctrl%02d", 1:6), sprintf("trt%02d", 1:6))
groups <- rep(c("control", "treatment"), each = 6)

# --- Phylogeny over the MAGs -------------------------------------------------
gut_tree <- ape::rcoal(n_mags, tip.label = mags)
gut_tree$edge.length <- gut_tree$edge.length / max(
  ape::node.depth.edgelength(gut_tree)
)

# --- Count table (MAGs x samples) --------------------------------------------
# Group-dependent mean abundances so beta diversity is non-trivial: a block of
# MAGs is enriched in the treatment group, plus a shared lognormal background.
base <- rlnorm(n_mags, meanlog = 4, sdlog = 1.2)
enriched <- sample(seq_len(n_mags), 8L)
lambda <- matrix(base, nrow = n_mags, ncol = n_samples)
lambda[enriched, groups == "treatment"] <-
  lambda[enriched, groups == "treatment"] * 4

gut_counts <- matrix(rpois(n_mags * n_samples, lambda),
                     nrow = n_mags, dimnames = list(mags, samples))
# Sparsify: drop low-abundance occurrences to create realistic zeros.
gut_counts[gut_counts < 3] <- 0L
storage.mode(gut_counts) <- "integer"

# --- Functional traits (MAGs x traits) ---------------------------------------
gut_traits <- data.frame(
  genome_size   = round(rnorm(n_mags, 3.2, 0.6), 2),     # Mbp
  gc_content    = round(runif(n_mags, 0.35, 0.65), 3),   # proportion
  oxygen        = factor(sample(c("aerobe", "anaerobe", "facultative"),
                                n_mags, replace = TRUE)),
  motility      = sample(c(0L, 1L), n_mags, replace = TRUE),
  row.names = mags
)

usethis::use_data(gut_counts, gut_tree, gut_traits,
                  overwrite = TRUE, compress = "xz")
