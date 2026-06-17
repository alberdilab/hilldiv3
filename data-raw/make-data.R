# Generates the bundled example data: `gut_counts`, `gut_tree`, `gut_traits`.
# These are SIMULATED data, designed to exercise every hilldiv3 path (neutral,
# phylogenetic and functional) on a small, realistic-looking gut microbiome
# MAG table. Run with: source("data-raw/make-data.R")
#
# Re-running regenerates data/*.rda; the fixed seed keeps them reproducible.
#
# Design (the use case is built on this structure): the 24 MAGs fall into two
# deep, reciprocally monophyletic clades of 12. Clade A dominates control
# samples, clade B dominates treatment samples -- so the intervention *swaps
# which lineages dominate* without changing richness. Crucially, clade B is a
# FUNCTIONAL MIRROR of clade A: the i-th genome of clade B carries the same trait
# profile as the i-th genome of clade A. Every functional role is therefore held
# by one control-clade genome AND one treatment-clade genome, so swapping which
# clade dominates leaves the community's FUNCTIONAL composition unchanged even as
# its neutral and PHYLOGENETIC composition turns over completely -- the signature
# of high functional redundancy, which `hillred()` then quantifies.

set.seed(22)

n_mags <- 24L
n_samples <- 12L
mags <- sprintf("mag%02d", seq_len(n_mags))
# Two host groups (control vs. an intervention), 6 samples each.
samples <- c(sprintf("ctrl%02d", 1:6), sprintf("trt%02d", 1:6))
groups <- rep(c("control", "treatment"), each = 6)

clade_A <- 1:12    # dominant in control
clade_B <- 13:24   # enriched (dominant) under treatment

# --- Phylogeny: two deep, reciprocally monophyletic clades -------------------
# Build a shallow coalescent for clade A, then make clade B an exact structural
# mirror of it (same within-clade topology and branch lengths, relabelled tips)
# and join the two under a deep split. The clades are thus genuinely different
# parts of the tree, yet internally identical -- so a sample dominated by either
# clade has the same within-sample (alpha) phylogenetic diversity, and all the
# between-group phylogenetic signal comes from the deep split they straddle.
sub_A <- ape::rcoal(length(clade_A), tip.label = mags[clade_A])
sub_A$edge.length <- sub_A$edge.length * 0.3
sub_B <- sub_A
sub_B$tip.label <- mags[clade_B]
gut_tree <- ape::read.tree(text = sprintf(
  "(%s:0.7,%s:0.7);",
  sub(";$", "", ape::write.tree(sub_A)),
  sub(";$", "", ape::write.tree(sub_B))
))
gut_tree$edge.length <- gut_tree$edge.length / max(
  ape::node.depth.edgelength(gut_tree)
)

# --- Functional traits: clade B mirrors clade A ------------------------------
# Twelve distinct functional profiles span the trait space; profile i is carried
# by genome i of clade A and, identically (up to tiny measurement noise), by
# genome i of clade B. Function is thus orthogonal to the deep phylogenetic
# split, and each profile is held redundantly across the two clades.
n_prof <- length(clade_A)                              # 12 profiles
profile <- integer(n_mags)
profile[clade_A] <- seq_len(n_prof)
profile[clade_B] <- seq_len(n_prof)

prof_traits <- data.frame(
  genome_size = seq(1.6, 6.0, length.out = n_prof),                # Mbp
  gc_content  = seq(0.33, 0.67, length.out = n_prof),              # proportion
  oxygen      = rep(c("anaerobe", "microaerophile", "facultative", "aerobe"),
                    length.out = n_prof),
  motility    = rep(c(0L, 1L), length.out = n_prof)
)
gut_traits <- data.frame(
  genome_size = round(prof_traits$genome_size[profile] + rnorm(n_mags, 0, 0.06), 2),
  gc_content  = round(pmin(pmax(prof_traits$gc_content[profile] +
                                  rnorm(n_mags, 0, 0.008), 0.30), 0.70), 3),
  oxygen      = factor(prof_traits$oxygen[profile]),
  motility    = prof_traits$motility[profile],
  row.names   = mags
)

# --- Count table (MAGs x samples) --------------------------------------------
# Shared lognormal background keeps every MAG present (richness ~ constant);
# the resident clade is boosted in each group so the *dominant* lineages -- and
# hence neutral and phylogenetic composition -- swap between control and
# treatment. Functional composition does not, because clade B mirrors clade A.
#
# Each sample also gets a random preference over the 12 functional profiles,
# drawn independently of group and applied to BOTH clade members of a profile.
# This makes per-sample neutral and functional diversity covary (giving the
# redundancy curve real signal) WITHOUT differing systematically between control
# and treatment (so functional beta and the functional ordination stay flat).
prof_pref <- matrix(exp(rnorm(n_prof * n_samples, 0, 0.85)), nrow = n_prof,
                    dimnames = list(NULL, samples))
pref <- prof_pref[profile, ]                          # n_mags x n_samples

base <- rlnorm(n_mags, meanlog = 3.0, sdlog = 0.5)
boost <- 7
lambda <- matrix(base, nrow = n_mags, ncol = n_samples,
                 dimnames = list(mags, samples)) * pref
lambda[clade_A, groups == "control"]   <- lambda[clade_A, groups == "control"]   * boost
lambda[clade_B, groups == "treatment"] <- lambda[clade_B, groups == "treatment"] * boost

gut_counts <- matrix(rpois(n_mags * n_samples, lambda),
                     nrow = n_mags, dimnames = list(mags, samples))
# Sparsify: drop low-abundance occurrences to create realistic zeros.
gut_counts[gut_counts < 3] <- 0L
storage.mode(gut_counts) <- "integer"

usethis::use_data(gut_counts, gut_tree, gut_traits,
                  overwrite = TRUE, compress = "xz")
