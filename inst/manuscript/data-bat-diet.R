# Simulated dietary-metabarcoding data set for hilldiv3 use case 1.
#
# Insectivorous bats sampled non-invasively (guano), with insect prey
# characterised by COI metabarcoding and collapsed into amplicon sequence
# variants (ASVs). The design is nested: individual bats within roosts within
# two contrasting foraging habitats (forest vs. farmland).
#
# Forest bats are dietary SPECIALISTS that concentrate on a single, dominant
# moth (Lepidoptera) clade; farmland bats are GENERALISTS spread evenly across
# Diptera, Coleoptera and Hemiptera. Prey ASVs are organised into four insect
# orders, and the phylogeny has shallow within-order clades on a deep
# order-level backbone -- so a broad (multi-order) diet is phylogenetically far
# more diverse than a narrow (single-clade) one, and neutral vs. phylogenetic
# Hill numbers tell genuinely different stories.
#
# These data are SIMULATED purely to exercise the hilldiv3 paths; they are not
# real diet records. Sourcing this script defines `bat_counts` (ASVs x samples),
# `bat_tree` (prey phylogeny) and `bat_metadata` (sample design) in the caller.

suppressPackageStartupMessages(library(ape))
set.seed(2025)

# --- Prey ASVs, grouped into four insect orders -----------------------------
orders  <- c("Lepidoptera", "Diptera", "Coleoptera", "Hemiptera")
n_per   <- c(18L, 16L, 10L, 6L)            # 50 prey ASVs in total
asv     <- sprintf("asv%03d", seq_len(sum(n_per)))
asv_ord <- rep(orders, n_per)
names(asv_ord) <- asv

# --- Prey phylogeny: shallow per-order clades on a deep backbone ------------
scale_depth <- function(tr, depth) {
  tr$edge.length <- tr$edge.length /
    max(node.depth.edgelength(tr)) * depth
  tr
}
backbone <- scale_depth(rcoal(length(orders)), 1.0)   # deep order-level splits
backbone$tip.label <- orders

bat_tree <- backbone
for (o in orders) {
  tips  <- asv[asv_ord == o]
  clade <- scale_depth(rcoal(length(tips)), 0.15)     # shallow within order
  clade$tip.label <- tips
  bat_tree <- bind.tree(bat_tree, clade,
                        where = which(bat_tree$tip.label == o))
}
bat_tree <- multi2di(bat_tree)
bat_tree$edge.length[bat_tree$edge.length <= 0] <- 1e-6
bat_tree <- scale_depth(bat_tree, 1.0)                # unit total depth

# --- Nested sampling design -------------------------------------------------
# 2 habitats x 3 roosts x 5 individual bats = 30 guano samples.
bat_metadata <- data.frame(
  habitat = rep(c("forest", "farmland"), each = 15L),
  roost   = rep(sprintf("roost%d", 1:6), each = 5L),
  stringsAsFactors = FALSE
)
bat_metadata$bat <- sprintf("bat%02d", seq_len(nrow(bat_metadata)))
rownames(bat_metadata) <- bat_metadata$bat
forest <- bat_metadata$habitat == "forest"

# --- Prey detection and abundance -------------------------------------------
# Diet is zero-inflated: each prey ASV occurs in a sample with an
# order-by-habitat probability, then contributes Poisson read counts. A steep
# abundance distribution within Lepidoptera makes the forest specialists'
# diet dominance-skewed (a few moths dominate), depressing their evenness.
occ_forest   <- c(Lepidoptera = 0.85, Diptera = 0.06,
                  Coleoptera = 0.18, Hemiptera = 0.06)
occ_farmland <- c(Lepidoptera = 0.10, Diptera = 0.60,
                  Coleoptera = 0.50, Hemiptera = 0.60)

abund <- rlnorm(length(asv), meanlog = 3,
                sdlog = ifelse(asv_ord == "Lepidoptera", 1.6, 0.8))

bat_counts <- matrix(0L, nrow = length(asv), ncol = nrow(bat_metadata),
                     dimnames = list(asv, bat_metadata$bat))
for (j in seq_len(nrow(bat_metadata))) {
  p <- if (forest[j]) occ_forest[asv_ord] else occ_farmland[asv_ord]
  roost_effect <- exp(rnorm(1, 0, 0.3))               # roost-level variation
  present <- rbinom(length(asv), 1L, p)
  bat_counts[, j] <- as.integer(present * rpois(length(asv), abund * roost_effect))
}
bat_counts[bat_counts < 2L] <- 0L            # drop singletons -> realistic zeros

# Keep only ASVs detected in at least one sample, and align the tree.
keep <- rowSums(bat_counts) > 0
bat_counts <- bat_counts[keep, , drop = FALSE]
bat_tree   <- keep.tip(bat_tree, rownames(bat_counts))

invisible(list(counts = bat_counts, tree = bat_tree, metadata = bat_metadata))
