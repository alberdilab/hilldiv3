# hilldiv3 manuscript -- Use case 1: dietary metabarcoding of insectivorous bats
#
# Reproduces every figure and number cited in use-case-1-bat-diet.md. Run from
# the package root with:  Rscript inst/manuscript/use-case-1-bat-diet.R
#
# Showcased capabilities: neutral and phylogenetic alpha diversity (hilldiv),
# diversity profiles (hillprof), evenness (hilleven) and multi-scale nested
# partitioning (hillpart with a hierarchy formula).

suppressPackageStartupMessages({
  library(hilldiv3)
  library(ggplot2)
})

here   <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))
if (length(here) == 0) here <- "inst/manuscript"
figdir <- file.path(here, "figures")
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

# --- Data -------------------------------------------------------------------
sim          <- source(file.path(here, "data-bat-diet.R"))$value
bat_counts   <- sim$counts
bat_tree     <- sim$tree
bat_metadata <- sim$metadata
hab <- function(s) bat_metadata$habitat[match(s, bat_metadata$bat)]
pal <- c(forest = "#2E7D32", farmland = "#C9A227")

# hilldiv3 results are long-format data.frame subclasses; drop the class so
# base `$<-` / `rbind` work without coercion (as.data.frame() would pivot wide).
as_df <- function(x) { class(x) <- "data.frame"; x }

theme_set(theme_bw(base_size = 12) +
            theme(panel.grid.minor = element_blank(),
                  legend.position = "top"))

# --- 1. Neutral vs. phylogenetic alpha diversity ----------------------------
# One call with `tree =` returns neutral and phylogenetic diversity together,
# tagged by a `type` column.
labs <- c(neutral = "Neutral (qD)", phylogenetic = "Phylogenetic (qPD)")
alpha <- as_df(hilldiv(bat_counts, q = c(0, 1, 2), tree = bat_tree))
alpha$flavour <- factor(labs[alpha$type], labs)
alpha$habitat <- hab(alpha$sample)

p_alpha <- ggplot(alpha, aes(factor(q), value, fill = habitat)) +
  geom_boxplot(outlier.size = 0.6, width = 0.7) +
  facet_wrap(~ flavour, scales = "free_y") +
  scale_fill_manual(values = pal, name = NULL) +
  labs(x = "Diversity order (q)", y = "Effective number of prey lineages")
ggsave(file.path(figdir, "uc1-alpha.png"), p_alpha, width = 7, height = 3.4, dpi = 200)

# --- 2. Diversity profiles (neutral) by habitat -----------------------------
prof <- as_df(hillprof(bat_counts, q = seq(0, 3, by = 0.1)))
prof$habitat <- hab(prof$sample)
prof_mean <- aggregate(value ~ q + habitat, prof, mean)

p_prof <- ggplot(prof, aes(q, value, group = sample, colour = habitat)) +
  geom_line(alpha = 0.25, linewidth = 0.3) +
  geom_line(data = prof_mean, aes(group = habitat), linewidth = 1.3) +
  scale_colour_manual(values = pal, name = NULL) +
  labs(x = "Diversity order (q)", y = "Neutral diversity (qD)")
ggsave(file.path(figdir, "uc1-profile.png"), p_prof, width = 5, height = 3.6, dpi = 200)

# --- 3. Evenness ------------------------------------------------------------
even <- as_df(hilleven(bat_counts, q = 2))
even$habitat <- hab(even$sample)

# --- 4. Multi-scale (nested) partitioning -----------------------------------
# individuals < roosts < habitats < total, for neutral and phylogenetic.
hier_n <- hillpart(bat_counts, q = c(0, 1, 2),
                   hierarchy = ~ habitat / roost, metadata = bat_metadata)
hier_p <- hillpart(bat_counts, q = c(0, 1, 2), tree = bat_tree,
                   hierarchy = ~ habitat / roost, metadata = bat_metadata)

beta_n <- as_df(hier_n)[as_df(hier_n)$scale != "sample", ]
beta_p <- as_df(hier_p)[as_df(hier_p)$scale != "sample", ]
beta_n$flavour <- "Neutral"; beta_p$flavour <- "Phylogenetic"
betas <- rbind(beta_n, beta_p)
betas$scale <- factor(betas$scale, c("roost", "habitat", "total"),
                      labels = c("among individuals\n(within roost)",
                                 "among roosts\n(within habitat)",
                                 "among habitats"))

p_beta <- ggplot(betas, aes(scale, beta, fill = factor(q))) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  facet_wrap(~ flavour) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey40") +
  scale_fill_brewer(palette = "Blues", name = "q") +
  labs(x = NULL, y = expression("Turnover gained at scale (" * beta[k] * ")"))
ggsave(file.path(figdir, "uc1-partition.png"), p_beta, width = 7.5, height = 3.6, dpi = 200)

# --- Numbers cited in the prose ---------------------------------------------
cat("\n--- Use case 1 summary ---\n")
cat("ASVs x samples:", paste(dim(bat_counts), collapse = " x "), "\n")
cat("per-sample prey richness:", paste(range(colSums(bat_counts > 0)), collapse = "-"), "\n\n")
cat("Mean alpha by habitat:\n")
print(aggregate(value ~ flavour + q + habitat, alpha, function(x) round(mean(x), 2)))
cat("\nMean q2 evenness by habitat:\n")
print(aggregate(value ~ habitat, even, function(x) round(mean(x), 3)))
cat("\nNeutral multi-scale partition:\n");      print(hier_n)
cat("\nPhylogenetic multi-scale partition:\n"); print(hier_p)
cat("\nFigures written to", figdir, "\n")
