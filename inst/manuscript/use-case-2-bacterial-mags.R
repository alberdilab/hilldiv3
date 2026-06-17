# hilldiv3 manuscript -- Use case 2: gut bacterial MAGs under an intervention
#
# Reproduces every figure and number cited in use-case-2-bacterial-mags.md. Run
# from the package root with:  Rscript inst/manuscript/use-case-2-bacterial-mags.R
#
# Uses the bundled simulated data: gut_counts (24 MAGs x 12 host gut samples,
# control vs. treatment), gut_tree (genome phylogeny) and gut_traits (mixed
# genomic/physiological traits). The intervention swaps which of two deep clades
# dominates; the clades are functional mirrors, so the community turns over
# neutrally and PHYLOGENETICALLY while its FUNCTIONAL make-up is conserved.
# Showcased capabilities: neutral, phylogenetic AND functional diversity from one
# call (hilldiv); trait-to-distance conversion (traits2dist); between-group
# partitioning per flavour (hillpart); pairwise dissimilarity + ordination in all
# three spaces (hillpair); and functional redundancy (hillred).

suppressPackageStartupMessages({
  library(hilldiv3)
  library(ggplot2)
})

here   <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))
if (length(here) == 0) here <- "inst/manuscript"
figdir <- file.path(here, "figures")
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

# long-format hilldiv3 result -> plain data.frame (avoids wide pivot / coercion)
as_df <- function(x) { class(x) <- "data.frame"; x }

# --- Data and design --------------------------------------------------------
group <- rep(c("control", "treatment"), each = 6)
names(group) <- colnames(gut_counts)
metadata <- data.frame(group = group, row.names = colnames(gut_counts))
grp <- function(s) group[s]
pal <- c(control = "#1F77B4", treatment = "#D62728")

# Functional distance from the mixed trait table (Gower handles the mix of
# continuous, categorical and binary traits).
fdist <- traits2dist(gut_traits)

theme_set(theme_bw(base_size = 12) +
            theme(panel.grid.minor = element_blank(),
                  legend.position = "top"))

# --- 1. The same call, three flavours of diversity --------------------------
neu <- as_df(hilldiv(gut_counts, q = c(0, 1, 2)))
phy <- as_df(hilldiv(gut_counts, q = c(0, 1, 2), tree = gut_tree))
fun <- as_df(hilldiv(gut_counts, q = c(0, 1, 2), dist = fdist))
neu$flavour <- "Neutral (qD)"
phy$flavour <- "Phylogenetic (qPD)"
fun$flavour <- "Functional (qFD)"
alpha <- rbind(neu, phy, fun)
alpha$group   <- grp(alpha$sample)
alpha$flavour <- factor(alpha$flavour,
                        c("Neutral (qD)", "Phylogenetic (qPD)", "Functional (qFD)"))

p_alpha <- ggplot(alpha, aes(factor(q), value, fill = group)) +
  geom_boxplot(outlier.size = 0.6, width = 0.7) +
  facet_wrap(~ flavour, scales = "free_y") +
  scale_fill_manual(values = pal, name = NULL) +
  labs(x = "Diversity order (q)", y = "Effective number of lineages")
ggsave(file.path(figdir, "uc2-alpha.png"), p_alpha, width = 8, height = 3.4, dpi = 200)

# --- 2. Partitioning between groups, three flavours -------------------------
part <- function(type, ...) {
  m <- hillpart(gut_counts, q = c(0, 1, 2), hierarchy = ~ group,
                metadata = metadata, ...)
  m <- as_df(m)[as_df(m)$scale == "total", c("q", "beta")]
  m$flavour <- type
  m
}
beta <- rbind(part("Neutral"),
              part("Phylogenetic", tree = gut_tree),
              part("Functional",   dist = fdist))
beta$flavour <- factor(beta$flavour, c("Neutral", "Phylogenetic", "Functional"))

p_beta <- ggplot(beta, aes(factor(q), beta, fill = flavour)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey40") +
  scale_fill_brewer(palette = "Set2", name = NULL) +
  labs(x = "Diversity order (q)",
       y = expression("Between-group turnover (" * beta * ")"))
ggsave(file.path(figdir, "uc2-partition.png"), p_beta, width = 5.5, height = 3.6, dpi = 200)

# --- 3. Ordination in neutral, phylogenetic and functional space ------------
ord_one <- function(d, lab) {
  pc <- cmdscale(d, k = 2, eig = TRUE)
  o  <- data.frame(pc$points, mag = rownames(pc$points))
  names(o)[1:2] <- c("PCoA1", "PCoA2")
  o$group <- grp(o$mag); o$space <- lab
  attr(o, "ev") <- round(100 * pc$eig[1:2] / sum(pc$eig[pc$eig > 0]), 1)
  o
}
on <- ord_one(hillpair(gut_counts, q = 1, metric = "C"),               "Neutral (q=1)")
op <- ord_one(hillpair(gut_counts, q = 1, metric = "C", tree = gut_tree), "Phylogenetic (q=1)")
of <- ord_one(hillpair(gut_counts, q = 1, metric = "C", dist = fdist),  "Functional (q=1)")
ord <- rbind(on, op, of)
ord$space <- factor(ord$space, c("Neutral (q=1)", "Phylogenetic (q=1)", "Functional (q=1)"))

p_ord <- ggplot(ord, aes(PCoA1, PCoA2, colour = group)) +
  geom_point(size = 2.4) +
  stat_ellipse(level = 0.68) +
  facet_wrap(~ space, scales = "free") +
  scale_colour_manual(values = pal, name = NULL) +
  labs(x = "PCoA 1", y = "PCoA 2")
ggsave(file.path(figdir, "uc2-ordination.png"), p_ord, width = 8.5, height = 3.4, dpi = 200)

# --- 4. Functional redundancy -----------------------------------------------
red <- hillred(gut_counts, q = c(1, 2), dist = fdist)
png(file.path(figdir, "uc2-redundancy.png"), width = 5, height = 4,
    units = "in", res = 200)
plot(red)
dev.off()

# --- Numbers cited in the prose ---------------------------------------------
cat("\n--- Use case 2 summary ---\n")
cat("MAGs x samples:", paste(dim(gut_counts), collapse = " x "), "\n")
cat("functional distance range:", paste(round(range(fdist), 3), collapse = "-"), "\n\n")
cat("Mean alpha by group:\n")
print(aggregate(value ~ flavour + q + group, alpha, function(x) round(mean(x), 2)))
cat("\nBetween-group beta (total scale):\n")
print(transform(beta, beta = round(beta, 3)))
cat("\nFunctional redundancy:\n");   print(as_df(red)[, c("q", "redundancy")])
cat("\nPCoA1 group centroids (neutral / phylogenetic / functional):\n")
print(aggregate(PCoA1 ~ group, on, function(x) round(mean(x), 3)))
print(aggregate(PCoA1 ~ group, op, function(x) round(mean(x), 3)))
print(aggregate(PCoA1 ~ group, of, function(x) round(mean(x), 3)))
cat("\nFigures written to", figdir, "\n")
