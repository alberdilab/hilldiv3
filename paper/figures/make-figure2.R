# Fig. 2 (benchmarks) for the hilldiv3 JOSS paper.
# Run from the package root:  Rscript paper/figures/make-figure2.R
# Reads inst/benchmarks/results/performance-summary.csv (read-only), draws three
# svglite panels at their final width, composes them with compose.py and exports
# PDF (text outlined) and PNG (600 dpi) with Inkscape.
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(patchwork); library(svglite)})

root <- normalizePath(getwd())
stopifnot(file.exists(file.path(root, "inst/benchmarks/results/performance-summary.csv")))
fig  <- file.path(root, "paper/figures")
for (d in c("raw", "src", "final")) dir.create(file.path(fig, d), FALSE, TRUE)
compose <- Sys.getenv("COMPOSE_PY", "paper/figures/compose.py")
ink     <- Sys.getenv("INKSCAPE", "/Applications/Inkscape.app/Contents/MacOS/inkscape")
stopifnot(file.exists(compose), file.exists(ink))

W <- 180; FS <- 6.5; FAM <- "Helvetica"       # mm, pt
H <- c(a = 42, b = 52, c = 42)                  # panel heights, mm

# ---- data --------------------------------------------------------------------
s <- read.csv(file.path(root, "inst/benchmarks/results/performance-summary.csv"), stringsAsFactors = FALSE)
pk <- c("hilldiv3", "hilldiv2", "hillR", "entropart", "vegan", "BAT")
# Okabe-Ito; hilldiv3 = darkest/most saturated. None of the reserved #66C2A5/#FC8D62/#8DA0CB.
pal <- c(hilldiv3 = "#0072B2", hilldiv2 = "#E69F00", hillR = "#009E73",
         entropart = "#CC79A7", vegan = "#56B4E9", BAT = "#7F7F7F")
shp <- c(hilldiv3 = 16, hilldiv2 = 17, hillR = 15, entropart = 18, vegan = 4, BAT = 3)
ops <- c(B1 = "Alpha", B2 = "Flat partition", B3 = "Nested partition", B4 = "Pairwise matrix")
fcl <- c(neutral = "Neutral", phylogenetic = "Phylo.", functional = "Func.")
s <- s |> mutate(package = factor(package, pk),
                 op  = factor(ops[section], ops),
                 typ = factor(facet, names(fcl), fcl),
                 support = factor(support, c("yes", "partial", "no")))
stopifnot(nrow(s) == 72)

th <- theme_classic(base_size = FS, base_family = FAM) +
  theme(axis.line = element_line(linewidth = .25), axis.ticks = element_line(linewidth = .25),
        axis.text = element_text(size = FS, colour = "black"), axis.title = element_text(size = FS),
        strip.background = element_blank(), strip.text = element_text(size = FS, face = "bold"),
        legend.position = "none", plot.margin = margin(5, 2, 1, 2, "mm"),
        panel.spacing = unit(1.5, "mm"))

# ---- a: support matrix -------------------------------------------------------
cellfill <- c(yes = 1, partial = .35, no = 0)
a_df <- s |> mutate(pkg_rev = factor(package, rev(pk)), al = unname(cellfill[as.character(support)]),
                    col = ifelse(support == "no", "#EDEDED", pal[as.character(package)]))
pa <- ggplot(a_df, aes(typ, pkg_rev)) +
  geom_tile(aes(fill = I(col), alpha = I(ifelse(support == "no", 1, ifelse(support == "partial", .3, 1)))),
            colour = "white", linewidth = .6) +
  geom_tile(data = filter(a_df, support == "partial"), aes(colour = I(pal[as.character(package)])),
            fill = NA, linewidth = .35, width = .86, height = .86) +
  facet_grid(. ~ op, switch = "x") +
  labs(x = NULL, y = NULL) + th +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.y = element_text(face = "plain"),
        strip.placement = "outside",
        strip.text = element_text(face = "bold", margin = margin(3, 0, 3, 0)),
        axis.text.x = element_text(size = FS - .5, margin = margin(1.5, 0, 0, 0)))
leg <- data.frame(x = c(1, 2.9, 7.2), lab = c("Supported", "Partial (related method or different definition)", "Not supported"),
                  f = c("#4D4D4D", "#4D4D4D", "#EDEDED"), a = c(1, .3, 1))
pl <- ggplot(leg) +
  geom_tile(aes(x, 1, fill = I(f), alpha = I(a)), width = .35, height = .8) +
  geom_tile(data = leg[2, ], aes(x, 1), fill = NA, colour = "#4D4D4D", linewidth = .35, width = .3, height = .7) +
  geom_text(aes(x + .3, 1, label = lab), hjust = 0, size = FS / .pt, family = FAM) +
  scale_x_continuous(limits = c(.8, 9.5), expand = c(0, 0)) + theme_void() +
  theme(plot.margin = margin(3, 2, 0, 2, "mm"))
panel_a <- pa / pl + plot_layout(heights = c(1, .13))

# ---- shared offsets for b / c ------------------------------------------------
off <- setNames(seq(-2.5, 2.5) * .125, pk)
s <- s |> mutate(xn = as.integer(typ) + off[as.character(package)], pt = ifelse(package == "hilldiv3", 1.5, 1.1))
xs <- scale_x_continuous(breaks = 1:3, labels = unname(fcl), limits = c(.45, 3.55), expand = c(0, 0))
cols <- list(scale_colour_manual(values = pal, breaks = pk, drop = FALSE, name = NULL),
             scale_shape_manual(values = shp, breaks = pk, drop = FALSE, name = NULL))
band <- function(lo, hi) annotate("rect", xmin = 1.5, xmax = 2.5, ymin = lo, ymax = hi, fill = "grey93")

# ---- b: median time ----------------------------------------------------------
b_df <- filter(s, !is.na(median_seconds))
pb <- ggplot(b_df, aes(xn, median_seconds, colour = package, shape = package)) + band(1e-4, 1e3) +
  geom_linerange(aes(ymin = min_seconds, ymax = max_seconds), linewidth = .35) +
  geom_point(aes(size = I(pt))) +
  facet_wrap(~op, nrow = 1) + xs + cols +
  scale_y_log10(limits = c(1e-4, 1e3), expand = c(0, 0), breaks = 10^seq(-4, 2), labels = scales::label_log()) +
  annotation_logticks(sides = "l", linewidth = .2, short = unit(.6, "mm"), mid = unit(.9, "mm"), long = unit(1.2, "mm")) +
  labs(x = NULL, y = "Median time (s)") + th
# ---- c: memory (pairwise matrices) -------------------------------------------
c_df <- filter(s, section == "B4", !is.na(memory_bytes))
fmt <- function(b) ifelse(b >= 1e9, sprintf("%.1f GB", b / 1e9), sprintf("%.0f MB", b / 1e6))
lab_df <- filter(c_df, package %in% c("hilldiv3", "hilldiv2", "hillR")) |>
  mutate(lab_h = c(hilldiv3 = .5, hilldiv2 = .5, hillR = 0)[as.character(package)],   # fan out: left / centre / right
         lab_v = c(hilldiv3 = 1.9, hilldiv2 = -.9, hillR = 1.9)[as.character(package)]) # below / above / below
pc <- ggplot(c_df, aes(xn, memory_bytes, colour = package, shape = package)) + band(1e5, 1e12) +
  geom_point(data = mutate(s, memory_bytes = NA_real_), na.rm = TRUE) +
  geom_point(aes(size = I(pt))) +
  geom_text(data = lab_df, aes(label = fmt(memory_bytes), hjust = lab_h, vjust = lab_v), size = 5.5 / .pt,
            family = FAM, show.legend = FALSE) +
  xs + cols +
  scale_y_log10(limits = c(1e5, 1e12), expand = c(0, 0), breaks = c(1e6, 1e8, 1e10),
                labels = c("1 MB", "100 MB", "10 GB")) +
  annotation_logticks(sides = "l", linewidth = .2, short = unit(.6, "mm"), mid = unit(.9, "mm"), long = unit(1.2, "mm")) +
  labs(x = "Pairwise beta diversity matrix (q = 1)", y = "Allocated memory (bytes)") +
  th + theme(legend.position = "right", legend.key.size = unit(3, "mm"),
             legend.text = element_text(size = FS), legend.margin = margin(0, 0, 0, 0)) +
  guides(colour = guide_legend(override.aes = list(size = 1.4)))
panel_c <- pc + guide_area() + plot_layout(widths = c(1.45, 1), guides = "collect")

save <- function(p, name, h) {
  f <- file.path(fig, "raw", name)
  ggsave(f, p, width = W, height = h, units = "mm", device = svglite, fix_text_size = FALSE)
  f
}
fa <- save(panel_a, "fig2a.svg", H[["a"]]); fb <- save(pb, "fig2b.svg", H[["b"]]); fc <- save(panel_c, "fig2c.svg", H[["c"]])

# ---- compose + export --------------------------------------------------------
src <- file.path(fig, "src/fig2.svg")
system2("python3", c("-I", compose, src, "--width", W, "--cols", 1, "--gap", 2, "--label-size", 8, fa, fb, fc))
pdf <- file.path(fig, "final/fig2.pdf"); png <- file.path(fig, "final/fig2.png")
system2(ink, c(src, "--export-type=pdf", "--export-text-to-path", paste0("--export-filename=", pdf)), stderr = FALSE)
system2(ink, c(src, "--export-type=png", "--export-dpi=600", "--export-background=white", paste0("--export-filename=", png)), stderr = FALSE)
cat("Wrote", src, pdf, png, "\n")
