#' Attribute Hill dissimilarity to taxa or sets
#'
#' Show which taxa contribute most to dissimilarity between every pair of
#' samples or among samples in a group. Default singleton-taxon contributions
#' are nonnegative and add up to the corresponding value from [hillpair()] or
#' [hilldiss()]. A taxon's
#' `share` is its fraction of that dissimilarity; it is `NA` when the
#' dissimilarity is zero.
#'
#' @inheritParams hilldiss
#' @param by `"pairwise"` (default) compares every pair of samples;
#'   `"collective"` compares all samples together or each specified group.
#' @param groups For `by = "collective"`, an optional named list of character
#'   vectors of sample names. Each group must contain at least two distinct
#'   samples. Defaults to one group named `"all"` containing every sample.
#' @param sets An optional named list of character vectors of taxon names,
#'   as in [hillshare()]. Defaults to one set per taxon. Sets may overlap.
#'
#' @details
#' For each taxon (or phylogenetic branch), the contribution kernel is the
#' Jensen gap between its pooled abundance and its abundances across samples.
#' The gap is zero if its abundance is constant. At `q = 0` the gap measures
#' uneven occurrence; at `q = 1` it is the entropy gap; at other orders it is
#' the power-sum gap. For neutral and functional diversity the kernel belongs
#' directly to a taxon. For phylogenetic diversity each branch's gap is
#' assigned to descendant taxa in proportion to their own neutral Jensen
#' gaps at the same order. This is one explicit allocation of shared branches,
#' rather than a unique taxon-level decomposition of phylogenetic turnover.
#'
#' The kernel fractions are multiplied by the selected dissimilarity metric.
#' Thus singleton-taxon contributions add exactly to that metric, while their
#' relative ranking is the same for `S`, `C`, `U`, and `V` at a given order and
#' comparison.
#' Zero-length branches and taxa absent from every selected sample receive
#' zero. Every selected sample must have a positive total. Functional beta
#' uses the package's raw-count partitioning, so unequal sample totals can
#' contribute to functional dissimilarity even when relative composition is
#' unchanged.
#'
#' @return A data frame with columns `first`, `second`, `group`, `q`, `type`,
#'   `metric`, `set`, `contribution`, `share`, and `rank`. Pairwise rows use
#'   `first` and `second`; collective rows use `group`. Rows are ranked from
#'   largest to smallest contribution within each comparison, order, and
#'   metric. Sets that omit taxa sum to less than the total dissimilarity;
#'   overlapping sets can sum to more.
#' @seealso [hillshare()], [hillpair()], [hilldiss()]
#' @examples
#' counts <- matrix(c(5, 0, 2, 3, 4, 0, 0, 2, 6), nrow = 3,
#'                  dimnames = list(c("A", "B", "C"), c("s1", "s2", "s3")))
#' hillcontrib(counts, q = 1, metric = "C")
#' hillcontrib(counts, q = 1, metric = "C", by = "collective")
#' @export
hillcontrib <- function(data, q = c(0, 1, 2), metric = c("S", "C", "U", "V"),
                        tree = NULL, dist = NULL, tau = NULL,
                        type = c("auto", "neutral", "phylogenetic",
                                 "functional"),
                        by = c("pairwise", "collective"), groups = NULL,
                        sets = NULL) {
  if (!is.numeric(q) || length(q) == 0L || anyNA(q) ||
      any(!is.finite(q)) || any(q < 0)) {
    cli::cli_abort("{.arg q} must contain finite, nonnegative numbers.")
  }
  metric <- match.arg(metric, c("S", "C", "U", "V"), several.ok = TRUE)
  type <- match.arg(type)
  by <- match.arg(by)
  x <- as_hill_input(data, tree = tree, dist = dist)
  if (!is.numeric(x$counts) || any(!is.finite(x$counts)) ||
      any(x$counts < 0)) {
    cli::cli_abort("{.arg data} must contain finite, nonnegative counts.")
  }
  if (is.null(rownames(x$counts))) {
    rownames(x$counts) <- paste0("taxon", seq_len(nrow(x$counts)))
  }
  if (anyDuplicated(rownames(x$counts))) {
    cli::cli_abort("Taxon names in {.arg data} must be unique.")
  }
  if (is.null(colnames(x$counts))) {
    colnames(x$counts) <- paste0("sample", seq_len(ncol(x$counts)))
  }
  if (anyDuplicated(colnames(x$counts))) {
    cli::cli_abort("Sample names in {.arg data} must be unique.")
  }
  x <- prep_data(x, q, type)
  type <- attr(x, "type")
  if (ncol(x$counts) < 2L) {
    cli::cli_abort("Dissimilarity attribution needs at least two samples.")
  }
  if (type == "phylogenetic" &&
      (is.null(x$tree$edge.length) ||
       any(!is.finite(x$tree$edge.length)) ||
       any(x$tree$edge.length < 0))) {
    cli::cli_abort("The {.arg tree} needs finite, nonnegative branch lengths.")
  }
  sets <- .share_sets(sets, rownames(x$counts))
  comparisons <- .contrib_comparisons(by, groups, colnames(x$counts))
  prep <- part_prep(x$counts, type, tree = x$tree, dist = x$dist, tau = tau)
  descendants <- if (type == "phylogenetic") {
    .contrib_descendants(x$tree)
  } else {
    NULL
  }

  out <- vector("list", length(comparisons) * length(q) * length(metric))
  row <- 0L
  for (nm in names(comparisons)) {
    cols <- comparisons[[nm]]
    if (any(colSums(x$counts[, cols, drop = FALSE]) == 0)) {
      cli::cli_abort("Every sample in comparison {.val {nm}} must have a
                      positive total.")
    }
    if (type == "phylogenetic" &&
        sum(prep$Li * rowSums(prep$aij[, cols, drop = FALSE])) == 0) {
      cli::cli_abort("Comparison {.val {nm}} has no positive phylogenetic
                      branch length.")
    }
    fractions <- .contrib_fractions(prep, cols, q, x$counts,
                                    descendants, rownames(x$counts))
    betas <- part_eval(prep, cols, q)[, "beta"]
    for (k in seq_along(q)) {
      values <- beta_to_dissim(unname(betas[k]), length(cols), q[k])
      for (met in metric) {
        row <- row + 1L
        contributions <- vapply(sets, function(members) {
          sum(fractions[match(members, rownames(x$counts)), k]) * values[[met]]
        }, numeric(1))
        shares <- if (values[[met]] > 0) contributions / values[[met]] else
          rep(NA_real_, length(sets))
        ord <- order(-contributions, seq_along(contributions))
        result <- data.frame(
          first = if (by == "pairwise")
            colnames(x$counts)[cols[1]] else NA_character_,
          second = if (by == "pairwise")
            colnames(x$counts)[cols[2]] else NA_character_,
          group = if (by == "collective") nm else NA_character_,
          q = rep(q[k], length(sets)), type = rep(type, length(sets)),
          metric = rep(met, length(sets)), set = names(sets),
          contribution = unname(contributions), share = unname(shares),
          stringsAsFactors = FALSE
        )[ord, , drop = FALSE]
        result$rank <- if (values[[met]] > 0) seq_along(sets) else
          rep(NA_integer_, length(sets))
        out[[row]] <- result
      }
    }
  }
  ans <- do.call(rbind, out)
  rownames(ans) <- NULL
  ans
}

.contrib_comparisons <- function(by, groups, samples) {
  if (by == "pairwise") {
    if (!is.null(groups)) {
      cli::cli_abort("{.arg groups} is available only with
                      {.code by = \"collective\"}.")
    }
    pairs <- utils::combn(samples, 2, simplify = FALSE)
    names(pairs) <- vapply(pairs, paste, collapse = " / ", FUN.VALUE = "")
    return(lapply(pairs, match, table = samples))
  }
  if (is.null(groups)) groups <- list(all = samples)
  if (!is.list(groups) || length(groups) == 0L || is.null(names(groups)) ||
      anyNA(names(groups)) || any(names(groups) == "") ||
      anyDuplicated(names(groups))) {
    cli::cli_abort("{.arg groups} must be a nonempty named list with unique
                    names.")
  }
  for (nm in names(groups)) {
    members <- groups[[nm]]
    if (!is.character(members) || length(members) < 2L || anyNA(members) ||
        anyDuplicated(members) || any(!members %in% samples)) {
      cli::cli_abort("Group {.val {nm}} needs at least two distinct known
                      sample names.")
    }
  }
  lapply(groups, match, table = samples)
}

# Per-category nonnegative Jensen gaps. Columns are samples; weights are
# branch lengths, attribute weights, or ones for neutral taxa.
.contrib_gap <- function(values, weights, q) {
  pooled <- rowMeans(values)
  if (q == 1) {
    xlogx <- function(x) {
      out <- numeric(length(x))
      keep <- x > 0
      out[keep] <- x[keep] * log(x[keep])
      out
    }
    gap <- rowMeans(matrix(xlogx(values), nrow(values))) - xlogx(pooled)
  } else {
    powered <- if (q == 0) (values > 0) * 1 else values^q
    pooled_power <- if (q == 0) (pooled > 0) * 1 else pooled^q
    gap <- rowMeans(powered) - pooled_power
    if (q < 1) gap <- -gap
  }
  pmax(0, weights * gap)
}

.contrib_fractions <- function(prep, cols, q, counts, descendants, taxa) {
  N <- length(cols)
  tax_values <- if (prep$type == "phylogenetic")
    tss(counts[, cols, drop = FALSE]) else NULL
  if (prep$type == "neutral") {
    values <- prep$pi[, cols, drop = FALSE]
    weights <- rep(1, nrow(values))
  } else if (prep$type == "phylogenetic") {
    aij <- prep$aij[, cols, drop = FALSE]
    Tval <- sum(prep$Li * rowSums(aij))
    values <- N * aij / Tval
    weights <- prep$Li
  } else {
    aik <- prep$aik[, cols, drop = FALSE]
    p <- prep$p[, cols, drop = FALSE]
    nplus <- sum(p)
    values <- N * aik / nplus
    aiplus <- rowSums(aik)
    weights <- numeric(length(aiplus))
    keep <- aiplus > 0
    weights[keep] <- rowSums(p)[keep] / aiplus[keep]
  }
  out <- matrix(0, length(taxa), length(q), dimnames = list(taxa, NULL))
  for (k in seq_along(q)) {
    gap <- .contrib_gap(values, weights, q[k])
    if (prep$type == "phylogenetic") {
      tax_gap <- .contrib_gap(tax_values, rep(1, length(taxa)), q[k])
      allocated <- numeric(length(taxa))
      for (edge in seq_along(gap)) {
        if (gap[edge] == 0) next
        members <- descendants[[edge]]
        w <- tax_gap[members]
        if (sum(w) == 0) w <- rep(1, length(members))
        allocated[members] <- allocated[members] + gap[edge] * w / sum(w)
      }
      gap <- allocated
    }
    if (sum(gap) > 0) out[, k] <- gap / sum(gap)
  }
  out
}

.contrib_descendants <- function(tree) {
  n_tip <- length(tree$tip.label)
  nodes <- vector("list", n_tip + tree$Nnode)
  for (tip in seq_len(n_tip)) nodes[[tip]] <- tip
  po <- ape::reorder.phylo(tree, "postorder")
  for (edge in seq_len(nrow(po$edge))) {
    parent <- po$edge[edge, 1]
    child <- po$edge[edge, 2]
    nodes[[parent]] <- c(nodes[[parent]], nodes[[child]])
  }
  lapply(tree$edge[, 2], function(node) nodes[[node]])
}
