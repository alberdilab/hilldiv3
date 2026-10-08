#' Allocate Hill power-sum weight to taxa or sets
#'
#' Calculate the share of a community's neutral, phylogenetic, or functional
#' Hill power sum assigned to each taxon or a specified set of taxa. Shares
#' range from zero to one within a nonempty sample. They describe an additive
#' allocation of the power sum, not a percentage of the effective Hill number.
#'
#' @param data Counts: a numeric vector or a taxa-by-samples matrix/data frame.
#'   The input adapters accepted by [hilldiv()] are also supported.
#' @param q Nonnegative diversity orders. Defaults to c(0, 1, 2).
#' @param tree Optional rooted phylogenetic tree of class phylo.
#' @param dist Optional functional distance matrix over the taxa.
#' @param tau Functional distance threshold. Defaults to max(dist).
#' @param type Diversity type(s). "auto" (default) computes neutral shares and
#'   also phylogenetic and functional shares when their inputs are supplied.
#'   Explicit types restrict the output as in [hilldiv()].
#' @param sets A named list whose elements are character vectors of taxon names.
#'   Each element can contain one taxon or several. By default, every taxon is
#'   returned as a singleton set. Sets may overlap, and each set's share is the
#'   sum of its members' allocations. An unnamed character vector is accepted
#'   as a request for singleton taxa.
#' @param allocation Allocation rule. The first iteration supports "qpower".
#'
#' @details
#' Let p_i be the relative abundance of taxon i. Neutral weight is p_i^q,
#' normalised by its sum over present taxa.
#'
#' For a phylogenetic branch b with length L_b and descendant abundance a_b,
#' branch weight is L_b a_b^q. That weight is divided among its present
#' descendants in proportion to p_i^q. At q = 0, present descendants share
#' the branch equally. Branch weights and taxon allocations are normalised by
#' the sum of all branch weights.
#'
#' For functional diversity, similarity is s_ij = 1 - min(d_ij, tau) / tau and
#' a_i = sum_j s_ij p_j. The existing functional Hill power sum in hilldiv3 is
#' sum_i p_i a_i^(q-1), so taxon i receives
#' p_i a_i^(q-1) / sum_j p_j a_j^(q-1). With identity similarity this reduces
#' to the neutral rule. At q = 1, functional shares equal p_i, although the
#' effective functional Hill number still depends on similarities.
#'
#' An empty sample has no power sum to allocate, so its shares are NA.
#'
#' @return A data frame with columns q, sample, type, set, share, and
#'   allocation. Shares for distinct singleton taxa sum to one in each
#'   nonempty sample, diversity type, and order.
#'
#' @references
#' Chao, A., Chiu, C.-H. & Jost, L. (2010). Phylogenetic diversity measures
#' based on Hill numbers. Phil. Trans. R. Soc. B, 365, 3599-3609.
#'
#' Leinster, T. & Cobbold, C.A. (2012). Measuring diversity: the importance
#' of species similarity. Ecology, 93, 477-489.
#'
#' @seealso [hilldiv()], [hillprof()]
#' @examples
#' counts <- c(A = 394, B = 202, C = 354, D = 50)
#' tree <- ape::read.tree(
#'   text = "(((A:0.300,B:0.300):0.185,C:0.485):0.009,D:0.494);"
#' )
#' hillshare(counts, q = 3, tree = tree,
#'           sets = list(C = "C", ABC = c("A", "B", "C")))
#' d <- as.matrix(stats::dist(c(A = 0, B = 0.5, C = 1, D = 2)))
#' hillshare(counts, q = c(0, 1, 2), dist = d, type = "functional",
#'           sets = list(C = "C", ABC = c("A", "B", "C")))
#' @export
hillshare <- function(data, q = c(0, 1, 2), tree = NULL, dist = NULL,
                      tau = NULL,
                      type = c("auto", "neutral", "phylogenetic", "functional"),
                      sets = NULL, allocation = "qpower") {
  if (!is.numeric(q) || length(q) == 0L || anyNA(q) ||
      any(!is.finite(q)) || any(q < 0)) {
    cli::cli_abort("{.arg q} must contain finite, nonnegative numbers.")
  }
  if (!identical(allocation, "qpower")) {
    cli::cli_abort("The available {.arg allocation} is {.val qpower}.")
  }
  type <- match.arg(type, several.ok = TRUE)

  xin <- as_hill_input(data, tree = tree, dist = dist)
  if (!is.numeric(xin$counts)) {
    cli::cli_abort("{.arg data} must contain numeric counts.")
  }
  if (any(!is.finite(xin$counts)) || any(xin$counts < 0)) {
    cli::cli_abort("{.arg data} must contain finite, nonnegative counts.")
  }
  if (is.null(rownames(xin$counts))) {
    rownames(xin$counts) <- paste0("taxon", seq_len(nrow(xin$counts)))
  }
  if (anyDuplicated(rownames(xin$counts))) {
    cli::cli_abort("Taxon names in {.arg data} must be unique.")
  }
  if (is.null(colnames(xin$counts))) {
    colnames(xin$counts) <- paste0("sample", seq_len(ncol(xin$counts)))
  }
  taxa <- rownames(xin$counts)
  sets <- .share_sets(sets, taxa)
  types <- .hilldiv_types(xin, type)

  out <- vector("list", length(types) * ncol(xin$counts) * length(q))
  row <- 0L
  for (ty in types) {
    xt <- xin
    if (ty != "phylogenetic") xt$tree <- NULL
    if (ty != "functional") xt$dist <- NULL
    x <- prep_data(xt, q, ty)
    p <- tss(x$counts)
    scores <- switch(ty,
      neutral = .share_neutral(p, q),
      phylogenetic = .share_phylogenetic(p, q, x$tree),
      functional = .share_functional(p, q, x$dist, tau)
    )
    for (j in seq_len(ncol(p))) {
      for (k in seq_along(q)) {
        row <- row + 1L
        shares <- vapply(sets, function(members) {
          sum(scores[match(members, rownames(p)), j, k])
        }, numeric(1))
        out[[row]] <- data.frame(
          q = rep(q[k], length(sets)),
          sample = rep(colnames(p)[j], length(sets)),
          type = rep(ty, length(sets)),
          set = names(sets),
          share = unname(shares),
          allocation = rep("qpower", length(sets)),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  ans <- do.call(rbind, out)
  rownames(ans) <- NULL
  ans
}

.share_sets <- function(sets, taxa) {
  if (is.null(sets)) {
    sets <- as.list(taxa)
    names(sets) <- taxa
    return(sets)
  }
  if (is.character(sets) && is.null(names(sets))) {
    sets <- as.list(sets)
    names(sets) <- unlist(sets, use.names = FALSE)
  }
  if (!is.list(sets) || length(sets) == 0L ||
      is.null(names(sets)) || anyNA(names(sets)) ||
      any(names(sets) == "") || anyDuplicated(names(sets))) {
    cli::cli_abort(
      "{.arg sets} must be a nonempty named list with unique names."
    )
  }
  for (nm in names(sets)) {
    members <- sets[[nm]]
    if (!is.character(members) || length(members) == 0L ||
        anyNA(members) || anyDuplicated(members)) {
      cli::cli_abort("Set {.val {nm}} needs unique, nonempty taxon names.")
    }
    missing <- setdiff(members, taxa)
    if (length(missing) > 0L) {
      cli::cli_abort("Unknown taxa in set {.val {nm}}: {.val {missing}}.")
    }
  }
  sets
}

.share_array <- function(p, q) {
  array(NA_real_, dim = c(nrow(p), ncol(p), length(q)),
        dimnames = list(rownames(p), colnames(p), paste0("q", q)))
}

.share_neutral <- function(p, q) {
  out <- .share_array(p, q)
  for (j in seq_len(ncol(p))) {
    for (k in seq_along(q)) {
      w <- if (q[k] == 0) as.numeric(p[, j] > 0) else p[, j]^q[k]
      z <- sum(w)
      if (z > 0) out[, j, k] <- w / z
    }
  }
  out
}

.share_phylogenetic <- function(p, q, tree) {
  if (is.null(tree$edge.length) ||
      any(!is.finite(tree$edge.length)) ||
      any(tree$edge.length < 0)) {
    cli::cli_abort("The {.arg tree} needs finite, nonnegative branch lengths.")
  }
  tree <- ape::reorder.phylo(tree, order = "cladewise")
  p <- p[tree$tip.label, , drop = FALSE]
  ba <- branch_abundance(tree, p)
  Li <- ba$Li
  out <- .share_array(p, q)
  n_tip <- nrow(p)
  parent <- tree$edge[, 1]
  child <- tree$edge[, 2]

  for (k in seq_along(q)) {
    pq <- if (q[k] == 0) (p > 0) * 1 else p^q[k]
    d_all <- branch_abundance(tree, pq)$ai
    for (j in seq_len(ncol(p))) {
      if (sum(p[, j]) == 0) next
      a <- ba$ai[, j]
      d <- d_all[, j]
      w <- numeric(length(Li))
      present <- a > 0
      w[present] <- Li[present] * a[present]^q[k]
      z <- sum(w)
      if (z == 0) next
      coef <- numeric(length(Li))
      keep <- d > 0
      coef[keep] <- w[keep] / d[keep]
      path <- numeric(n_tip + tree$Nnode)
      for (e in seq_along(parent)) {
        path[child[e]] <- path[parent[e]] + coef[e]
      }
      out[, j, k] <- pq[, j] * path[seq_len(n_tip)] / z
    }
  }
  out
}

.share_functional <- function(p, q, dist, tau) {
  sim <- .functional_similarity(dist, tau)
  out <- .share_array(p, q)
  for (j in seq_len(ncol(p))) {
    pj <- p[, j]
    present <- pj > 0
    if (!any(present)) next
    a <- as.vector(sim %*% pj)
    if (any(a[present] <= 0)) {
      cli::cli_abort("Functional similarity is zero for a present taxon.")
    }
    for (k in seq_along(q)) {
      w <- numeric(length(pj))
      w[present] <- pj[present] * a[present]^(q[k] - 1)
      z <- sum(w)
      if (z > 0) out[, j, k] <- w / z
    }
  }
  out
}
