#' Tidy result objects for hilldiv3
#'
#' Internal helpers that turn the engine's wide matrices into the long-format
#' (tidy) data frames returned by the user-facing `hill*` functions, and the S3
#' `print()` / `plot()` / `autoplot()` methods for those results. Every result
#' carries a common parent class `hill_result` plus a specific subclass, so the
#' shared machinery (printing, line plots) is written once.
#'
#' @name hill_result
#' @keywords internal
NULL

# Turn a (q x key) matrix into a long data.frame with columns `q`, <key>,
# `value`. `as.vector()` walks the matrix column-major, so q varies fastest
# within each column -> rep(q, times) for q and rep(cols, each) for the key.
.hill_longify <- function(mat, q, key) {
  cols <- colnames(mat)
  df <- data.frame(
    q = rep(q, times = length(cols)),
    key = rep(cols, each = length(q)),
    value = as.vector(mat),
    stringsAsFactors = FALSE
  )
  names(df)[2] <- key
  df
}

# Attach the result classes and metadata used by the print/plot methods.
new_hill_result <- function(df, subclass, type = NULL, group = NULL,
                            value_label = "Hill number (qD)") {
  attr(df, "hill_type") <- type
  attr(df, "hill_group") <- group
  attr(df, "hill_value_label") <- value_label
  class(df) <- c(subclass, "hill_result", "data.frame")
  df
}

#' @export
#' @noRd
print.hill_result <- function(x, ...) {
  type <- attr(x, "hill_type")
  header <- "<hilldiv3 result>"
  if (!is.null(type)) header <- sprintf("<hilldiv3 result: %s>", type)
  cat(header, sprintf("%d rows x %d cols\n", nrow(x), ncol(x)), sep = "\n")
  df <- x
  attributes(df)[c("hill_type", "hill_group", "hill_value_label")] <- NULL
  class(df) <- "data.frame"
  print(df, ...)
  invisible(x)
}

# Shared base-graphics line plot: one line per level of the grouping column,
# value against the diversity order q.
.hill_lineplot <- function(x, group_col, ylab, ...) {
  groups <- unique(x[[group_col]])
  cols <- grDevices::hcl.colors(max(length(groups), 2), "Dark 3")
  plot(NA, xlim = range(x$q), ylim = range(x$value, na.rm = TRUE),
       xlab = "Diversity order (q)", ylab = ylab, ...)
  for (i in seq_along(groups)) {
    g <- x[x[[group_col]] == groups[i], ]
    g <- g[order(g$q), ]
    graphics::lines(g$q, g$value, col = cols[i], lwd = 2)
    graphics::points(g$q, g$value, col = cols[i], pch = 16)
  }
  graphics::legend("topright", legend = groups, col = cols[seq_along(groups)],
                   lwd = 2, bty = "n")
  invisible(x)
}

#' @exportS3Method graphics::plot hill_diversity
plot.hill_diversity <- function(x, ...) {
  # Combined neutral/phylo/functional output carries a `type` column; draw one
  # line per type-and-sample so the flavours do not collapse onto each other.
  if ("type" %in% names(x)) {
    x$.group <- paste(x$type, x$sample, sep = " / ")
    return(.hill_lineplot(x, ".group", attr(x, "hill_value_label"), ...))
  }
  .hill_lineplot(x, "sample", attr(x, "hill_value_label"), ...)
}

#' @exportS3Method graphics::plot hill_evenness
plot.hill_evenness <- function(x, ...) {
  .hill_lineplot(x, "sample", "Evenness", ...)
}

#' @exportS3Method graphics::plot hill_partition
plot.hill_partition <- function(x, ...) {
  .hill_lineplot(x, "component", "Diversity", ...)
}

#' @exportS3Method graphics::plot hill_dissimilarity
plot.hill_dissimilarity <- function(x, ...) {
  .hill_lineplot(x, "metric", "Dissimilarity", ...)
}

#' @exportS3Method graphics::plot hill_hierarchy
plot.hill_hierarchy <- function(x, ...) {
  # One beta line per nesting level; the finest (sample) scale has no beta.
  b <- x[!is.na(x$beta), ]
  b$value <- b$beta
  b$scale <- factor(b$scale)
  .hill_lineplot(b, "scale", "Beta (turnover)", ...)
}

#' @exportS3Method graphics::plot hill_similarity
plot.hill_similarity <- function(x, ...) {
  .hill_lineplot(x, "metric", "Similarity", ...)
}

# ggplot2 autoplot worker: a single grouped profile geom. Registered for each
# result subclass in .onLoad() so ggplot2 stays a Suggests-only dependency.
.hill_autoplot <- function(object, group_col, ylab, ...) {
  rlang::check_installed("ggplot2", "for `autoplot()` methods.")
  ggplot2::ggplot(
    object,
    ggplot2::aes(x = .data$q, y = .data$value,
                 colour = .data[[group_col]], group = .data[[group_col]])
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point() +
    ggplot2::labs(x = "Diversity order (q)", y = ylab, colour = NULL)
}

autoplot.hill_diversity <- function(object, ...) {
  # Facet the combined output by diversity type; colour stays per sample.
  if ("type" %in% names(object)) {
    rlang::check_installed("ggplot2", "for `autoplot()` methods.")
    return(
      ggplot2::ggplot(
        object,
        ggplot2::aes(x = .data$q, y = .data$value,
                     colour = .data$sample, group = .data$sample)
      ) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::geom_point() +
        ggplot2::facet_wrap(ggplot2::vars(.data$type), scales = "free_y") +
        ggplot2::labs(x = "Diversity order (q)",
                      y = attr(object, "hill_value_label"), colour = NULL)
    )
  }
  .hill_autoplot(object, "sample", attr(object, "hill_value_label"), ...)
}
autoplot.hill_profile <- function(object, ...) {
  .hill_autoplot(object, "sample", "Hill number (qD)", ...)
}
autoplot.hill_evenness <- function(object, ...) {
  .hill_autoplot(object, "sample", "Evenness", ...)
}
autoplot.hill_partition <- function(object, ...) {
  .hill_autoplot(object, "component", "Diversity", ...)
}
autoplot.hill_hierarchy <- function(object, ...) {
  rlang::check_installed("ggplot2", "for `autoplot()` methods.")
  b <- object[!is.na(object$beta), ]
  b$scale <- factor(b$scale)
  ggplot2::ggplot(
    b, ggplot2::aes(x = .data$q, y = .data$beta,
                    colour = .data$scale, group = .data$scale)
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point() +
    ggplot2::labs(x = "Diversity order (q)", y = "Beta (turnover)",
                  colour = NULL)
}
autoplot.hill_dissimilarity <- function(object, ...) {
  .hill_autoplot(object, "metric", "Dissimilarity", ...)
}
autoplot.hill_similarity <- function(object, ...) {
  .hill_autoplot(object, "metric", "Similarity", ...)
}
autoplot.hill_redundancy <- function(object, ...) {
  rlang::check_installed("ggplot2", "for `autoplot()` methods.")
  fit <- attr(object, "hill_fit")
  if (is.null(fit)) {
    cli::cli_abort("This {.cls hill_redundancy} object carries no fit data.")
  }
  # One fitted saturating curve per order, evaluated over each order's x-range.
  curve <- do.call(rbind, lapply(split(object, object$q), function(r) {
    pts <- fit[fit$q == r$q, ]
    if (nrow(pts) == 0 || !all(is.finite(c(r$a, r$b, r$c)))) return(NULL)
    xx <- seq(min(pts$neutral), max(pts$neutral), length.out = 100)
    data.frame(q = r$q, neutral = xx, diversity = -r$a * 2^(-xx / r$b) + r$c)
  }))
  fit$q <- factor(fit$q)
  curve$q <- factor(curve$q)
  ggplot2::ggplot(
    mapping = ggplot2::aes(x = .data$neutral, y = .data$diversity,
                           colour = .data$q)
  ) +
    ggplot2::geom_point(data = fit) +
    ggplot2::geom_line(data = curve, linewidth = 1) +
    ggplot2::labs(x = "Neutral diversity (qD)",
                  y = .hill_red_ylab(attr(object, "hill_type")),
                  colour = "Order (q)")
}
