q_set <- c(0, 0.5, 1, 2)

# Fixed count table and a symmetric functional distance matrix (zero diagonal,
# values in (0, 1)) so expectations are deterministic.
counts <- matrix(
  c(8, 0, 4, 2,
    0, 5, 3, 1,
    2, 2, 0, 6),
  nrow = 4,
  dimnames = list(c("t1", "t2", "t3", "t4"), c("s1", "s2", "s3"))
)
dist <- matrix(
  c(0.0, 0.4, 0.9, 0.7,
    0.4, 0.0, 0.6, 0.5,
    0.9, 0.6, 0.0, 0.3,
    0.7, 0.5, 0.3, 0.0),
  nrow = 4,
  dimnames = list(rownames(counts), rownames(counts))
)

# Verbatim port of hilldiv2's hillpart.functional (R/hillpart.R), used as the
# reference implementation for cross-checking the engine.
hillpart_functional_ref <- function(data, q = c(0, 1, 2), dist, tau) {
  if (missing(tau)) tau <- max(dist)
  N <- ncol(data)
  dij <- as.matrix(dist)
  dij[which(dij > tau, arr.ind = TRUE)] <- tau
  aik <- apply(data, 2, function(col) as.vector((1 - dij / tau) %*% col))
  aiplus <- apply(aik, 1, sum)
  vi <- apply(data, 1, sum) / aiplus
  alpha_v <- rep(vi, N)
  nplus <- sum(data)
  aik <- as.vector(aik)
  alpha_v <- alpha_v[aik != 0]
  aik <- aik[aik != 0]
  results <- matrix(0, nrow = length(q), ncol = 3)
  for (r in seq_along(q)) {
    qvalue <- q[r]
    if (qvalue == 1) {
      alpha <- 1 / N * exp(sum(-alpha_v * aik / nplus * log(aik / nplus)))
      gamma <- exp(sum(-vi * aiplus / nplus * log(aiplus / nplus)))
      beta <- gamma / alpha
    } else {
      alpha <- 1 / N * (sum(alpha_v * (aik / nplus)^qvalue))^(1 / (1 - qvalue))
      gamma <- (sum(vi * (aiplus / nplus)^qvalue))^(1 / (1 - qvalue))
      beta <- gamma / alpha
    }
    results[r, 1] <- alpha
    results[r, 2] <- gamma
    results[r, 3] <- beta
  }
  rownames(results) <- paste0("q", q)
  colnames(results) <- c("alpha", "gamma", "beta")
  results
}

test_that("functional partition matches hilldiv2 hillpart.functional", {
  ref <- hillpart_functional_ref(counts, q = q_set, dist = dist)
  got <- suppressMessages(hillpart(counts, q = q_set, dist = dist,
                                   out = "matrix"))
  expect_equal(got, ref, tolerance = 1e-12)
})

test_that("functional partition matches hilldiv2 with a custom tau", {
  tau <- 0.5
  ref <- hillpart_functional_ref(counts, q = q_set, dist = dist, tau = tau)
  got <- suppressMessages(hillpart(counts, q = q_set, dist = dist, tau = tau,
                                   out = "matrix"))
  expect_equal(got, ref, tolerance = 1e-12)
})

test_that("functional partition has the expected structure", {
  part <- suppressMessages(hillpart(counts, q = q_set, dist = dist,
                                    out = "matrix"))
  N <- ncol(counts)

  # beta lies in [1, N]; gamma >= alpha.
  expect_true(all(part[, "beta"] >= 1 - 1e-9))
  expect_true(all(part[, "beta"] <= N + 1e-9))
  expect_true(all(part[, "gamma"] >= part[, "alpha"] - 1e-9))
})

test_that("identical samples give no functional turnover (beta = 1)", {
  one <- counts[, 1, drop = FALSE]
  twin <- cbind(s1 = one[, 1], s2 = one[, 1])
  rownames(twin) <- rownames(counts)
  part <- suppressMessages(hillpart(twin, q = q_set, dist = dist,
                                    out = "matrix"))
  expect_equal(unname(part[, "beta"]), rep(1, length(q_set)))
  expect_equal(unname(part[, "alpha"]), unname(part[, "gamma"]))
})

test_that("q = 1 functional partition is the limit of nearby q", {
  lim <- suppressMessages(hillpart(counts, q = 1, dist = dist, out = "matrix"))
  near <- suppressMessages(hillpart(counts, q = 1 + 1e-6, dist = dist,
                                    out = "matrix"))
  expect_equal(unname(lim), unname(near), tolerance = 1e-4)
})

test_that("hilldiss/hillsim wire through the functional path", {
  d <- suppressMessages(hilldiss(counts, q = q_set, dist = dist,
                                 out = "matrix"))
  s <- suppressMessages(hillsim(counts, q = q_set, dist = dist,
                                out = "matrix"))
  expect_equal(d + s, matrix(1, nrow(d), ncol(d), dimnames = dimnames(d)))
  expect_true(all(d >= -1e-9 & d <= 1 + 1e-9))

  # beta = 1 (identical samples) implies zero dissimilarity on every metric.
  twin <- cbind(counts[, 1], counts[, 1])
  rownames(twin) <- rownames(counts)
  d0 <- suppressMessages(hilldiss(twin, q = q_set, dist = dist,
                                  out = "matrix"))
  expect_equal(unname(d0), matrix(0, length(q_set), 4))
})

test_that("functional partitioning requires a distance matrix", {
  expect_error(hilldiv3:::part_prep(counts, "functional", dist = NULL))
})
