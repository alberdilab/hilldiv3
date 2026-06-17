# hilldiv3 (development version)

## hilldiv3 3.0.0

Complete redesign of the package built on a tested, isolated compute engine.

### Architecture
* Math separated into an internal, unit-tested engine (`hill_alpha()`,
  `hill_partition()`, `hill_beta_to_dissim()`); the user-facing `hill*`
  functions are thin wrappers around it.
* Single validation/alignment layer (`prep_data()`) that checks names with
  `setequal()` and **reorders** data to match the tree/distance matrix,
  fixing the silent-misalignment behaviour of hilldiv2.
* New `as_hill_input()` adapter accepting matrices, data frames, tibbles,
  `phyloseq` and `TreeSummarizedExperiment` objects.
* Fast `ape`-based post-order tree traversal replaces `geiger::tips()`.

### New features
* `hillprof()` — diversity profiles across a sweep of q values.
* `hilleven()` — evenness from Hill numbers.
* Real `match_data()` helper (previously only referenced in the docs).
* **Tidy output by default.** `hilldiv()`, `hillpart()`, `hilldiss()`,
  `hillsim()`, `hilleven()`, `hillprof()` and `hillred()` now return
  long-format `data.frame`s with `print()`, `plot()` and (when `ggplot2` is
  installed) `autoplot()` methods. Pass `out = "matrix"` for the legacy shape.
* **Explicit `type` argument** — `type = c("auto", "neutral", "phylogenetic",
  "functional")` on every entry point; `"auto"` keeps input-based detection,
  an explicit value asserts and validates the diversity type.
* `hillpair()` computes the type-specific structure once over all samples and
  reuses it per pair (no more full re-partition per pair), and reports a
  `progressr` bar when available.
* Bundled, documented example data: `gut_counts`, `gut_tree`, `gut_traits`.

### Infrastructure
* Added test-coverage (Codecov) and lint (lintr) GitHub Actions workflows.
* Added golden-value tests cross-checking the engine against `vegan` and
  hand-computed constants, plus edge cases (single taxon, empty sample, q = 1).

### Compatibility
* The `hilldiv()`, `hillpart()`, `hilldiss()`, `hillsim()`, `hillpair()`,
  `hillred()`, `tss()` and `traits2dist()` names from hilldiv2 are preserved.
  The default return shape is now a tidy `data.frame`; use `out = "matrix"`
  for the hilldiv2-style matrix.
