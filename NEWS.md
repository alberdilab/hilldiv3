# hilldiv3 (development version)

* **Breaking:** the `out = "matrix"` form of the per-sample functions
  (`hilldiv()`, `hillprof()`, `hilleven()`) now returns samples in rows and
  diversity orders (`q0`, `q1`, ...) in columns, instead of the previous
  orders-in-rows/samples-in-columns layout. This matches the usual
  observations-in-rows convention and joins cleanly with per-sample metadata.
  Code that indexed these matrices by `["q0", ]` should switch to `[, "q0"]`
  (or transpose). The component/metric matrices from `hillpart()`,
  `hilldiss()` and `hillred()` are unchanged.
* `hilldiv()` and friends now accept count tables whose first column holds
  taxa names: a leading non-numeric column is automatically promoted to row
  names (instead of failing with "Count data must be numeric"). Any other
  non-numeric column raises a clear error asking the user to fix the input.
* Verbose output now reports the number of taxa and samples being analysed.

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
* **`reference` argument on `hilldiv()`** — `reference = c("pool", "sample")`
  selects the reference tree depth for *phylogenetic* Hill numbers. `"pool"`
  (default) reads every sample at one common depth `T = mean(T_j)` so values
  are mutually comparable (hilldiv2's `multi` behaviour); `"sample"` reads each
  sample at its own depth `T_j`. The two coincide on ultrametric trees. The
  option is deliberately absent from `hillpart()`, where `T` is fixed at the
  Chiu et al. (2014) mean per-sample depth — the unique value for which
  `gamma / alpha` is a valid decomposition with `beta` in `[1, N]`.

### Bug fixes
* **Phylogenetic `hilldiv()` alpha corrected.** The per-sample phylogenetic
  Hill number previously raised branch length `L_i` to the power `q` along with
  abundance (`(L_i a_i / T)^q`), which matched no standard quantity and
  disagreed with the partition engine. Branch length is now a linear weight
  (`(L_i / T) a_i^q`, Chao et al. 2010), so `q = 0` recovers Faith's PD / T and
  the per-sample value matches `hillpart()` alpha at a single sample. **This
  changes phylogenetic `hilldiv()` and `hillprof()` outputs** relative to
  earlier 3.0.0 development snapshots.

### Infrastructure
* Added test-coverage (Codecov) and lint (lintr) GitHub Actions workflows.
* Added golden-value tests cross-checking the engine against `vegan` and
  hand-computed constants, plus edge cases (single taxon, empty sample, q = 1).

### Compatibility
* The `hilldiv()`, `hillpart()`, `hilldiss()`, `hillsim()`, `hillpair()`,
  `hillred()`, `tss()` and `traits2dist()` names from hilldiv2 are preserved.
  The default return shape is now a tidy `data.frame`; use `out = "matrix"`
  for the hilldiv2-style matrix.
