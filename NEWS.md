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

### Compatibility
* The `hilldiv()`, `hillpart()`, `hilldiss()`, `hillsim()`, `hillpair()`,
  `hillred()`, `tss()` and `traits2dist()` names and core behaviour from
  hilldiv2 are preserved.
