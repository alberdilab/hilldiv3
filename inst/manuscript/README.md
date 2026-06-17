# hilldiv3 manuscript use cases

Two worked use cases for the manuscript introducing `hilldiv3`, each designed
to exercise a broad slice of the package on a realistic study design. The prose
(`.md`) is written to drop into the paper; the figures and every cited number
are produced by the companion `.R` scripts, which are fully reproducible.

| Use case | System | Diversity flavours | Capabilities exercised |
|---|---|---|---|
| **1** | Insectivorous **bat diet** — insect prey ASVs (COI metabarcoding) | neutral, **phylogenetic** | `hilldiv`, `hillprof`, `hilleven`, multi-scale `hillpart` (nested `hierarchy`) |
| **2** | Gut **bacterial MAGs** under an intervention | neutral, phylogenetic, **functional** | `hilldiv` (3 flavours from one call), `traits2dist`, `hillpart` (between-group, per flavour), `hillpair` + ordination (3 spaces), `hillred` |

Together they touch every exported `hill*` function and all three diversity
types.

## Files

```
README.md                       this overview
data-bat-diet.R                 simulates the use-case-1 data set (ASVs, tree, design)
use-case-1-bat-diet.R           use-case-1 analysis + figures
use-case-1-bat-diet.md          use-case-1 manuscript prose
use-case-2-bacterial-mags.R     use-case-2 analysis + figures
use-case-2-bacterial-mags.md    use-case-2 manuscript prose
figures/                        generated PNGs
```

Use case 2 uses the bundled `gut_counts` / `gut_tree` / `gut_traits` data;
use case 1 simulates its own data in `data-bat-diet.R`. Both data sets are
**simulated** to illustrate the analysis paths, not real observations.

## Reproduce

From the package root, with `hilldiv3` and `ggplot2` installed:

```sh
Rscript inst/manuscript/use-case-1-bat-diet.R
Rscript inst/manuscript/use-case-2-bacterial-mags.R
```

Each script writes its figures to `inst/manuscript/figures/` and prints the
summary numbers cited in the corresponding prose.
