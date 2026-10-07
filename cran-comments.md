## Initial CRAN release (3.0.0)

These are the submission notes for version 3.0.0, published on CRAN on
6 October 2026. They should be updated for a later submission.

## Resubmission

This is a resubmission. The previous source archive was built more than a
month before it was uploaded, which caused the incoming-feasibility NOTE
"This build time stamp is over a month old." The archive for this submission
has been rebuilt immediately before upload.

In response to the earlier CRAN reviewer's comments:

* The redundant "Tools for" at the start of the 'Description' field has
  been removed.

* All acronyms in the 'Description' field are now written out on first
  use: operational taxonomic units (OTUs), amplicon sequence variants
  (ASVs), metagenome-assembled genomes (MAGs) and Faith's phylogenetic
  diversity (PD).

* Regarding suggested packages not available from a mainstream
  repository: the only non-CRAN packages in 'Suggests' are the three
  Bioconductor packages used by optional input adapters ('phyloseq',
  'SummarizedExperiment', 'TreeSummarizedExperiment'). These are
  available from the standard Bioconductor software repository, which
  the incoming checks already treat as a mainstream repository, so the
  previous 'Additional_repositories' entry matched no dependency and was
  reported as an unused entry ("?  ?  <URL>") in the availability table.
  The field has therefore been removed, and the three packages are
  instead named in the 'Description' field together with the
  Bioconductor URL. Every use of them is guarded with
  `rlang::check_installed()` / `requireNamespace()`, so the package
  builds, checks and runs without them.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Notes

* The check reports a NOTE for "New submission". This is the first
  submission of 'hilldiv3' to CRAN.

* The same NOTE may flag possibly misspelled words in DESCRIPTION
  ('ASVs', 'MAGs', 'OTUs', 'Rao', 'Sorensen', 'UniFrac', 'Jost', 'Chao',
  'Chiu', 'Alberdi', 'et', 'al', and the package names 'phyloseq',
  'SummarizedExperiment', 'TreeSummarizedExperiment'). These are domain
  terms, package names, author surnames and the abbreviation "et al." in
  the cited references; they are spelled correctly.

## Downstream dependencies

* There are currently no downstream dependencies (new package).

## Test environments

* local macOS, R 4.3.3 -- 0 errors | 0 warnings | 1 note
* win-builder, R-devel (2026-06-26 r90195 ucrt) -- 0 errors | 0 warnings |
  1 note (the "New submission" NOTE described above)
