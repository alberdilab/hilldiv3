# Tidy result objects for hilldiv3

Internal helpers that turn the engine's wide matrices into the
long-format (tidy) data frames returned by the user-facing `hill*`
functions, and the S3 [`print()`](https://rdrr.io/r/base/print.html) /
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) / `autoplot()`
methods for those results. Every result carries a common parent class
`hill_result` plus a specific subclass, so the shared machinery
(printing, line plots) is written once.
