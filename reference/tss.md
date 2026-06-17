# Total Sum Scaling normalisation

Normalise a numeric vector or count matrix so that each sample (column)
sums to one. Columns that sum to zero are returned as all-zero (the
`0/0 = NaN` case is mapped to `0`).

## Usage

``` r
tss(abund)
```

## Arguments

- abund:

  A numeric vector or a matrix/data.frame of counts with taxa
  (OTUs/ASVs/MAGs) in rows and samples in columns.

## Value

A normalised object of the same shape as `abund` (vector in, vector out;
matrix/data.frame in, matrix out).

## Examples

``` r
tss(c(a = 1, b = 3))
#>    a    b 
#> 0.25 0.75 
tss(matrix(c(1, 0, 3, 0, 0, 2), nrow = 3))
#>      [,1] [,2]
#> [1,] 0.25    0
#> [2,] 0.00    0
#> [3,] 0.75    1
```
