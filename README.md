
<!-- README.md is generated from README.Rmd. Please edit that file -->

# flipped

The goal of flipped is to explore models for coin flippings beyond the
standard binomial.

## Installation

You can install the development version of this package from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("bomeara/flipped")
```

## Example

This is a basic example:

``` r
library(flipped)
nflips <- 100
pheads <- 0.3
nheads <- stats::rbinom(n=1, size=nflips, prob=pheads)
```
