
<!-- README.md is generated from README.Rmd. Please edit that file -->

# funestim

<!-- badges: start -->
<!-- badges: end -->

`funestim` is an `R`-package that allows users to estimate the mean and
the covariance of a functional dataset using an estimation of the
regularity of the curves. The curves can be irregularly sampled.

## Installation

To install the latest version directly from
[GitHub](https://github.com/), please use

``` r
# install.packages("devtools")
devtools::install_github("StevenGolovkine/funestim")
```

To build the vignette as well, please use

``` r
# install.packages("devtools")
devtools::install_github("StevenGolovkine/funestim", build_vignettes = TRUE)
```

## Dependencies

The `funestim` package depends on the `R`-packages
[`dplyr`](https://CRAN.R-project.org/package=dplyr),
[`fdapace`](https://CRAN.R-project.org/package=fdapace),
[`funData`](https://CRAN.R-project.org/package=funData),
[`gss`](https://CRAN.R-project.org/package=gss),
[`KernSmooth`](https://CRAN.R-project.org/package=KernSmooth),
[`magrittr`](https://CRAN.R-project.org/package=magrittr),
[`purrr`](https://CRAN.R-project.org/package=purrr),
[`Rcpp`](https://CRAN.R-project.org/package=Rcpp) and
[`RcppArmadillo`](https://CRAN.R-project.org/package=RcppArmadillo),
[`tidyr`](https://CRAN.R-project.org/package=tidyr),

## References

The theoretical foundations of the estimation of the mean and covariance
functions are described in:

Golovkine S., Klutchnikoff N., Patilea V. (2021) - Adaptive estimation
of irregular mean and covariance functions.

## Bug reports

Please use [GitHub
issues](https://github.com/StevenGolovkine/simulater/issues) for
reporting bugs or issues.
