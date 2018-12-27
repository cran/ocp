<!-- README.md is generated from README.Rmd. Please edit that file -->
ocp
===

The goal of ocp is to implement Bayesian Online Changepoint Detection, as described: <https://arxiv.org/abs/0710.3742>

Example
-------

This is a basic example of how to use the function "onlineCPD" on simulated univariate Gaussian data as input.

``` r
library(ocp)
 # the true changepoint locations including the first and last point
truecps<- c(1, 51, 71, 121)
#simulate the data
set.seed(1)
uvg<- c(rnorm(n=diff(truecps)[1], mean=0, sd=2), 
        rnorm(n=diff(truecps)[2], mean=20, sd=4),
        rnorm(n=diff(truecps)[3], mean=10, sd=3))
ocpd_output<- onlineCPD(uvg) 
```
