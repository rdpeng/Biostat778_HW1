Biostat778_HW1
==============

Homework 1 for Biostat 778

Due: November 13, 2013

## Fast Multivariate Normal Density

Write a function called dmvnorm that evaluates the multivariate Normal
density with mean $\mu$ and covariance $S$. The density function is

Assume the the covariance matrix $S$ is always symmetric and
full-rank. The function should have the form

```r
dmvnorm <- function(x, mu, S, log = TRUE) {
        ## Your code here
}
```

where $x$ is a $n\times k$ matrix of points to be evaluated, $\mu$ is
a $k \times 1$ vector of means for the $k$-dimensional Normal, and $S$
is a $k\times k$ covariance matrix. Your function should be written
using R code only and you MAY NOT use any additional packages beyond
the standard base packages that come with R.


## Fast Linear Regression
