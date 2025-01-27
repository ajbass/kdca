# kdca package

<!-- badges: start -->

[![R-CMD-check](https://github.com/ajbass/kdca/actions/workflows/R-CMD-CHECK.yml/badge.svg)](https://github.com/ajbass/kdca/actions/workflows/R-CMD-CHECK.yml)

<!---[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/kdca)](https://cran.r-project.org/package=kdca)--->
<!-- badges: end -->

## Overview

The `kdca` package implements a statistical test to assess whether a pathway is differentially co-expressed. The inputs into `kdca` are pathway expression values and a risk factor(s) of interest. The package can model multiple risk factors, covariates, and/or batch effects. `kdca` also accounts for variance-specific effects from the risk factors (or other variables) to avoid false discoveries. The risk factors can be continuous, categorical, or a mixture of both. Finally, we assume an appropriate normalization has been applied to the expression values (e.g., logCPM). For RNA-seq data, we also recommend treating library size as a covariate in the model (see below).

`kdca`  captures complex differential co-expression using the linear, projection, and Gaussian kernels. Because pathway architecture is unknown, it is difficult to select the optimal kernel.  As such, we provide an aggregate version (`Combined') that maximizes power across kernels. In future updates, users will be able to add their own kernel functions. 

## Citing this package 

The methods implemented in this package are described in:

> Bass AJ, Cutler DJ, Epstein MP. A powerful framework for differential co-expression analysis of continuous and categorical risk factors. Submitted; 2024.

## Getting help

To report any bugs or issues related to usage please report it on GitHub at https://github.com/ajbass/kdca.

## Installation

To install the development version of the package:

```{r, eval = FALSE}
# install development version of package
install.packages("devtools")
library("devtools")
devtools::install_github("ajbass/kdca")
```

## Quick start

We demonstrate the functionality of the `kdca` package using simulated data with no differential co-expression. First, load the `kdca` package and simulate data:

```{r}
library(kdca)
# set seed
set.seed(123)

# Generate a categorical risk factor (sample size 100)
x <- matrix(rbinom(100, size = 1, prob = 0.5), ncol = 1)

# Generate a continuous covariate
z <-  matrix(rnorm(100 * 1), ncol = 1)

# Generate a pathway with 10 genes (no differential co-expression)
# with a variance-specific effect
ve <- rep(x + 1, 5)
error <- matrix(rnorm(100 * 10, sd = ve), ncol = 10)
beta <- matrix(1, ncol = 10) # effect sizes
y <- z %*% beta + x %*% beta + error # pathway
```

We generated a risk factor `x` that is categorical, a covariate `z` that is continuous, and a pathway with 10 genes (100 samples). The risk factor and covariate impact the expression values but there is no differential co-expression. Additionally, the risk factor has a strong variance-specific effect. We can then run the `kdca` function on this data set:

```{r, eval = F}
kdca_out <- kdca(x, y, mean_adjust = z, type_x = "categorical")
```

Note that the argument `type_x` specifies the risk factor type and the covariates/batch effects can be inputted into the argument `mean_adjust`.

We can apply `kdca` to conituous risk factors as well. Let us now suppose that `z` is the risk factor of interest and `x` is a covariate:

```{r, eval = F}
kdca_out <- kdca(z, y, mean_adjust = x, type_x = "continuous")
```
Finally, if `x` and `z` are risk factors of interest, then we can apply `kdca` to multiple risk factors:

```{r, eval = F}
kdca_out <- kdca(cbind(x,z), y, type_x = "continuous")
```

When inputting multiple risk factors, you have to set `type_x = "continuous"`. See `?kdca` for additional arguments.
