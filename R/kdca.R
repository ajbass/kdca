#' Kernel-based Differential Co-expression Analysis
#'
#' @description
#' Given a risk factor of interest, the \code{kdca} function tests for
#' differential co-expression of a pathway. 
#' This function implements three kernels (linear, projection, and Guassian) and 
#' then aggregates information across the kernels to maximize power.
#'
#' @details
#' Please see our pre-print for additional details. 
# The user can specify the type of kernel for inference in the \code{kernel} input, namely, the "linear",
# "projection", or "gaussian" kernel, each of which may be more suitable depending on the data set.
# A user-specified kernel function can also be inputted but following the structure shown in the vignette to avoid
# issues. Currently, the  \code{kernel} input specifies the kernel of  \code{y} and \code{x} is assumed to be the linear kernel.
# A future update will provide an option to choose a kernel function for x.
#'
#' @param x The risk factor(s) of interest (n observations by p variables).
#' @param y A matrix of genes in a pathway (n observations by k genes).
#' @param type_x Specify the variable type for the risk factor x. Either "categorical" or "continuous". If there are multiple risk factors then set to "continuous".
#' @param mean_adjust A matrix of covariates to adjust for mean effects.
#' @param var_adjust A matrix of covariates to adjust on the variance scale.
# @param cor_adjust A matrix of covariates to adjust on the correlation scale.
#' @param perm.its The number of permutations. The default is 1000.
# @param adaptive.its.num The number of permutations to check if the p-values are small enough to continue generating null statistics. The default is \code{NA} which means all the permutations will run and there is no adaptive checking.
# @param adaptive.threshold.num The threshold number of tests when adaptive checking is implemented. If the observed statistic is greater than this threshold then the permutation algorithm continues. Otherwise, it stops early.  
#' @param RINT Apply rank-based inverse normal transform to the residuals. Default is FALSE.
#' @param max.rank Specify rank of kernel matrices. Default is NULL. While setting a low rank can decrease the computational time, the signal may be missed.
#'
#' @return
#' A data frame containing:
#' \item{kernel}{Kernels tested}
#' \item{pvalues}{Empirical p-values.}
#'
#' @examples
#' # set seed
#' set.seed(123)
#'
#' # Generate random categorical risk factor and pathway
#' x <- matrix(rbinom(50, size = 1, prob = 0.5), ncol = 1)
#' y <- matrix(rnorm(50 * 10), ncol = 10)
#'
#' kdca.out <- kdca(x, y, type_x = "categorical")
#' 
#' @export
kdca <- function(x,
                 y,
                 type_x = c("categorical", "continuous"),
                 mean_adjust = matrix(rep(1, nrow(x), ncol = 1)),
                 var_adjust = matrix(rep(1, nrow(x), ncol = 1)),
                 perm.its = 1000,
                 RINT = FALSE,
                 max.rank = NULL) {
  
  # input checks
  if (!is.matrix(x) & !is.matrix(y)) stop("Inputs x, y must be matrices")
  if (!is.matrix(mean_adjust) & !is.matrix(var_adjust)) stop("Input covariates must be matrices")
  if (nrow(x) != nrow(y)) stop("Inputs x, y must have the same number of observations")
  if ((nrow(mean_adjust) != nrow(y)) & (nrow(var_adjust) != nrow(y))) stop("Inputs mean_adjust and var_adjust must have the same number of observations")
  if (ncol(y) <= 1) stop("Pathway size must be greater than 1")
  
  # check for intercept in var adjustments
  s <- matrixStats::colSds(var_adjust)
  if (!any(s < 1e-14)) {
    # add intercept
    var_adjust = cbind(1, var_adjust)
  }
  
  # check for intercept in mean adjustments
  m <- matrixStats::colSds(mean_adjust)
  if (!any(m < 1e-14)) {
    # add intercept
    mean_adjust = cbind(1, mean_adjust)
  }
 
  out <- apply_permutation(x = x,
                           y = y,
                           type_x = type_x,
                           mean_adjust =  mean_adjust,
                           var_adjust = var_adjust,
                           perm.its = perm.its,
                           RINT = RINT,
                           max.rank = max.rank)
 
  return(out)
}

# internal use
run_dkat <- function(x,
                     y,
                     type_x = c("categorical", "continuous"),
                     mean_adjust = matrix(rep(1, nrow(x), ncol = 1)),
                     var_adjust = matrix(rep(1, nrow(x), ncol = 1)),
                     perm.its = 1000,
                     RINT = FALSE,
                     max.rank = NULL) {
  
  # input checks
  if (!is.matrix(x) & !is.matrix(y)) stop("Inputs x, y must be matrices")
  if (!is.matrix(mean_adjust) & !is.matrix(var_adjust)) stop("Input covariates must be matrices")
  if (nrow(x) != nrow(y)) stop("Inputs x, y must have the same number of observations")
  if ((nrow(mean_adjust) != nrow(y)) & (nrow(var_adjust) != nrow(y))) stop("Inputs mean_adjust and var_adjust must have the same number of observations")
  if (ncol(y) <= 1) stop("Gene set size must be greater than 1")
  
  # check for intercept in var adjustments
  s <- matrixStats::colSds(var_adjust)
  if (!any(s < 1e-14)) {
    # add intercept
    var_adjust = cbind(1, var_adjust)
  }
  
  # check for intercept in mean adjustments
  m <- matrixStats::colSds(mean_adjust)
  if (!any(m < 1e-14)) {
    # add intercept
    mean_adjust = cbind(1, mean_adjust)
  }

  out <- kdca_dkat(x = x,
                   y = y,
                   type_x = type_x,
                   mean_adjust = mean_adjust,
                   var_adjust = var_adjust,
                   RINT = RINT,
                   max.rank = max.rank)

  return(out)
}

#' @import stats
#' @import dglm
#' @import RSpectra
#' @importFrom Rcpp sourceCpp
#' @importFrom matrixStats colSds
#' @useDynLib kdca
NULL
