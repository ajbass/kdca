#' Inference of differential co-expression using the eigengene approach
#'
#' @description
#' Given a risk factor of interest, the \code{eigengene} function tests for
#' differential co-expression of a pathway using the top eigengene of the correlation matrix.
#'
#' @param x The primary variable(s) of interest (n observations by p variables)
#' @param y A matrix of genes  in a gene set (n observations by k genes)
#' @param type_x Specify type of primary variable. "categorical" or "continuous.
#' @param mean_adjust A matrix of covariates to adjust for mean effects.
#' @param var_adjust A matrix of covariates to adjust for variance effects.
#' @param perm.its The number of permutations. The default is 1000.
#' @param return.stats If TRUE then the permutation statistics are returned. The default is FALSE.
#'
#' @return
#' A list of object type "kdca" containing:
#' \item{call}{Function call.}
#' \item{pvalues}{An empirical p-value.}
#' \item{null.stats}{If return.stats is TRUE, then a vector of null statistics is returned.}
#' \item{obs.stats}{If return.stats is TRUE, then the observed test statistic is returned.}
#'
#' @examples
#' # set seed
#' set.seed(123)
#'
#' # Generate random study design and gene set
#' x <- matrix(rbinom(50, size = 1, prob = 0.5), ncol = 1)
#' y <- matrix(rnorm(50 * 10), ncol = 10)
#'
#' eg.out <- eigengene(x, y, type_x = "categorical")
#' @export
eigengene <- function(x,
                      y,
                      type_x,
                      mean_adjust = matrix(rep(1, length(x), ncol = 1)),
                      var_adjust = matrix(rep(1, length(x), ncol = 1)),
                      perm.its = 1000,
                      return.stats = FALSE) {

  # input checks
  if (!is.matrix(x) & !is.matrix(y)) stop("Inputs x, y must be matrices")
  if (!is.matrix(mean_adjust) & !is.matrix(var_adjust)) stop("Input covariates must be matrices")
  if (nrow(x) != nrow(y)) stop("Inputs x, y must have the same number of observations")
  if ((nrow(mean_adjust) != nrow(y)) & (nrow(var_adjust) != nrow(y)) ) stop("Inputs mean_adjust and var_adjust must have the same number of observations")
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

  # Run permutation algorithm
  out <- apply_permutation_eg(x = x,
                              y = y,
                              type_x = type_x,
                              mean_adjust = mean_adjust,
                              var_adjust = var_adjust,
                              perm.its = perm.its,
                              return.stats = return.stats)

  if (return.stats) {
    retval <- list(call = match.call(),
                   pvalues =  out$p,
                   null.stats = out$null.stats,
                   obs.stats = out$obs.stats)
  } else {
    retval <- list(call = match.call(),
                   pvalues = out$p)
  }

  class(retval) <- "eigengene"
  return(retval)
}
