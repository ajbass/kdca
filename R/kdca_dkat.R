#' @import PearsonDS
kdca_dkat <- function(x,
                      y,
                      type_x,
                      mean_adjust,
                      var_adjust,
                      var_adjust_off = FALSE,
                      RINT = FALSE,
                      max.rank = NULL) {

  if (type_x == "continuous") {
    continuous = TRUE
    xk <- xc <- scale(x)
    xkernel <- kernel_linear
  } else {
    continuous = FALSE
    xk <- model.matrix(~1 + as.factor(x))[,-1, drop = F]
    x <- as.numeric(as.factor(x))
    xc <- as.factor(x)
    xkernel = kernel_categorical
  }

  # adjust for mean/variance effects under null data
  if (var_adjust_off) {
    rmod <- apply_regression_dkat(x, xk, y, mean_adjust, var_adjust = NULL, continuous)
  } else {
    rmod <- apply_regression_dkat(x, xk, y, mean_adjust, var_adjust, continuous)
  }
  k <-  ncol(mean_adjust) + 1


  if (RINT) {
    y.scaled <- apply(rmod$y.scaled, 2, RNOmni::RankNorm)
  } else {
    y.scaled <- rmod$y.scaled
  }
  
  cp <- scale(.pairwise_prod(y.scaled), scale = FALSE)
  
  if (is.null(max.rank)) max.rank <- pmin(ncol(cp), nrow(cp))

  # kernel function to measure pair-wise similarity
  Kx <-  xkernel(xc, mod.rank = ncol(xk))
  Ky <- kernel_linear(cp, mod.rank = k, max.rank = max.rank)

  p_lin <- .dkat_speed_cpp(Ux = Kx$vectors,
                           dx =  as.matrix(diag(as.matrix(Kx$values))),
                           Uy = Ky$vectors,
                           dy = as.matrix(diag(Ky$values)))
  levs <-  length(Ky$values)
  Ky$values <- rep(1 / levs, levs)  
  p_proj <- .dkat_speed_cpp(Ux = Kx$vectors,
                            dx =  as.matrix(diag(as.matrix(Kx$values))),
                            Uy = Ky$vectors,
                            dy = as.matrix(diag(Ky$values)))
  Ky <- kernel_gaussian(cp, mod.rank = k, max.rank = max.rank)
  p_gauss <- .dkat_speed_cpp(Ux = Kx$vectors,
                             dx =  as.matrix(diag(as.matrix(Kx$values))),
                             Uy = Ky$vectors,
                             dy = as.matrix(diag(Ky$values)))


  # return p-values
  ret <- data.frame(kernel = c("Linear", "Gaussian", "Projection" ),
                    pvalues = c(p_lin, p_gauss, p_proj))
  return(ret)
}
