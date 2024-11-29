apply_permutation <- function(x,
                              y,
                              type_x,
                              mean_adjust,
                              var_adjust,
                              perm.its,
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
  rmod <- apply_regression(x, xk, y, mean_adjust, var_adjust, continuous)
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
  
  obs.stats.linear <- .dkat_statistic(Ux = Kx$vectors,
                                      dx =  as.matrix(diag(Kx$values)),
                                      Uy = Ky$vectors,
                                      dy = as.matrix(diag(Ky$values)))
  levs <-  length(Ky$values)
  max.rank.lin <- pmin(levs, max.rank)
  Ky$values <- rep(1 / levs, levs)
  obs.stats.projection <- .dkat_statistic(Ux = Kx$vectors,
                                          dx =  as.matrix(diag(Kx$values)),
                                          Uy = Ky$vectors,
                                          dy = as.matrix(diag(Ky$values)))
  Ky <- kernel_gaussian(cp, mod.rank = k, max.rank = max.rank)
  max.rank.gaus <- pmin(length(Ky$values), max.rank)
  obs.stats.gaussian <- .dkat_statistic(Ux = Kx$vectors,
                                        dx =  as.matrix(diag(Kx$values)),
                                        Uy = Ky$vectors,
                                        dy = as.matrix(diag(Ky$values)))
  
  obs.stats <- max(c(obs.stats.gaussian, obs.stats.projection, obs.stats.linear))
  # generate data under the null
  stopit <- FALSE
  i = 1
  null.stats <- null.stats.gauss <- null.stats.linear <- null.stats.proj <- rep(NA, perm.its)
  while (i <= perm.its) {
    # printing
    cat("\r", "Null iteration: ", i)
    if (i == perm.its) cat("\n")
    
    # generate null data set
    id <- sample(1:nrow(y), replace = FALSE)
    y.null <- rmod$y.fit + rmod$y.scaled[id,]  * rmod$w
    
    # adjust for mean/variance effects under null data
    rmod.null <- apply_regression(x, xk, y.null, mean_adjust, var_adjust, continuous)
    
    if (any(is.infinite(rmod.null$y.scaled))) next # dGLM gave poor variance estimate
    
    if (RINT) rmod.null$y.scaled <- apply(rmod.null$y.scaled, 2, RNOmni::RankNorm)
    
    # cross products
    cp.null <- scale(.pairwise_prod(rmod.null$y.scaled), scale = FALSE)
    
    # kernel functions to measure pair-wise similarity
    Ky <- kernel_linear(cp.null, mod.rank = k, max.rank = max.rank.lin)
    null.stats.linear[i]  <- .dkat_statistic(Ux = Kx$vectors,
                                             dx = as.matrix(diag(Kx$values)),
                                             Uy = Ky$vectors,
                                             dy = as.matrix(diag(Ky$values)))
    
    levs <-  length(Ky$values)
    Ky$values <- rep(1 / levs, levs)
    null.stats.proj[i] <- .dkat_statistic(Ux = Kx$vectors,
                                          dx = as.matrix(diag(Kx$values)),
                                          Uy = Ky$vectors,
                                          dy = as.matrix(diag(Ky$values)))
    
    Ky <- kernel_gaussian(cp.null, mod.rank = k, max.rank = max.rank.gaus)
    null.stats.gauss[i] <- .dkat_statistic(Ux = Kx$vectors,
                                           dx = as.matrix(diag(Kx$values)),
                                           Uy = Ky$vectors,
                                           dy = as.matrix(diag(Ky$values)))
    
    i = i + 1
  }
  
  # Empirical p-values
  p_lin <- empirical_pvalues(stat = obs.stats.linear, stat0 = null.stats.linear)
  p_gauss <- empirical_pvalues(stat = obs.stats.gaussian, stat0 = null.stats.gauss)
  p_proj <- empirical_pvalues(stat = obs.stats.projection, stat0 = null.stats.proj)
  
  # Aggregate information across kernels using Fisher's statistic
  lin2 <- sapply(c(obs.stats.linear, null.stats.linear), FUN = function(x) empirical_pvalues(x, null.stats.linear) )
  proj2 <- sapply(c(obs.stats.projection, null.stats.proj), FUN = function(x) empirical_pvalues(x, null.stats.proj) )
  gaus2 <- sapply(c(obs.stats.gaussian, null.stats.gauss), FUN = function(x) empirical_pvalues(x, null.stats.gauss) )
  df2 <- cbind(lin2, proj2, gaus2)
  stat4 <- apply(df2, 1, FUN = function(x) (-2) * sum(log(x)))
  
  p4 <- empirical_pvalues(stat = stat4[1], stat0 = stat4[-1])
  
  # return p-values
  ret <- data.frame(kernel = c("Combined", "Linear", "Gaussian", "Projection" ),
                    pvalues = c(p4, p_lin, p_gauss, p_proj))
  return(ret)
}

apply_permutation_eg <- function(x,
                                 y,
                                 type_x,
                                 mean_adjust,
                                 var_adjust,
                                 perm.its,
                                 return.stats) {

  # check for categorical/continuous
  if (type_x == "continuous") {
    continuous = TRUE
    xk <- xm <- scale(x)
  } else {
    continuous = FALSE
    xm <- as.factor(x)
    xk <- model.matrix(~1 + xm)[,-1, drop = F]
    x <- as.numeric(xm)
  }

  # adjust for mean/variance effects under null data
  rmod <- apply_regression(x, xk, y, mean_adjust, var_adjust, continuous)

  # cross products
  cp <- scale(.pairwise_prod(rmod$y.scaled), scale = FALSE)

  # kernel function to measure pair-wise similarity
  svd.cp <- RSpectra::svds(cp, 1)

  # observed statistic
  PC <- svd.cp$u[, 1]
  obs.stats <- anova(lm(PC ~ 1 + xm))$`F value`[1]

  # generate statistics under the null hypothesis
  null.stats <- rep(NA, perm.its)
  for (i in 1:perm.its) {
    # printing
    cat("\r", "Null iteration: ", i)
    if (i == perm.its) cat("\n")

    # generate null data set
    id <- sample(1:nrow(y), replace = FALSE)
    y.null <- rmod$y.fit + rmod$y.scaled[id,]  * rmod$w

    # adjust for mean/variance effects under null data
    rmod.null <- apply_regression(x, xk, y.null, mean_adjust, var_adjust, continuous)

    if (any(is.infinite(rmod.null$y.scaled))) next

    # cross products
    cp.null <- scale(.pairwise_prod(rmod.null$y.scaled), scale = FALSE)

    # null statistics
    svd.cp <- RSpectra::svds(cp.null, 1)
    PC <- svd.cp$u[, 1]
    null.stats[i] <- anova(lm(PC ~ 1 + xm))$`F value`[1]
  }
  
  # empirical p-value
  ret <- list()
  ret$p <- empirical_pvalues(stat = obs.stats, stat0 = null.stats)

  if (return.stats) {
    ret$null.stats <- null.stats
    ret$obs.stats <- obs.stats
  }

  return(ret)
}
 