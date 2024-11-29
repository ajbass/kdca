apply_regression <- function(x, xk, y, mean_adjust, var_adjust, continuous) {
  scaled <- w <- y.fit <- y
  
  # Continuous risk factors (or multiple): double GLM
  if (continuous | ncol(var_adjust) > 1) {
    for (i in 1:ncol(y)) {
      rmod <- double_GLM(x, y[, i], mean_adjust, var_adjust)
      if (is.null(rmod$dispersion.fit)) {
        w[, i] <- rep(sd(residuals(rmod)), nrow(w))
      } else {
        w[, i] <- sqrt(predict.glm(rmod$dispersion.fit, type = 'response'))
      }
      h <- hatvalues(rmod)
      y.fit[, i] <- fitted(rmod)
      scaled[, i] <- residuals(rmod) / sqrt(1 - h) / w[, i]
    }

    ret <- list(y.scaled = scaled,
                y.fit = y.fit,
                w = w)
  } else {
    # categorical risk factors: standardize by group
    rmod <- .quick_lm_boot_cpp(cbind(mean_adjust, xk), y) 
    res <- rmod$res 
    h <- rmod$hatvals 
    grps <- unique(x)
    for (i in 1:length(grps)) {
      id <- x == grps[i]
      res_grp <- res[id,]
      h_grp <- h[id]
      tmp <- res_grp / sqrt(1 - h_grp) 
      df.res <- ncol(mean_adjust) + 1
      df.tot <- length(x)
      obs.df <-  (df.tot - df.res) / df.tot  
      sd_grp <-  sqrt(colSums((tmp) ^ 2) / (obs.df * sum(id)))
      res_grp <-  t(t(tmp) / sd_grp)
      scaled[id,] <- res_grp
      w[id,] <- rep(sd_grp, each=sum(id)) 
    }
 
    ret <- list(y.scaled = scaled,
                y.fit = rmod$fitted, 
                w = w)
  }
  return(ret)
}

double_GLM <- function(x, y, mean_adjust, var_adjust) {
  tryCatch(dglm(y ~ -1 + x + mean_adjust, ~ -1 + x + var_adjust, family = gaussian),
           warning = function(e) {
             error = lm(y ~ -1 + x + mean_adjust)},
           error = function(e) {
             error = lm(y ~ -1 + x + mean_adjust)})
}

apply_regression_dkat <- function(x, xk, y, mean_adjust, var_adjust, continuous) {
  scaled <- w <- y.fit <- y

  nc <- ifelse(!is.null(var_adjust), ncol(var_adjust), 0)
  if (continuous | nc > 1) {
    for (i in 1:ncol(y)) {
      rmod <- double_GLM2(x, y[, i], mean_adjust, var_adjust)
      if (is.null(rmod$dispersion.fit)) {
        w[, i] <- rep(sd(residuals(rmod)), nrow(w))
      } else {
        w[, i] <- sqrt(predict.glm(rmod$dispersion.fit, type = 'response'))
      }
      h <- hatvalues(rmod)
      y.fit[, i] <- fitted(rmod)
      scaled[, i] <- residuals(rmod) 
    }

    ret <- list(y.scaled = scaled,
                y.fit = y.fit,
                w = w)
  } else {
    rmod <- .quick_lm_boot_cpp(cbind(mean_adjust, xk), y) #
    res <- rmod$res 
    h <- rmod$hatvals 
    grps <- unique(x)
    for (i in 1:length(grps)) {
      id <- x == grps[i]
      res_grp <- res[id,]
      h_grp <- h[id]
      tmp <- res_grp / sqrt(1 - h_grp) 
      df.res <- ncol(mean_adjust) + 1
      df.tot <- length(x)
      obs.df <-  (df.tot - df.res) / df.tot  
      sd_grp <-  sqrt(colSums((tmp) ^ 2) / (obs.df * sum(id)))
      res_grp <- if (!is.null(var_adjust)) {
      res_grp = t(t(tmp) / sd_grp)
      } else {
        res_grp = tmp
      }
      scaled[id,] <- res_grp
      w[id,] <- rep(sd_grp, each = sum(id)) 
    } 
    ret <- list(y.scaled = scaled,
                y.fit = rmod$fitted, 
                w = w)
  }
  return(ret)
}

double_GLM2 <- function(x, y, mean_adjust, var_adjust) {
  if (!is.null(var_adjust)) {
    tryCatch(dglm(y ~ -1 + x + mean_adjust, ~ -1 + x + var_adjust, family = gaussian),
             warning = function(e) {
               error = lm(y ~ -1 + x + mean_adjust)},
             error = function(e) {
               error = lm(y ~ -1 + x + mean_adjust)})
  } else {
    lm(y ~ -1 + x + mean_adjust)
  }
}
