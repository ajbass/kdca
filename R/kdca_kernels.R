kernel_categorical <- function(x, mod.rank = 2) {
  m <- length(x)
  K = outer(x, x, function(x, y) as.integer(x == y))
  K <- t(t(K - colSums(K) / m) -  rowSums(K) / m) + sum(K) / m ^ 2

  K <- RSpectra::eigs_sym(K, k = length(unique(x)) + 1)
  id <- K$values > 1e-10
  K$values <- K$values[id] / sum(K$values[id])
  K$vectors <- K$vectors[, id, drop = F]
  return(K)
}

kernel_gaussian <- function(x, mod.rank = 2, max.rank = 5000) {
  s = 1/10000
  rbf <- kernlab::rbfdot(sigma =  s)
  K <- kernlab::kernelMatrix(rbf, x)
  m <- nrow(x)
  K <- t(t(K - colSums(K) / m) -  rowSums(K) / m) + sum(K) / m ^ 2
  K <- RSpectra::eigs_sym(K, k = pmin(ncol(x), pmin(max.rank,nrow(x) - mod.rank)))
  id <- K$values > 1e-10
  K$values <- K$values[id] / sum(K$values[id])
  K$vectors <- K$vectors[, id, drop = F]
  return(K)
}

kernel_linear <- function(x, mod.rank = 2, max.rank = 5000, svd.scale = TRUE) {
  p <- ncol(x)
  k = pmin(p, nrow(x) - mod.rank)
  svd.out <- .quick_svd(x)
  k <- pmin(k, max.rank, length(svd.out$evalues))
  evals <- svd.out$evalues
  K <- list(vectors = svd.out$evectors[, 1:k, drop = F],
            values = as.numeric(evals)[1:k] / sum(evals[1:k]))
  return(K)
}

kernel_projection <- function(x, mod.rank = 2, max.rank = 5000) {
  p <- ncol(x)
  k = pmin(p, nrow(x) - mod.rank)
  svd.out <- .quick_svd(x)
  k <- pmin(k, max.rank, length(svd.out$evalues))
  evals <- rep(1, k)
  svd.out$evalues <- as.numeric(evals) / sum(evals)
  K <- list(vectors = svd.out$evectors[, 1:k, drop = F], values = svd.out$evalues)
  return(K)
}
