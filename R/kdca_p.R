empirical_pvalues <- function(stat, stat0) {
  r <- sum(stat0 >= stat, na.rm = T)
  (r + 1) / (length(stat0) + 1)
}
