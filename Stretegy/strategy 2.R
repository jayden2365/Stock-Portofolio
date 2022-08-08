portfolio_fun <- function(data) {
  X <- as.matrix(diff(log(data$adjusted))[-1])
  mu <- colMeans(X)
  Sigma <- cov(X)
  sigma <- sqrt(diag(Sigma))
  w <- 1/sigma
  w <- w/sum(w)
  return(w)
}
