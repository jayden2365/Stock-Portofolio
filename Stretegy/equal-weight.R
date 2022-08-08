portfolio_fun <- function(data) {
  X <- as.matrix(diff(log(data$adjusted))[-1])  # compute log returns
  mu <- colMeans(X)  # compute mean vector
  Sigma <- cov(X)  # compute the SCM
  # design mean-variance portfolio
  i <- sort(mu, decreasing = TRUE, index.return = TRUE)$ix
  N <- length(mu) - 1
  w <- rep(0, N)
  w[i[1:N]] <- 1/N
  names(w) <- colnames(X)
  return(w)
}