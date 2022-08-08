portfolio_fun <- function(data) {
  X <- as.matrix(diff(log(data$adjusted))[-1])
  N <- ncol(X)
  T <- nrow(X)
  mu_sm <- colMeans(X)
  Sigma_scm <- cov(X)
  lambdas <- eigen(Sigma_scm)$values
  lmd_mean <- mean(lambdas)
  lmd_max <- max(lambdas)
  JS_t <- rep(sum(solve(Sigma_scm, mu_sm))/sum(solve(Sigma_scm, rep(1, N))), N)
  JS_rho <- (1/T)*(N*lmd_mean - 2*lmd_max)/norm(mu_sm - JS_t, "2")^2
  mu_JS <- (1-JS_rho)*mu_sm + JS_rho*JS_t
  ranked_mu <- sort(mu_JS, decreasing = TRUE, index.return = TRUE)$ix
  w <- rep(0, N)
  w[ranked_mu[1:floor(0.99*N)]] <- 1/floor(0.99*N)
  return(w)
}
