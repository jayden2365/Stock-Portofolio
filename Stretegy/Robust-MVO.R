library(CVXR)
library(ICSNP)

portfolio_fun <- function(data) {
  
  portfolioMarkowitzRobust <- function(mu_hat, Sigma_hat, kappa, delta_, lmd = 0.5) {
    N <- length(mu_hat)
    S12 <- chol(Sigma_hat)  # t(S12) %*% S12 = Sigma
    w <- Variable(N)
    prob <- Problem(Maximize(t(w) %*% mu_hat - kappa*norm2(S12 %*% w) 
                             - lmd*(norm2(S12 %*% w) + delta_*norm2(w))^2),
                    constraints = list(w >= 0, sum(w) == 1))
    result <- solve(prob)
    return(as.vector(result$getValue(w)))
  }
  
  # Ref from Code # using T instead of T-1 # shrinkage
  cov_LedoitWolf <- function(X) { 
    T <- nrow(X) 
    N <- ncol(X) 
    T_eff <- T # T for biased SCM and T-1 for unbiased 
    
    Xc <- X - matrix(colMeans(X), T, N, byrow = TRUE) # center just in case (it does not destroy anything if it was already centered more sophisticately) 
    S <- crossprod(Xc)/(T_eff-1) # SCM #### T-1
    Sigma_T <- mean(diag(S)) * diag(N) # target 
    d2 <- sum((S - Sigma_T)^2) 
    b2_ <- (1/T^2) * sum(rowSums(Xc^2)^2) - (2*T_eff/T - 1)/T * sum(S^2) 
    b2 <- min(b2_, d2) #a2 = d2 - b2 
    rho <- b2/d2 
    Sigma_clean <- (1-rho)*S + rho*Sigma_T #rho=b2/d2, 1-rho=a2/d2 
    return(Sigma_clean) 
  }
  
  ## Run Code and paramaters
  X <- as.matrix(diff(log(data$adjusted))[-1])
  T <- nrow(X)
  sigma <- cov_LedoitWolf(X)
  
  # multiple robust solutions
  kappa <- 1.0
  delta <- 0.1*sqrt(T)
  mu <- ICSNP::spatial.median(X)
  w <- portfolioMarkowitzRobust(mu, sigma, kappa, delta/sqrt(T-1))
  return(w)
}
