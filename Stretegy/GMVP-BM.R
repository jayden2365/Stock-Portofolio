portfolio_fun <- function(data) {
  
  # Form Code # using T instead of T-1
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
  
  GMVP <- function(Sigma) {
    ones <- rep(1, nrow(Sigma))
    Sigma_inv_1 <- solve(Sigma, ones)  #same as: inv(Sigma) %*% ones
    w <- (1/as.numeric(ones %*% Sigma_inv_1)) * Sigma_inv_1
    return(w)
  }
  
  ## Run Code and paramaters
  X <- as.matrix(diff(log(data$adjusted))[-1])
  w <- GMVP(cov_LedoitWolf(X))
  return(w)
}
