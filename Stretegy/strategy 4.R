library(ICSNP)

portfolio_fun <- function(data) {
  
  ## Support function
  # Tyler estimator with shrinkage
  estimateTylerShrinkage <- function(X, R_target, alpha) {
    max_iter <- 100
    error_th_Sigma <- 1e-3
    
    # Gaussian initial point
    Sigma <- 1/(1+alpha)*cov(X) + alpha/(1+alpha)*R_target 
    Sigma <- Sigma/sum(diag(Sigma))
    # Loop
    obj_value_record <- Sigma_diff_record <- rep(NA, max_iter)
    for (k in 1:max_iter) {
      Sigma_prev <- Sigma
      
      #Tyler update
      weights <- 1/rowSums(X * (X %*% solve(Sigma)))   # 1/diag( X %*% inv(Sigma) %*% t(X) )
      obj_value_record[k] <- - (N/2)*sum(log(weights)) + (T/2)*sum(log(eigen(Sigma)$values))
      Sigma <- (N/T) * crossprod( sqrt(weights)*X )  # (N/T) * t(X) %*% diag(weights) %*% X
      Sigma <- 1/(1+alpha)*Sigma + alpha/(1+alpha)*R_target 
      Sigma <- Sigma/sum(diag(Sigma))
      
      #stopping criterion
      Sigma_diff_record[k] <- norm(Sigma - Sigma_prev, "F")/norm(Sigma_prev, "F")
      if (Sigma_diff_record[k] < error_th_Sigma)
        break
    }
    obj_value_record <- obj_value_record[1:k]
    Sigma_diff_record <- Sigma_diff_record[1:k]
    
    # Recover missing scaling factor
    sigma2 <- apply(X, 2, var)
    d <- diag(Sigma)
    kappa <- sum(sigma2*d)/sum(d*d)
    Sigma <- kappa*Sigma
    
    return(Sigma)
  }
  
  
  ## Run Code and paramaters
  X <- as.matrix(diff(log(data$adjusted))[-1])
  N <- ncol(X)
  T <- nrow(X)
  
  ## Median and upadate the data
  mu_gmedian <- ICSNP::spatial.median(X)
  X_ <- X - rep(mu_gmedian, each = T)
  
  ## shrinkage
  R_target <- mean(diag(cov(X))) * diag(N)  # shrinkage target
  alpha <- min(0.01*N/T, sqrt(2)/2)  #shrinkage factor
  
  ## portfolio construction
  R_Tyler_shrinked <- estimateTylerShrinkage(X_, R_target, alpha)
  sigma <- sqrt(diag(R_Tyler_shrinked))
  w <- 1/sigma
  w <- w/sum(w)
  return(w)
}
