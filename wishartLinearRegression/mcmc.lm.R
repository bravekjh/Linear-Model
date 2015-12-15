##
## Wishart Multivariate Response Linear Regression
##
## John Tipton - created 12.214.2015
##

##
## libraries and functions
##

mcmc <- function(Y, X, params){
  n_mcmc <- params$n_mcmc
  mu_beta <- params$mu_beta
  s2_beta <- params$s2_beta
  s2_lower <- params$s2_lower
  s2_lower <- params$s2_upper
  nu <- params$nu
  S <- params$S
    
  ##
  ## Initialize variables
  ##
  
  n <- dim(Y)[1]
  d <- dim(Y)[2]
  p <- dim(X)[2]

  s <- runif(s2_lower, s2_upper)
  s2 <- s^2
  Q_inv <- solve(matrix(rWishart(1, nu, S), d, d))

  ##
  ## save variables
  ##
  
  beta_save <- array(0, dim=c(n_mcmc, p, d))
  s2_save <- rep(0, n_mcmc)
  Q_inv_save <- array(0, dim=c(n_mcmc, d, d))

  ##
  ## Start MCMC
  ##
  
  for(k in 1:n_mcmc){
    if(k %% 100 == 0){
      cat(" ", k)
    }
    
    ##
    ## sample beta
    ##
#     for (j in 1:d) {
#       
#     }
    
    ##
    ## sample mu.beta
    ##
    
    A.chol <- chol(Sigma.beta.inv + Sigma.0.inv)
    b <- (Sigma.beta.inv %*% beta + Sigma.0.inv %*% mu.0)
    mu.beta <- rMVN(A.chol, b)
    
    ##
    ## sample sigma.squared.beta
    ##
    
    sigma.squared.beta <- 1 / rgamma(1, alpha.beta + tau / 2, beta.beta + 1 / 2 * t(beta - mu.beta) %*% (beta - mu.beta))
    Sigma.beta <- sigma.squared.beta * I.beta
    Sigma.beta.inv <- 1 / sigma.squared.beta * I.beta
    
    ##
    ## sample sigma.squared.epsilon
    ##
    
    sigma.squared.epsilon <- 1 / rgamma(1, alpha.epsilon + n / 2, beta.epsilon + 1 / 2 * t(Y - X %*% beta) %*% (Y - X %*% beta))
    Sigma.epsilon <- sigma.squared.epsilon * I.epsilon
    Sigma.epsilon.inv <- 1 / sigma.squared.epsilon* I.epsilon
    
    ##
    ## DIC calculations
    ##
    
    Dbar.save[k] <- - 2 * sum(dnorm(Y, X %*% beta, sqrt(sigma.squared.epsilon), log = TRUE))
    
    ##
    ## save variables
    ##
    
    beta.save[, k] <- beta
    sigma.squared.epsilon.save[k] <- sigma.squared.epsilon
    sigma.squared.beta.save[k] <- sigma.squared.beta
  }
  
  ##
  ## output
  ##
  
  list(beta.save = beta.save, sigma.squared.beta.save = sigma.squared.beta.save, sigma.squared.epsilon.save = sigma.squared.epsilon.save)
}