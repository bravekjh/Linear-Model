##
## Simple linear regression model 
##
## John Tipton - created 01.25.2014
##

##
## model: y = X %*% beta + epsilon
##

mcmc <- function(y, X, H, n_mcmc, alpha_epsilon, beta_epsilon, alpha_lambda, beta_lambda, pca = FALSE){
  
  ##
  ## libraries and functions
  ##
  
  library(myFunctions)
  
  ##
  ## Initialize variables
  ##
  
  n <- length(y)
  tau <- dim(X)[2]
  if(pca == TRUE){
    X_pca = makePCA(X)
    Lambda_inv = diag(1/ as.vector(X_pca$sdev^2))
    X = X_pca$X_pca 
  } else {
    Lambda_inv = diag(tau)
  }
  
  # s2_epsilon
  s2_epsilon <- max(10, 1 / rgamma(1, alpha_epsilon, beta_epsilon))
  sqrt_s2_epsilon = sqrt(s2_epsilon)
  # lambda
  lambda2 <- rgamma(1, alpha_lambda, beta_lambda)
  # gamma
  gamma2_inv <- 1 / rgamma(tau, 1, lambda2 / 2)
  D_gamma_inv <- diag(gamma2_inv)
  D_gamma_Lambda_inv = D_gamma_inv %*% Lambda_inv
  # beta
  ybar = mean(y)
  J_mu = rep(1, length(y))
  mu <- J_mu * ybar
  y_minus_mu = y - mu
  mu_beta <- rep(0, tau)
  beta <- as.vector(mvrnormArma(1, mu_beta, s2_epsilon * solve(D_gamma_inv)))
  
  ## variables for speed
  HX = X[H, ]
  tHX = t(X[H, ])
  tHXX = tHX %*% X[H, ]
  HXbeta = HX %*% beta
  sqrt_n = sqrt(n)
  
  ##
  ## save variables
  ##
  
  n_burn=n_mcmc / 2
  n_save = n_mcmc - n_burn
  mu_save = rep(0, n_save)
  beta_save <- matrix(nrow = tau, ncol = n_save)
  s2_epsilon_save <- vector(length = n_save)
  gamma2_save <-  matrix(nrow = tau, ncol = n_save)
  lambda2_save <- vector(length = n_save)
  Dbar_save <- vector(length = n_save)
  
  ##
  ## Start MCMC
  ##
  
  for(k in 1:n_mcmc){
    if(k %% 100 == 0){
      cat(" ", k)
    }
    
    ##
    ## sample mu
    ##
    
    a = 1 / s2_epsilon * n
    b = 1 / s2_epsilon * sum(y - HXbeta)
    mu = rMVNArmaScalar(a, b)
    y_minus_mu = y - mu
    
    ##
    ## sample beta
    ##
    
    A <- 1 / s2_epsilon * (tHXX + D_gamma_Lambda_inv) ## check if I need variance term here
    b <- 1 / s2_epsilon * (tHX %*% y_minus_mu) ## check if I need variance term here
    beta <- rMVNArma(A, b)
    HXbeta = HX %*% beta
    
    ##
    ## sample s2_epsilon
    ##
    
    deviance = y_minus_mu - HXbeta
    s2_epsilon <- 1 / rgamma(1, alpha_epsilon + n / 2 + tau / 2, beta_epsilon + 1 / 2 * t(deviance) %*% deviance + 1 / 2 * t(beta) %*% D_gamma_Lambda_inv %*% beta)
    sqrt_s2_epsilon = sqrt(s2_epsilon)
    
    ##
    ## sample gamma2
    ##
    
    mu_tilde <- sqrt(lambda2 * s2_epsilon / beta^2)
    lambda_tilde <- lambda2
    gamma2_inv <- rinvgauss(tau, mu_tilde, lambda_tilde)
    D_gamma_inv <- diag(as.vector(gamma2_inv))
    D_gamma_Lambda_inv = D_gamma_inv %*% Lambda_inv
    gamma2 <- 1 / gamma2_inv
    
    ##
    ## lambda2
    ##
    
    lambda2 <- rgamma(1, alpha_lambda + tau, beta_lambda + sum(gamma2) / 2)
    
    ##
    ## save variables
    ##
    if(k > n_burn){
      k_tmp = k - n_burn
      ## DIC calculations
      Dbar_save[k_tmp] <- - 2 * sum(dnorm(y, mu + HXbeta, sqrt_s2_epsilon, log = TRUE))
      mu_save[k_tmp] <- mu[1]
      beta_save[, k_tmp] <- beta
      s2_epsilon_save[k_tmp] <- s2_epsilon
      gamma2_save[, k_tmp] <- 1 / gamma2_inv
      lambda2_save[k_tmp] <- lambda2
    }
  }
  
  ##
  ## output
  ##
  
  list(mu = mu_save, beta = beta_save, s2_epsilon = s2_epsilon_save, gamma2 = gamma2_save, lambda2 = lambda2_save, X = X)
}