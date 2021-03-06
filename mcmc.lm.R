##
## Simple linear regression model 
##
## John Tipton - created 01.25.2014
##

##
## model: Y = X %*% beta + epsilon
##

##
## libraries and functions
##

mcmc.lm <- function(Y, X, n.mcmc, mu.0, sigma.squared.beta.0, alpha.beta, beta.beta, alpha.epsilon, beta.epsilon){

##
## Initialize variables
##

n <- length(Y)
tau <- dim(X)[2]
n.burn <- floor(n.mcmc / 5) + 1
#
I.beta <- diag(tau)
sigma.squared.beta <- 1 / rgamma(1, alpha.beta, beta.beta)
Sigma.beta <- sigma.squared.beta * I.beta
Sigma.beta.inv <- 1 / sigma.squared.beta * I.beta 
#
I.epsilon <- diag(n)
sigma.squared.epsilon <- 1 / rgamma(1, alpha.epsilon, beta.epsilon) 
Sigma.epsilon <- sigma.squared.epsilon * I.epsilon
Sigma.epsilon.inv <- 1 / sigma.squared.epsilon * I.epsilon
#
Sigma.0 <- sigma.squared.0 * I.beta
Sigma.0.inv <- 1 / sigma.squared.0 * I.beta

mu.beta <- rMVN(A.chol = chol(Sigma.0), b = Sigma.0 %*% mu.0)
beta <- rMVN(A.chol = chol(Sigma.beta), b = Sigma.beta %*% mu.beta) 

##
## save variables
##

beta.save <- matrix(nrow = tau, ncol = n.mcmc)
sigma.squared.beta.save <- vector(length = n.mcmc)
sigma.squared.epsilon.save <- vector(length = n.mcmc)
Dbar.save <- vector(length = n.mcmc)
Y.pred.mn <- vector(length = n)

##
## Start MCMC
##

for(k in 1:n.mcmc){
  if(k %% 100 == 0){
    cat(" ", k)
  }

  ##
  ## sample beta
  ##
  
  A.chol <- chol(t(X) %*% Sigma.epsilon.inv %*% X + Sigma.beta.inv)
  b <- (t(X) %*% Sigma.epsilon.inv %*% Y + Sigma.beta.inv %*% mu.beta)
  beta <- rMVN(A.chol, b)
  
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
  ## DIC
  ##
  
  Dbar.save[k] <- -2 * sum(dnorm(Y, X %*% beta, sqrt(sigma.squared.epsilon), log = TRUE))
  ##
  ## Posterior Predictive Calculations 
  ##
  
  if(k > n.burn){
    Y.pred <- rMVN(chol(Sigma.epsilon.inv), Sigma.epsilon.inv %*% X %*% beta)
    Y.pred.mn <- Y.pred.mn + Y.pred / (n.mcmc - n.burn)
  }
  
  ##
  ## save variables
  ##
  
  beta.save[, k] <- beta
  sigma.squared.epsilon.save[k] <- sigma.squared.epsilon
  sigma.squared.beta.save[k] <- sigma.squared.beta
  }
###
###  Calculate DIC and Print to Screen
###

if(dim(X)[2] == 1){
  postbetamn <- mean(beta.save[, - (1:n.burn)])
}
if(dim(X)[2] > 1){
  postbetamn <- apply(beta.save[, -(1:n.burn)], 1, mean)
}
posts2mn <- mean(sigma.squared.epsilon.save[ - (1:n.burn)])
cat("\n", "Posterior Mean for Beta:", "\n")
print(postbetamn)
cat("Posterior Mean for s2:", "\n")
print(posts2mn)
Dhat <- -2 * (sum(dnorm(Y, X %*% postbetamn, sqrt(posts2mn), log=TRUE)))
Dbar <- mean(Dbar.save[ - (1:n.burn)])
pD <- Dbar - Dhat
DIC <- Dhat + 2 * pD

cat("Dhat:", Dhat, "Dbar:", Dbar, "pD:", pD, "DIC:", DIC, "\n")

##
## output
##

list(beta.save = beta.save, sigma.squared.beta.save = sigma.squared.beta.save, sigma.squared.epsilon.save = sigma.squared.epsilon.save, Y.pred.mn = Y.pred.mn)
}