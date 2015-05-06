rm(list = ls())
##
## libraries and functions
##

library(statmod)
setwd('~/Linear-Model/')
source('dinvgamma.R')
source('rMVN.R')

setwd("~/Linear-Model/bayesianLassoRegression/fixedEffectModel/")
source("mcmc.lm.lasso.R")

make.model.plot <- function(out){
  layout(matrix(1:6, 3, 2))
  matplot(out$mu, type = 'l')
  matplot(t(out$beta), type = 'l')
  hist(out$beta[2,], main = 'Posterior of Beta2')
  abline(v = beta[2], col = 'red')
  plot(out$s2_epsilon, type = 'l')
  abline(h = s2_epsilon, col = 'red')
  plot(out$lambda2, type = 'l')
}


##
## Simulate some data
##

N <- 1000
n <- 100
beta <- -3:3
s2_epsilon <- 0.25
tau <- length(beta)

make.lm.data <- function(N, n, beta, s.sqaured_epsilon){
  tau <- length(beta)
  X <- matrix(nrow = N, ncol = tau)
  for(i in 1:tau){
    X[, i] <- rnorm(N, 0, 1)
  }
#   X <-matrix(c(rep(1, N), rep(seq(0, 1, length = N), tau - 1)), nrow = N, ncol = tau)
#   if(is.null(dim(X))){
#     Y <- X * beta + rnorm(N, 0, s2_epsilon)
#   } else {
    Y <- X %*% beta + rnorm(N, 0, s2_epsilon)
#   }
  #list(X = X, Y = Y, N = N, n = n, s2_epsilon = s2_epsilon)
  data.frame(Y, X)
}

data <- make.lm.data(N, n, beta, s2_epsilon)

samp <- sample(1:N, n)
data.samp <- data[samp, ]

lm(Y ~ . ,data = data)

##
## Setup priors
##

# hyperparameters for mu.beta and s2.beta
alpha_epsilon <- 0.1
beta_epsilon <- 0.1
alpha_lambda <- 0.1
beta_lambda <- 0.1
# 
# mu.0 <- rep(0, tau)
# s2.0 <- 100 
# # hyerparamters for s2.beta
# alpha.beta <- 2
# beta.beta <- 10
# curve(dinvgamma(x, alpha.beta, beta.beta), from = 0, to = 10)
# # hyperparameters for s2_epsilon
# alpha_epsilon <- 2
# beta_epsilon <- 10
# curve(dinvgamma(x, alpha_epsilon, beta_epsilon), from = 0, to = 10)

n_mcmc <- 5000

##
## Fit mcmc
##

Y <- data.samp[, 1]
X <- as.matrix(data.samp[, 2:(tau + 1)], ncol = tau)

out <- mcmc(Y, X, n_mcmc, alpha_epsilon, beta_epsilon, alpha_lambda, beta_lambda)
              

make.model.plot(out)
beta
rowMeans(out$beta)

