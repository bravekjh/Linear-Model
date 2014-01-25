##
## libraries and functions
##

source('dinvgamma.R')
source('mcmc.lm.R')
source('rMVN.R')

make.model.plot <- function(out){
  layout(matrix(1:4, 2))
  hist(out$beta.save[1,])
  abline(v = beta[1], col = 'red')
  hist(out$beta.save[2,])
  abline(v = beta[2], col = 'red')
  plot(out$sigma.squared.beta.save, type = 'l')
  plot(out$sigma.squared.epsilon.save, type = 'l')
  abline(h = sigma.squared.epsilon, col = 'red')
}


##
## Simulate some data
##

N <- 1000
n <- 100

X <- cbind(rep(1, N), seq(0, 1, length = N))
beta <- c(0, 2)
sigma.squared.epsilon <- 0.25

if(is.null(dim(X))){
  Y <- X * beta + rnorm(N, 0, sigma.squared.epsilon)
} else {
  Y <- X %*% beta + rnorm(N, 0, sigma.squared.epsilon)
}

plot(Y ~ X[, 2])
lm(Y~X[, 2])

samp <- sample(1:N, n)
Y.samp <- Y[samp]
X.samp <- X[samp, ]

plot(Y.samp ~ X.samp[, 2])
lm(Y.samp ~ X.samp[, 2])

##
## Setup priors
##

# hyperparameters for mu.beta and sigma.squared.beta
mu.0 <- c(0, 0)
sigma.squared.0 <- 100 
# hyerparamters for sigma.squared.beta
alpha.beta <- 2
beta.beta <- 10
curve(dinvgamma(x, alpha.beta, beta.beta), from = 0, to = 10)
# hyperparameters for sigma.squared.epsilon
alpha.epsilon <- 2
beta.epsilon <- 10
curve(dinvgamma(x, alpha.epsilon, beta.epsilon), from = 0, to = 10)
n.mcmc <- 5000

##
## Fit mcmc
##

out <- mcmc.lm(Y.samp, X.samp, mu.0, sigma.squared.beta.0, alpha.beta, beta.beta, alpha.epsilon, beta.epsilon)

make.model.plot(out)
