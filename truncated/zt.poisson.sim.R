###
### Simulate zero-truncated Poisson count data and 
### fit model using zt.poisson.mcmc.R
###

rm(list=ls())

##########################################################
### Simulate zero-truncated count data
##########################################################

n <- 1000  # number of obserations to simulate
X <- cbind(1,rnorm(n))  # design matrix
qX <- ncol(X)
beta <- c(1,2)  # coefficients
lambda <- exp(X%*%beta)  # intensity of Poisson process
z <- rpois(n,lambda)  # observed counts

sum(z>0)

##########################################################
### Fit model
##########################################################

source('~/Documents/git/zero/truncated/zt.poisson.mcmc.R')
priors <- list(sigma.beta=5)
tune <- list(beta=0.01)
start <- list(beta=beta)
out1 <- zt.poisson.mcmc(z[z>0],X[z>0,],priors,start,tune,adapt=TRUE,n.mcmc=5000)

matplot(out1$beta, type="l", lty=1);abline(h=beta,col=1:qX,lty=3)


