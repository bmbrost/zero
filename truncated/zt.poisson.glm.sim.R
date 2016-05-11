
rm(list=ls())

##########################################################
### Simulate data for zero-truncated Poisson GLM
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

source('~/Documents/git/GLM/poisson/zt.poisson.glm.mcmc.R')
priors <- list(sigma.beta=5)
tune <- list(beta=0.01)
start <- list(beta=beta)
out1 <- zt.poisson.glm.mcmc(z[z>0],X[z>0,],priors,start,tune,adapt=TRUE,n.mcmc=5000)

matplot(out1$beta, type="l", lty=1);abline(h=beta,col=1:qX,lty=3)


