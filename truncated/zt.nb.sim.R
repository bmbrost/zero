###
### Simulate zero-truncated negative binomial count data and 
### fit model using zt.nb.mcmc.R
###

rm(list=ls())

##########################################################
### Simulate zero-truncated count data
##########################################################

n <- 10000  # number of obserations to simulate
X <- cbind(1,rnorm(n))  # design matrix
qX <- ncol(X)
beta <- c(-1,3)  # coefficients
mu <- exp(X%*%beta)  # intensity of Poisson process
hist(mu,breaks=100)
alpha <- 2 # dispersion parameter; var(z)=mu+mu^2/alpha, i.e.,  alpha=\infty yields Poisson
z <- rnbinom(n,size=alpha,mu=mu)  # observed counts
hist(z,breaks=20)

sum(z>0)

##########################################################
### Fit model
##########################################################

source('~/Documents/git/zero/truncated/zt.nb.mcmc.R')
priors <- list(sigma.beta=10,a=1,b=0.1)
# hist(rgamma(1000,1,0.01),breaks=100)  # vague parameterization
tune <- list(beta=0.01,alpha=0.1)
start <- list(beta=beta,alpha=alpha)
out1 <- zt.nb.mcmc(z[z>0],X[z>0,],priors,start,tune,adapt=TRUE,n.mcmc=1000)

matplot(out1$beta, type="l", lty=1);abline(h=beta,col=1:qX,lty=3)
apply(out1$beta,2,mean)
matplot(out1$alpha, type="l", lty=1);abline(h=alpha,col=1,lty=3)
mean(out1$alpha)


# Compare to canned funcation for zero-truncated NB regression
install.packages("countreg", repos="http://R-Forge.R-project.org",type="source")
library(countreg)
zerotrunc(z[z>0]~X[z>0,2],dist="negbin")

install.packages("VGAM")
library("VGAM")
ztnb <- vglm(z[z>0]~X[z>0,2],family=posnegbinomial)
coef(ztnb)[-2]  # intercept and slopes
ztnb@dispersion
summary(ztnb)