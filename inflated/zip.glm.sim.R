
rm(list=ls())

##########################################################
### Simulate data for zero-inflated Poisson GLM
##########################################################

n <- 1000  # number of obserations to simulate
X <- cbind(1,rnorm(n))  # design matrix
qX <- ncol(X)
beta <- c(-2,3)  # coefficients

p <- 0.25  # probability of Poisson mixture component
z <- rbinom(n,1,p)  # latent variable indicating mixture component
table(z)

y <- numeric(n)
lambda <- exp(X[z==1,]%*%beta)  # intensity of Poisson process
y[z==1] <- rpois(sum(z==1),lambda)  # observed counts


##########################################################
### Fit model
##########################################################

source('~/Documents/git/GLM/poisson/zip.glm.mcmc.R')
priors <- list(sigma.beta=5,a=1,b=1)
tune <- list(beta=0.05)
start <- list(beta=coef(glm(z ~ 0+X, family=poisson())),p=p)
start <- list(beta=beta,p=p)
out1 <- zip.glm.mcmc(y,X,priors,start,tune,n.mcmc=5000)

matplot(out1$beta, type="l", lty=1);abline(h=beta,col=1:qX,lty=3)
apply(out1$beta,2,mean)
plot(out1$p,type="l");abline(h=p,lty=3)
boxplot(z,apply(out1$z,2,mean))

library(lattice)
bwplot(apply(out1$z,2,mean)~z,horizontal=FALSE)
bwplot(apply(out1$z,2,mean)~z,horizontal=FALSE,
	panel=panel.violin)