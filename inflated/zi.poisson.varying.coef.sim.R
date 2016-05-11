###
### Simulate grouped zero-inflated count data and 
### fit model using zi.poisson.varying.coef.mcmc.R
###

rm(list=ls())
library(mvtnorm)

I <- 1000  # number of observations per group
J <- 10  # number of groups
g <- rep(1:J,each=I)  # grouping variable


##########################################################
### Population-level process model parameters
##########################################################

mu.beta <- matrix(c(-0.5,1),,1)  # mean of betas
qX <- nrow(mu.beta)
rho <- -0.15  # correlation between betas
Sigma <- diag(qX)*0.75  # variance-covariance of betas
Sigma[1,2] <- Sigma[2,1] <- Sigma[1,1]*Sigma[2,2]*rho


##########################################################
### Simulate group-level process model parameters
##########################################################

beta <- t(rmvnorm(J,mu.beta,Sigma))  # betas for each group
plot(t(beta))

##########################################################
### Simulate count data
##########################################################

p <- 0.5  # probability of Poisson mixture component
z <- rbinom(I*J,1,p)  # latent variable indicating mixture component
table(z)

y <- numeric(I*J)
X <- cbind(1,rnorm(I*J,0,sd=1))  # design matrix
beta.tmp <- t(beta[,g[z==1]])
lambda <- exp(rowSums(X[z==1,]*beta.tmp))  # intensity of Poisson process
hist(lambda,breaks=100)
y[z==1] <- rpois(sum(z==1),lambda)  # observed counts

##########################################################
### Fit model
##########################################################

source('~/Documents/git/Multilevel/nested/poisson/zi.poisson.varying.coef.mcmc.R')
start <- list(beta=beta,mu.beta=mu.beta,Sigma=Sigma,p=rep(p,J))
priors <- list(sigma.beta=5,S0=diag(qX),nu=qX+1,a=1,b=1)
tune <- list(beta=rep(50,J))
# tune <- list(beta=out1$tune$beta)
out1 <- zi.poisson.varying.coef.mcmc(y,X,g,priors,start,tune,adapt=TRUE,1000)
out1$tune

# Examine estimates for beta_j
g.idx <- 5  # group idx for plotting beta_j
matplot(out1$beta[,,g.idx],type="l",lty=1);abline(h=beta[,g.idx],col=1:qX,lty=2)

# Examine estimates for mu.beta
matplot(out1$mu.beta,type="l");abline(h=mu.beta,col=1:qX,lty=2)

# Examine estimate for p
matplot(out1$p,type="l");abline(h=p,col=1:J,lty=2)

# Examine estimates for z
boxplot(z,apply(out1$z,2,mean))

library(lattice)
bwplot(apply(out1$z,2,mean)~z,horizontal=FALSE)
bwplot(apply(out1$z,2,mean)~z,horizontal=FALSE,
	panel=panel.violin)

# Examine estimates for Lambda
matplot(cbind(out1$Sigma[1,1,],out1$Sigma[1,2,]),type="l")
abline(h=c(Sigma[1,1],Sigma[1,2]),lty=2,col=1:qX)



