###
### Simulate zero-inflated Poisson distributed data collected in a three-level 
### hierarchy and fit model using zi.poisson.varying.coef.3.mcmc.R
###

rm(list=ls())

library(mvtnorm)

I <- 12  # number of groups
J <- 5  # number subgroups within each group
K <- 1000  # number of observations per subgroup

g1 <- rep(1:I,each=J*K)  # level 1 grouping variable (groups, indexing I)
g2 <- rep(rep(1:J,each=K),I)  # level 2 grouping variable (subgroups, indexing J)
cbind(g1,g2)


##########################################################
### Population-level process model parameters
##########################################################

mu.beta <- matrix(c(-0.5,1),,1)  # mean of betas
qX <- nrow(mu.beta)
rho <- -0.15  # correlation between betas
Sigma.beta <- diag(qX)*0.75  # variance-covariance of betas
Sigma.beta[1,2] <- Sigma.beta[2,1] <- Sigma.beta[1,1]*Sigma.beta[2,2]*rho

##########################################################
### Simulate group-level process model parameters
##########################################################

beta <- t(rmvnorm(I,mu.beta,Sigma.beta))  # betas for each group
plot(t(beta))

##########################################################
### Simulate subgroup-level process model parameters
##########################################################

Sigma.alpha <- diag(qX)*0.25  # variance-covariance of betas
Sigma.alpha[1,2] <- Sigma.alpha[2,1] <- Sigma.alpha[1,1]*Sigma.alpha[2,2]*rho
Sigma.alpha <- lapply(1:I,function(x) Sigma.alpha)
alpha <- lapply(1:I,function(x) t(rmvnorm(J,beta[,x],Sigma.alpha[[x]])))
	# coefficients for each subgroup

##########################################################
### Simulate count data
##########################################################

p <- 0.35  # probability of Poisson mixture component
z <- rbinom(I*J*K,1,p)  # latent variable indicating mixture component
table(z)

y <- numeric(I*J*K)
X <- cbind(1,rnorm(I*J*K,0,sd=1))  # design matrix

for(i in 1:length(y)){  # simulate counts
	lambda <- exp(t(X[i,])%*%alpha[[g1[i]]][,g2[i]])  # intensity of Poisson process
	y[i] <- z[i]*rpois(1,lambda)  # observed counts
}
hist(y,breaks=100)

##########################################################
### Fit model
##########################################################

source('~/Documents/git/Multilevel/nested/poisson/zi.poisson.varying.coef.3.mcmc.R')
start <- list(beta=beta,mu.beta=mu.beta,alpha=alpha,p=lapply(1:I,function(x) rep(p,J)),
	Sigma.beta=Sigma.beta,Sigma.alpha=Sigma.alpha)
priors <- list(sigma.mu.beta=5,S0=diag(qX),nu=qX+1,a=1,b=1)
tune <- list(alpha=lapply(1:I,function(x) rep(0.01,J)))
# tune <- list(beta=out1$tune$beta)
out1 <- zi.poisson.varying.coef.3.mcmc(y,X,g1,g2,priors,start,tune,adapt=TRUE,1000)

# Examine estimates for beta_i
g1.idx <- 5  # group idx for plotting beta_j
matplot(out1$beta[,,g1.idx],type="l",lty=1);abline(h=beta[,g1.idx],col=1:qX,lty=2)

# Examine estimates for alpha_ij
g2.idx <- 5  # group idx for plotting beta_j
matplot(out1$alpha[[g1.idx]][,,g2.idx],type="l",lty=1)
abline(h=alpha[[g1.idx]][,g2.idx],col=1:qX,lty=2)

# Examine estimates for p_ij
matplot(out1$p[[g1.idx]],type="l")
abline(h=p)

# Examine estimates for Sigma.alpha_ij
matplot(cbind(out1$Sigma.alpha[[g1.idx]][1,1,],out1$Sigma.alpha[[g1.idx]][1,2,]),type="l")
abline(h=c(Sigma.alpha[[g1.idx]][1,1],Sigma.alpha[[g1.idx]][1,2]),lty=2,col=1:qX)

# Examine estimates for mu.beta
matplot(out1$mu.beta,type="l");abline(h=mu.beta,col=1:qX,lty=2)

# Examine estimates for Sigma.beta
matplot(cbind(out1$Sigma.beta[1,1,],out1$Sigma.beta[1,2,]),type="l")
abline(h=c(Sigma.beta[1,1],Sigma.beta[1,2]),lty=2,col=1:qX)
