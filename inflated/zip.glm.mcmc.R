#
#
# Bayesian zero-inflated Poisson generalized linear model
#
# Function name: zip.glm.MCMC
#
# Author: Brian M. Brost
# Contact: bmbrost@gmail.com
#
# Last updated: 05 April 2016
#
# Model statement:
#	y_i ~ Pois(lambda_i), z_i=1
#	y_i ~ 0, z_i=0
# 	z_i ~ Bern(p)
#	log(lambda_i) = x_i%*%beta
#	beta ~ N(0,sigma.beta^2*I)
#	p ~ Beta(a,b)
#
# Reference:
#
# Required R packages: 
#
# Inputs:
#
# y - vector of length n containing the count of observations corresponding to 
# 	each row in the design matrix X. Note that the value of z[1] corresponds to
#	X[1,], z[2] corresponds to X[2,], etc.
# X - design matrix of dimension n x qX containing covariates (plus
#	intercept) for which inference is desired
# priors - dist of priors containing the following elements:
#	1. sigma.beta - Standard deviation of normal prior on beta
#	2. a - first scale parameter of beta prior on p
#	3. b - second scale parameter of beta prior on p
# start - list of starting values containing the following elements:
#	1. beta - vector of starting values for coefficients
#	2. p - probability associated with z, the latent mixture component indicator variable
# tune - list of tuning parameters containing the following elements:
#	1. beta - tuning parameter for Metropolis-Hasting update on beta
# adapt - switch to enable adapative tuning (TRUE/FALSE)
# n.mcmc - number of desired MCMC iterations
#
#

zip.glm.mcmc <- function(y,X,priors,start,tune,adapt=TRUE,n.mcmc=1000){

	###
	###  Libraries and Subroutines
	###

	# library(mvtnorm)

	get.tune <- function(tune,keep,k,target=0.44){  # adaptive tuning
		# a <- min(0.01,1/sqrt(k))
		a <- min(0.025,1/sqrt(k))
		exp(ifelse(keep<target,log(tune)-a,log(tune)+a))
	}
	
	###
	###  Setup variable, starting values, and priors
	###
# browser()
	qX <- ncol(X)  # number of betas
	y0 <- which(y==0)  # zero valued observations
	n.y0 <- length(y0)
	n <- length(y)
	z <- rep(1,n)  # latent mixture component indicator
	
	beta <- start$beta
	lambda <- exp(X%*%beta)  # intensity of Poisson process
	p <- start$p 
	
	mu.beta <- rep(0,qX)  # mean of normal prior on beta
	sigma.beta <- priors$sigma.beta
	
	# Sigma.beta <- solve(crossprod(X))
	# c <- sigma.beta/min(diag(Sigma.beta))
	# Sigma.beta <- c*Sigma.beta  # variance-covariance matrix of
		# # normal prior on beta

	###
	### Create receptacles for output
	###
  
	beta.save <- matrix(0,n.mcmc,qX)  # rsf coefficients
	p.save <- numeric(n.mcmc)
	z.save <- matrix(0,n.mcmc,n)
			
	keep <- list(beta=0)
	keep.tmp <- keep  # track MH accpetance rate for adaptive tuning
	Tb <- 50  # frequency of adaptive tuning

	###
	###  Begin MCMC loop
	###

	for(k in 1:n.mcmc){

		if(k%%1000==0) cat(k,"");flush.console()	

		if(adapt==TRUE & k%%Tb==0) {  # Adaptive tuning
			keep.tmp <- lapply(keep.tmp,function(x) x/Tb)
			tune$beta <- get.tune(tune$beta,keep.tmp$beta,k)
			keep.tmp <- lapply(keep.tmp,function(x) x*0)
	   	} 			
			
		###
		### Update z
		###
		
		p.tmp <- p*exp(-lambda[y0])
		p.tmp <- p.tmp/(p.tmp+1-p)  
		z[y0] <- rbinom(n.y0,1,p.tmp)

		###
		### Update p
		###

		sum.z <- sum(z)
		p <- rbeta(1,sum.z+1,n-sum.z+1)

		###
		### Update beta
		### 
		
		beta.star <- rnorm(qX,beta,tune$beta)
		lambda.star <- exp(X%*%beta.star)  # intensity of Poisson process
  		mh.0 <- sum(z*dpois(y,lambda,log=TRUE))+
			sum(dnorm(beta,mu.beta,sigma.beta,log=TRUE))
		mh.star <- sum(z*dpois(y,lambda.star,log=TRUE))+
			sum(dnorm(beta.star,mu.beta,sigma.beta,log=TRUE))
		if(exp(mh.star-mh.0)>runif(1)){
			beta <- beta.star
			lambda <- lambda.star
			keep$beta <- keep$beta+1
			keep.tmp$beta <- keep.tmp$beta+1
		}

		###
		###  Save samples 
	    ###

		beta.save[k,] <- beta
		p.save[k] <- p
		z.save[k,] <- z
	}
	
	keep$beta <- keep$beta/n.mcmc
	cat(paste("\nbeta acceptance rate:",round(keep$beta,2))) 

	###
	### Write output
	###

	list(beta=beta.save,p=p.save,z=z.save,keep=keep,
		y=y,X=X,priors=priors,start=start,tune=tune,n.mcmc=n.mcmc)
}