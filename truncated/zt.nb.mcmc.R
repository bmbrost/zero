#
#
# Zero-truncated negative binomial model for count data
#
# Function name: zt.nb.MCMC
#
# Author: Brian M. Brost
# Contact: bmbrost@gmail.com
#
# Last updated: 23 MAR 2016
#
# Model statement:
#	z_i|z_i>0 ~ ZTNB(lambda_i,alpha)
#	log(lambda_i) = x_i%*%beta
#	beta ~ N(0,sigma.beta^2*I)
# 	alpha ~ gamma(a,b)	
#	note: E[alpha]=a/b and Var[alpha]=a/(b^2)
#
# Reference:
#
# Required R packages: 
#
# Inputs:
#
# z - Vector of length n containing the count of observations corresponding to 
# 	each row in the design matrix X. Note that the value of z[1] corresponds to
#	X[1,], z[2] corresponds to X[2,], etc. Also note that z>0.
# X - Design matrix of dimension n x qX containing covariates (plus
#	intercept) for which inference is desired
# priors - List of priors containing the following elements:
#	1. sigma.beta - Standard deviation of normal prior on beta
#	2. a - Shape parameter of gamma prior for alpha 
#	3. b - Rate parameter of gamma prior for alpha
# start - List of starting values containing the following elements:
#	1. beta - Vector of starting values for coefficients
#	2. alpha - Over-dispersion parameter for negative binomial distribution 
# tune - list of tuning parameters containing the following elements:
#	1. beta - Tuning parameter for Metropolis-Hastings update on beta
#	2. alpha - Tuning parameter for Metropolis-Hastings update on alpha
# adapt - Switch to enable adapative tuning (TRUE/FALSE)
# n.mcmc - Number of desired MCMC iterations
#
#

zt.nb.mcmc <- function(z,X,priors,start,tune,adapt=TRUE,n.mcmc=1000){

	###
	###  Libraries and Subroutines
	###

	get.tune <- function(tune,keep,k,target=0.44){  # adaptive tuning
		# a <- min(0.01,1/sqrt(k))
		a <- min(0.025,1/sqrt(k))
		exp(ifelse(keep<target,log(tune)-a,log(tune)+a))
	}
	
	
	###
	###  Setup variable, starting values, and priors
	###

	qX <- ncol(X)  # number of betas
	
	beta <- start$beta
	lambda <- exp(X%*%beta)  # mean of negative binomial
	alpha <- start$alpha  # over-dispersion parameter in negative binomial

	mu.beta <- rep(0,qX)  # mean of normal prior on beta
	sigma.beta <- priors$sigma.beta
	a <- priors$a  # shape parameter of gamma prior for alpha
	b <- priors$b  # rate parameter of gamma prior for alpha
	
	# Sigma.beta <- solve(crossprod(X))
	# c <- sigma.beta/min(diag(Sigma.beta))
	# Sigma.beta <- c*Sigma.beta  # variance-covariance matrix of
		# # normal prior on beta

	###
	### Create receptacles for output
	###
  
	beta.save <- matrix(0,n.mcmc,qX)  # rsf coefficients
	alpha.save <- numeric(n.mcmc)  # over-dispersion parameter
	
	keep <- list(beta=0,alpha=0)
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
			tune$alpha <- get.tune(tune$alpha,keep.tmp$alpha,k)
			keep.tmp <- lapply(keep.tmp,function(x) x*0)
	   	} 	

		###
		### Update beta
		### 
		
		# See Zuur et al. 2007. Mixed effects models and 
		# extensions in ecology with R 
		beta.star <- rnorm(qX,beta,tune$beta)
		lambda.star <- exp(X%*%beta.star)  # mean of negative binomial	
		mh.0.beta <- sum(dnbinom(z,size=alpha,mu=lambda,log=TRUE))-
 			sum(log(1-(alpha/(lambda+alpha))^alpha))+
			sum(dnorm(beta,mu.beta,sigma.beta,log=TRUE))
		mh.star.beta <- sum(dnbinom(z,size=alpha,mu=lambda.star,log=TRUE))-
			sum(log(1-(alpha/(lambda.star+alpha))^alpha))+
			sum(dnorm(beta.star,mu.beta,sigma.beta,log=TRUE))
		if(exp(mh.star.beta-mh.0.beta)>runif(1)){
			beta <- beta.star
			lambda <- lambda.star
			keep$beta <- keep$beta+1
			keep.tmp$beta <- keep.tmp$beta+1
		}

		###
		### Update alpha (over-dispersion parameter)
		### 

		alpha.star <- rnorm(1,alpha,tune$alpha)
		if(alpha.star>0){
	  		mh.0.alpha <- sum(dnbinom(z,size=alpha,mu=lambda,log=TRUE))+
				dgamma(alpha,shape=a,rate=b,log=TRUE)
			mh.star.alpha <- sum(dnbinom(z,size=alpha.star,mu=lambda,log=TRUE))+
				dgamma(alpha.star,shape=a,rate=b,log=TRUE)
			if(exp(mh.star.alpha-mh.0.alpha)>runif(1)){
				alpha <- alpha.star
				keep$alpha <- keep$alpha+1
				keep.tmp$alpha <- keep.tmp$alpha+1
			}		
		}

		###
		###  Save samples 
	    ###

		beta.save[k,] <- beta
		alpha.save[k] <- alpha
	}
	
	keep$beta <- keep$beta/n.mcmc
	keep$alpha <- keep$alpha/n.mcmc
	cat(paste("\nbeta acceptance rate:",round(keep$beta,2))) 
	cat(paste("\nalpha acceptance rate:",round(keep$alpha,2))) 

	###
	### Write output
	###

	list(beta=beta.save,alpha=alpha.save,keep=keep,
		z=z,X=X,priors=priors,start=start,tune=tune,n.mcmc=n.mcmc)
}