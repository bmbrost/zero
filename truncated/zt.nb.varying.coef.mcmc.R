#
#
# Bayesian zero-truncated negative binomial generalized linear mixed model for count data
#
# Function name: zt.nb.varying.coef.MCMC
#
# Author: Brian M. Brost
# Contact: bmbrost@gmail.com
#
# Last updated: 23 MAR 2016
#
# Model statement:
#	z_ij|z_ij>0 ~ ZTNB(lambda_ij,alpha)
#	log(lambda_ij) = x_ij%*%beta_j
# 	beta_j ~ N(mu_beta,Sigma)
#	mu_beta ~ N(0,sigma.beta^2*I)
#	Sigma ~ Wish(S_0,nu)
#	alpha ~ Gamma(a,b)
#	note: E[z_ij]=lambda_ij and Var[z_ij]=lambda_ij+lambda_ij^2/alpha
#	note: E[alpha]=a/b and Var[alpha]=a/(b^2)
#
# Reference:
#
# Required R packages: mvtnorm (if using multivariate normal prior on beta)
#
# Inputs:
#
# z - Vector of length n containing the count of observations corresponding to 
# 	each row in the design matrix X. Note that the value of z[1] corresponds to
#	X[1,], z[2] corresponds to X[2,], etc. Also note that z>0.
# X - Design matrix of dimension n x qX containing covariates (plus
#	intercept) for which inference is desired
# g - variable that defines groups of observations in z
# priors - list of priors containing the following elements:
#	1. sigma.beta - Standard deviation of normal prior on mu.beta
#	2. S0 - Scale matrix for the inverse-Wishart prior on Sigma
#	3. nu - Degrees of freedom for the IW prior on Sigma
#	4. a - Shape parameter of gamma prior for alpha 
#	5. b - Rate parameter of gamma prior for alpha
# start - list of starting values containing the following elements:
#	1. beta - Vector of starting values for coefficients
#	2. mu.beta - Vector of starting values for mean of betas
#	3. Sigma - Variance-covariance matrix for betas
#	4. alpha - Over-dispersion parameter for negative binomial distribution
# tune - List of tuning parameters containing the following elements:
#	1. beta - Tuning parameter for Metropolis-Hastings update on beta
#	2. alpha - Tuning parameter for Metropolis-Hastings update on alpha
# adapt - Switch to enable adapative tuning (TRUE/FALSE)
# n.mcmc - Number of desired MCMC iterations
#
#

zt.nb.varying.coef.mcmc <- function(z,X,g,priors,start,tune,adapt=TRUE,n.mcmc=1000){

	###
	###  Libraries and Subroutines
	###

	library(mvtnorm)

	get.tune <- function(tune,keep,k,target=0.44){  # adaptive tuning
		# a <- min(0.01,1/sqrt(k))
		a <- min(0.025,1/sqrt(k))
		exp(ifelse(keep<target,log(tune)-a,log(tune)+a))
	}

	zt.dpois <- function(z,X,beta){  # log-likelihood for zero-truncated Poisson
		# See Eq. 11.1 in Hilbe. 2011. Negative binomial regression. 
		Xb <- X%*%beta
		z*Xb-exp(Xb)-lfactorial(z)-log(1-exp(-exp(Xb)))	
	}
	
	
	###
	###  Setup variable, starting values, and priors
	###

	J <- length(unique(g))  # number of groups
	g <- as.numeric(g)
	g.idx <- sapply(sort(unique(g)),function(x) which(g==x),simplify=FALSE)
		# indexes of observations in y by group
	n.j <- unlist(lapply(g.idx,length))  # number of observations per group
	n <- length(z)  # total number of observations
	qX <- ncol(X)

	###
	###  Starting values, and priors
	###

	beta <- matrix(start$beta,qX,J)
	lambda <- sapply(1:J,function(x) exp(X[g.idx[[x]],]%*%beta[,x]))  # mean of negative binomial

	mu.beta <- matrix(start$mu.beta,qX,1)  # mean of normal prior on beta
	Sigma <- start$Sigma
	Sigma.inv <- solve(Sigma)
	
	Sigma.beta <- priors$sigma.beta^2*diag(qX)
	Sigma.beta.inv <- solve(Sigma.beta)
	S0 <- priors$S0
	nu <- priors$nu
	
	# Sigma.beta <- solve(crossprod(X))
	# c <- sigma.beta/min(diag(Sigma.beta))
	# Sigma.beta <- c*Sigma.beta  # variance-covariance matrix of
		# # normal prior on beta

	###
	### Create receptacles for output
	###
  
	beta.save <- array(0,dim=c(n.mcmc,qX,J))
	mu.beta.save <- matrix(0,n.mcmc,qX)
	Sigma.save <- array(0,dim=c(qX,qX,n.mcmc))
	alpha.save <- numeric(n.mcmc)  # over-dispersion parameter

	keep <- list(beta=rep(0,J),alpha=0)
	keep.tmp <- keep  # track MH accpetance rate for adaptive tuning
	Tb <- 50  # frequency of adaptive tuning
	
	###
	###  Begin MCMC loop
	###

	for(k in 1:n.mcmc){

		if(k%%1000==0) cat(k,"");flush.console()	

		if(adapt==TRUE & k%%Tb==0) {  # Adaptive tuning
			keep.tmp <- lapply(keep.tmp,function(x) x/Tb)
			# tune$beta <- sapply(1:J,function(x) get.tune(tune$beta[x],keep.tmp$beta[x],k))
			tune$beta <- get.tune(tune$beta,keep.tmp$beta,k)
			tune$alpha <- get.tune(tune$alpha,keep.tmp$alpha,k)
			keep.tmp <- lapply(keep.tmp,function(x) x*0)
	   	} 	

		###
		### Update beta_j
		### 
		
		for(i in 1:J){
			idx <- g.idx[[i]]
			# beta.star <- rnorm(qX,beta[,i],tune$beta)
			tune.tmp <- (tune$beta[i]/n.j[i])*solve(crossprod(X[idx,]))
			beta.star <- c(rmvnorm(1,beta[,i],tune.tmp))
			lambda.star <- exp(X[idx,]%*%beta.star)  # mean of negative binomial	
	  		mh.0 <- sum(dnbinom(z[idx],size=alpha,mu=lambda[[i]],log=TRUE))-
				sum(log(1-(alpha/(lambda[[i]]+alpha))^alpha))+	  		 
				sum(dmvnorm(beta[,i],mu.beta,Sigma,log=TRUE))
			mh.star <- sum(dnbinom(z[idx],size=alpha,mu=lambda.star,log=TRUE))-
				sum(log(1-(alpha/(lambda.star+alpha))^alpha))+
				sum(dmvnorm(beta.star,mu.beta,Sigma,log=TRUE))
			if(exp(mh.star-mh.0)>runif(1)){
				beta[,i] <- beta.star
				lambda[[i]] <- lambda.star
				keep$beta[i] <- keep$beta[i]+1
				keep.tmp$beta[i] <- keep.tmp$beta[i]+1
			}
		}
		
	  	###
	  	### Sample mu_beta
	  	###

		beta.mean <- apply(beta,1,sum)
		A.inv <- solve(J*Sigma.inv+Sigma.beta.inv)
		b <- Sigma.inv%*%beta.mean
	    mu.beta <- A.inv%*%b+t(chol(A.inv))%*%matrix(rnorm(qX),qX,1)

	  	###
	  	### Sample Sigma
	  	###
		# browser()		
	  	Sn <- S0+crossprod(t(beta)-matrix(mu.beta,J,qX,byrow=TRUE))
		Sigma <- solve(rWishart(1,nu+J,solve(Sn))[,,1])
		Sigma.inv <- solve(Sigma)

		###
		### Update alpha (over-dispersion parameter)
		### 

# browser()
		alpha.star <- rnorm(1,alpha,tune$alpha)
		if(alpha.star>0){
			lambda.tmp <- unlist(lambda)
	  		mh.0.alpha <- sum(dnbinom(z,size=alpha,mu=lambda.tmp,log=TRUE))+
				dgamma(alpha,shape=priors$a,rate=priors$b,log=TRUE)
			mh.star.alpha <- sum(dnbinom(z,size=alpha.star,mu=lambda.tmp,log=TRUE))+
				dgamma(alpha.star,shape=priors$a,rate=priors$b,log=TRUE)
			if(exp(mh.star.alpha-mh.0.alpha)>runif(1)){
				alpha <- alpha.star
				keep$alpha <- keep$alpha+1
				keep.tmp$alpha <- keep.tmp$alpha+1
			}		
		}

		###
		###  Save samples 
	    ###

	  	beta.save[k,,] <- beta
		mu.beta.save[k,] <- mu.beta
		Sigma.save[,,k] <- Sigma
		alpha.save[k] <- alpha

	}
# browser()
	cat("\n")	
	keep$beta <- keep$beta/n.mcmc
	keep$alpha <- keep$alpha/n.mcmc
	cat("beta acceptance rate:",round(keep$beta,2)) 
	cat(paste("\nalpha acceptance rate:",round(keep$alpha,2))) 

	###
	### Write output
	###

	list(beta=beta.save,mu.beta=mu.beta.save,Sigma=Sigma.save,alpha=alpha.save,keep=keep,
		z=z,X=X,g,priors=priors,start=start,tune=tune,n.mcmc=n.mcmc)
}