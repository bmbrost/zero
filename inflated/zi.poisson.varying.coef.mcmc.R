#
#
# Bayesian zero-inlfated Poisson generalized linear mixed model for count data
#
# Function name: zi.poisson.varying.coef.MCMC
#
# Author: Brian M. Brost
# Contact: bmbrost@gmail.com
#
# Last updated: 05 April 2016
#
# Model statement:
#	y_ij ~ Pois(lambda_ij), z_ij=1
#	y_ij ~ 0, z_ij=0
#	z_ij ~ Bern(p_j)
#	log(lambda_ij) = x_ij%*%beta_j
# 	beta_j ~ N(mu_beta,Sigma)
#	mu_beta ~ N(0,sigma.beta^2*I)
#	Sigma ~ Wish(S_0,nu)
#	p_j ~ Beta(a,b)
#
# Reference:
#
# Required R packages: mvtnorm
#
# Inputs:
#
# y - vector of length n containing the count of observations corresponding to 
# 	each row in the design matrix X. Note that the value of z[1] corresponds to
#	X[1,], z[2] corresponds to X[2,], etc.
# X - design matrix of dimension n x qX containing covariates (plus
#	intercept) for which inference is desired
# g - variable that defines groups of observations in z
# priors - list of priors containing the following elements:
#	1. sigma.beta - standard deviation of normal prior on mu.beta
#	2. S0 - scale matrix for the inverse-Wishart prior on Sigma
#	3. nu - degrees of freedom for the IW prior on Sigma
#	4. a - first scale parameter of beta prior on p
#	5. b - second scale parameter of beta prior on p
# start - list of starting values containing the following elements:
#	1. beta - vector of starting values for coefficients
#	2. mu.beta - vector of starting values for mean of betas
#	3. Sigma - variance-covariance matrix for betas
#	4. p - probability associated with z, the latent mixture component indicator variable
# tune - list of tuning parameters containing the following elements:
#	1. beta - tuning parameter for Metropolis-Hasting update on beta
# adapt - switch to enable adapative tuning (TRUE/FALSE)
# n.mcmc - number of desired MCMC iterations
#
#

zi.poisson.varying.coef.mcmc <- function(y,X,g,priors,start,tune,adapt=TRUE,n.mcmc=1000){

	###
	###  Libraries and Subroutines
	###

	library(mvtnorm)

	get.tune <- function(tune,keep,k,target=0.44){  # adaptive tuning
		# a <- min(0.01,1/sqrt(k))
		a <- min(0.025,1/sqrt(k))
		exp(ifelse(keep<target,log(tune)-a,log(tune)+a))
	}
	
	
	###
	###  Setup variable
	###
# browser()
	J <- length(unique(g))  # number of groups
	g <- as.numeric(g)
	g.idx <- sapply(sort(unique(g)),function(x) which(g==x),simplify=FALSE)
		# indexes of observations in y by group
	n.j <- unlist(lapply(g.idx,length))  # number of observations per group
	y0 <- sapply(sort(unique(g)),function(x) which(g==x&y==0),simplify=FALSE)
		# zero valued observations by group
	n.y0 <- unlist(lapply(y0,length))
	n <- length(z)  # total number of observations
	qX <- ncol(X)
	z <- rep(1,n)  # latent mixture component indicator
	
	###
	###  Starting values, and priors
	###

	beta <- matrix(start$beta,qX,J)
	lambda <- sapply(1:J,function(x) exp(X[g.idx[[x]],]%*%beta[,x]))  # intensity of Poisson process
	p <- start$p
	
	mu.beta <- matrix(start$mu.beta,qX,1)
	Sigma <- start$Sigma
	Sigma.inv <- solve(Sigma)
	
	Sigma.beta <- priors$sigma.beta^2*diag(qX)
	Sigma.beta.inv <- solve(Sigma.beta)
	S0 <- priors$S0
	nu <- priors$nu

	###
	### Create receptacles for output
	###
  
	beta.save <- array(0,dim=c(n.mcmc,qX,J))
	p.save <- matrix(0,n.mcmc,J)
	z.save <- matrix(0,n.mcmc,n)
	mu.beta.save <- matrix(0,n.mcmc,qX)
	Sigma.save <- array(0,dim=c(qX,qX,n.mcmc))
	
	keep <- list(beta=rep(0,J))
	keep.tmp <- keep  # track MH accpetance rate for adaptive tuning
	Tb <- 50  # frequency of adaptive tuning

	###
	###  Begin MCMC loop
	###

	for(k in 1:n.mcmc){
# browser()
		if(k%%1000==0) cat(k,"");flush.console()	

		if(adapt==TRUE & k%%Tb==0) {  # Adaptive tuning
			# browser()
			keep.tmp <- lapply(keep.tmp,function(x) x/Tb)
			tune$beta <- get.tune(tune$beta,keep.tmp$beta,k)
			keep.tmp <- lapply(keep.tmp,function(x) x*0)
	   	} 	

		###
		### Update group-level parameters
		### 
		
		for(i in 1:J){  # loop through groups

			###
			### Update z_ij
			###

			idx <- y0[[i]]		
			p.tmp <- p[i]*exp(-lambda[idx])
			p.tmp <- p.tmp/(p.tmp+1-p[i])  
			z[idx] <- rbinom(n.y0[i],1,p.tmp)

			###
			### Update p_j
			###

			idx1 <- g.idx[[i]]
			z.tmp <- z[idx1]	
			sum.z <- sum(z.tmp)
			p[i] <- rbeta(1,sum.z+priors$a,n.j[i]-sum.z+priors$b)
	
			###
			### Update beta_j
			###	

			idx2 <- which(z.tmp==1)
			y.tmp <- y[idx1[idx2]]
			# beta.star <- rnorm(qX,beta[,i],tune$beta)
			tune.tmp <- (tune$beta[i]/n.j[i])*solve(crossprod(X[idx,]))
			beta.star <- c(rmvnorm(1,beta[,i],tune.tmp))
			lambda.star <- exp(X[idx1,]%*%beta.star)  # intensity of Poisson process
	  		mh.0 <- sum(dpois(y.tmp,lambda[idx2,i],log=TRUE))+
				sum(dmvnorm(beta[,i],mu.beta,Sigma,log=TRUE))
			mh.star <- sum(dpois(y.tmp,lambda.star[idx2],log=TRUE))+
				sum(dmvnorm(beta.star,mu.beta,Sigma,log=TRUE))
			if(exp(mh.star-mh.0)>runif(1)){
				beta[,i] <- beta.star
				lambda[,i] <- lambda.star
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
		###  Save samples 
	    ###

	  	beta.save[k,,] <- beta
		p.save[k,] <- p
		z.save[k,] <- z
		mu.beta.save[k,] <- mu.beta
		Sigma.save[,,k] <- Sigma
	}
	cat("\n")
	keep$beta <- keep$beta/n.mcmc
	cat("beta acceptance rate:",round(keep$beta,2)) 

	###
	### Write output
	###

	list(beta=beta.save,mu.beta=mu.beta.save,Sigma=Sigma.save,p=p.save,z=z.save,
		keep=keep,y=y,X=X,g=g,priors=priors,start=start,tune=tune,n.mcmc=n.mcmc)
}