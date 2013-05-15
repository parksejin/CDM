

###########################################
# General Diagnostic Model 
###########################################
#
#..........................#
# OUTLINE:				   #
#..........................#
#
# Input:
# ------
# data ... polytomous responses
# 	I Items with categories 0,1,...,K
# group ... 1, ... , G
# 	covariates and groups. Up to now, only use groups
# theta.k ... Multidimensional ability design vector
#			as input. It is of dimension D and therefore a
# 			a [TH,D] matrix (TH is the number of theta points)
# 	Design matrix for smoothing of theta.k distribution
#
# Model:
# ------
# P(X_{pi}=k) = g( b_ik + a_i1 q_i1k theta_i1 + ... + a_iD q_iD*k theta_iD )
# b_ik ... item-category difficulty
# a_id ... item slopes
# q_ik ... design matrix for slopes
#	* in general it will be (0,1,...,K)
#     but can be specified by the user
#   * The theta vectors can also be transformed in the R formula
#     framework. Therefore D* will in general be larger
#     than D. For example, an interaction \theta_1*\theta_2 or
#     theta^2 can be used in the estimation.
#
# Algorithm:
# ----------
# b ... b[1:I,1:K]  matrix of difficulties
# a ... a[1:I,1:D] array of item slopes
# n_tikg n.ik[1:TH,1:I,1:K,1:G] array of expected item responses
#        of ability pattern t at item i in category k and group g
# pi.k ... [1:TH,1:G] marginal ability distribution in groups 1,...,G
#
# Specifications:
# ---------------
# o Q matrix: allocate items to dimensions:
#      can be either a matrix or an array
# o matrix with fixed b or a parameters
# o itemgroups for a and b parameters / constraints
#    est.a and est.b
# o log-linear smoothing of ability distribution
#     -> design matrix of theta
#
# Details:
# --------
# o If there are different categories per item,
#	then use the maximal category K for all items and 
#   set b_{ik} to infinity if for item i categeory k
#   is not observed.
# o Calculations of probabilities:
#    P(X_pi=k|\theta) = exp( b_{ik} + a_i1 * q_{i1k} theta_1 +
#					    	... + a_{iD} * q_{i1D} theta_D ) / NN
#	The normalization constraint NN must be calculated
#   appropriately .
#
#################################################################
# begin GDM function here
# 
gdm <- function( data, theta.k, irtmodel="2PL", group=NULL, weights=rep(1, nrow(data)), 
	Qmatrix=NULL ,thetaDes = NULL , skillspace="loglinear",
    b.constraint=NULL, a.constraint=NULL, 
	mean.constraint=NULL , Sigma.constraint=NULL , 
    delta.designmatrix=NULL, standardized.latent=FALSE , 
	centered.latent=FALSE , 
    maxiter=1000, conv=10^(-5), globconv=10^(-5), 
	decrease.increments = FALSE , 
	...){	
# mean.constraint [ dimension , group , value ]
# Sigma.constraint [ dimension1 , dimension2 , group , value ]		
	#*************************
	# data preparation
	s1 <- Sys.time()
	e1 <- environment()	
	## prevent from warnings in R CMD check "no visible binding"
	## gdm: no visible binding for global variable 'TD'
	TD <- TP <- EAP.rel <- mean.trait <- sd.trait <- skewness.trait <- NULL
	K.item <- correlation.trait <- NULL 	
	data <- as.matrix(data)
	dat.resp <- 1 - is.na(data)
	dat <- data
	dat[ is.na(data) ] <- 0
	# maximal categories
	K <- max(dat)
	# list of indicator data frames
	dat.ind <- as.list( 1:(K+1) )
	for (ii in 0:K){
		dat.ind[[ii+1]] <- 1 * ( dat==ii )*dat.resp
				}
	I <- ncol(dat)	# number of items
	n <- nrow(dat)
	
	# arrange groups
	if ( is.null(group)){ 
		G <- 1 
		group <- rep(1,n)
				} 
			else {
		gr2 <- unique( group )
		G <- length(gr2)
		group <- match( group , gr2 )
							}
	Ngroup <- aggregate( 1+0*group , list(group) , sum )[,2] 							
	
	# theta design
	res <- .gdm.thetadesign( theta.k , thetaDes , Qmatrix )
	.attach.environment( res , envir=e1 )

	
#print("p200")	
		
	#****
	# arrange b parameters and starting values
	b <- matrix( 0 , nrow=I , ncol=K )	
	for (ii in 1:K){	# ii <- 1
		cm1 <- colMeans( ( dat >= ii )*dat.resp )
		b[,ii] <-  qlogis( ( cm1 + .01 ) / 1.02 )
				}	
	# define some item parameter constraints here
#	if ( is.null( b.constraint ) & (skillspace=="loglinear") &
#	     	(irtmodel=="1PL") ){
#		b.constraint <- matrix( 0 , nrow=1,ncol=3)
#		b.constraint[ 1 , 1:2 ] <- c( which.min( abs( b ) )[1]  ,1)
#		b.constraint[ 1 , 3 ] <- 0
#					}

#print("p300")						

	#****
	# item slope matrix
	# a[1:I,1:TD] ... Items x theta dimension
	# a <- matrix( 1 , nrow=I , ncol=TD )
	# item x category slopes are in principle also possible
	KK <- K	# if KK == 1 then a slope parameter for all items is estimated
	a <- array( 1 , dim=c(I,TD,KK) )

	# define Q matrix
	res <- .gdm.Qmatrix(Qmatrix,irtmodel,I,TD,K,a)
	.attach.environment( res , envir=e1 )

	# constraints on item parameters
	res <- .gdm.constraints.itempars( b.constraint , a.constraint , 
			K , TD , Qmatrix , a)	
	.attach.environment( res , envir=e1 )	


	# starting values for distributions
	Sigma <- diag(1,D)
	library(mvtnorm)
	pik <- dmvnorm( matrix( theta.k ,ncol=D) , 
				mean=rep(0,D) , sigma = Sigma )
	pi.k <- matrix( 0 , TP , G )
	for (gg in 1:G){ pi.k[,gg] <- pik }
	n.ik <- array( 0 , dim=c(TP,I,K+1,G) )	
	
	#***
	# response patterns
    resp.ind.list <- list( 1:I )
	for (ii in 1:I){ resp.ind.list[[ii]] <- which( dat.resp[,ii] == 1)  }	
	
    #***
	# extract number of skills per dimensions
	skill.levels <- rep(0,D)
	for (dd in 1:D){ skill.levels[dd] <- length( unique(theta.k[,dd] ) ) }
	
	#****
	# create thetaDes design matrix for loglinear smoothing
	res <- .gdm.create.delta.designmatrix( delta.designmatrix , 
			TP , D , theta.k , skill.levels,G)	
	.attach.environment( res , envir=e1 )
		
	se.a <- 0*a
	max.increment.a <- .3
	max.increment.b <- 3
	
	if ( standardized.latent ){
		mean.constraint <- rbind( mean.constraint , cbind( 1:D , 1 , 0 )  )
		Sigma.constraint <- rbind( Sigma.constraint , cbind( 1:D , 1:D , 1 , 1 ) )
		skillspace <- "normal"
		}
	if ( centered.latent ){
		mean.constraint <- rbind( mean.constraint , cbind( 1:D , 1 , 0 )  )
		skillspace <- "normal"
		}

	
	#***
	# set constraints for a and b parameters if the maximal 
	# item category differs from item to item
	res <- .gdm.constraints.itempars2( b.constraint , a.constraint , K , TD , I ,
				dat )	
	.attach.environment( res , envir=e1 )
	
	#---
	# initial values algorithm
	dev <- 0	; iter <- 0
	globconv1 <- conv1 <- 1000
	disp <- paste( paste( rep(".", 100 ) , collapse="") ,"\n", sep="")
	
	############################################
	# BEGIN MML Algorithm
	############################################
		
	while( ( iter < maxiter ) & ( ( globconv1 > globconv) | ( conv1 > conv) ) ){
		
		#****
		# collect old parameters
		b0 <- b ; a0 <- a ; dev0 <- dev
		delta0 <- delta ; pi.k0 <- pi.k

# a0 <- Sys.time()
 		
		#****
		#1 calculate probabilities
		probs <- .gdm.calc.prob( a,b,thetaDes,Qmatrix,I,K,TP,TD)
 
# cat("calc.prob") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1		
 
		#*****
		#2 calculate individual likelihood
		h1 <- matrix( 1 , nrow=n , ncol=TP )
		res.hwt <- calc_posterior.v2(rprobs= probs , gwt=h1 , 
					 resp=dat , 
					 nitems= I , 
					 resp.ind.list=resp.ind.list , normalization=FALSE , 
					 thetasamp.density= NULL , snodes=0 )	
		p.xi.aj <- res.hwt$hwt 	

# cat("calc.post.v2") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1				
		
		#*****
		#3 calculate posterior and marginal distributions
		res <- .gdm.calc.post(pi.k,group,p.xi.aj,weights,G)
		p.aj.xi <- res$p.aj.xi
		pi.k <- res$pi.k

# cat("gdm.calc.post") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1				
		
		#*****
		#4 calculate expected counts
		# n.ik [ 1:TP , 1:I , 1:(K+1) , 1:G ]
		res <- .gdm.calc.counts(G,weights,dat.ind,dat,
					dat.resp,p.aj.xi,K,n.ik,TP,I,group)		
		n.ik <- res$n.ik
		N.ik <- res$N.ik

# cat("calc.counts") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1		

	
		#*****
		#5 M step: b parameter estimation
		# n.ik  [1:TP,1:I,1:K,1:G]
		# probs[1:I,1:K,1:TP]
		res <- .gdm.est.b(probs, n.ik, N.ik, I, K, G,b,b.constraint,
				max.increment=max.increment.b )	
		b <- res$b
		se.b <- res$se.b
		if (decrease.increments){ 	max.increment.b <- res$max.increment.b/1.1	}

# cat("est.b") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1		
		
		#*****
		#6 M step: a parameter estimation	
		if ( irtmodel == "2PL"){
			res <- .gdm.est.a(probs, n.ik, N.ik, I, K, G,a,a.constraint,TD,
					Qmatrix,thetaDes,TP, max.increment = max.increment.a)
			a <- res$a
			se.a <- res$se.a
			if (decrease.increments){ 	max.increment.a <- res$max.increment.a/1.1	}
						}
		if ( irtmodel == "2PLcat"){
			res <- .gdm.est.a.cat(probs, n.ik, N.ik, I, K, G,a,a.constraint,TD,
					Qmatrix,thetaDes,TP, max.increment = max.increment.a)					
			a <- res$a
			se.a <- res$se.a
			if (decrease.increments){ 	
				max.increment.a <- res$max.increment.a / 1.1				
					}
						}	
			
# print("a600")							
						
		#*****
		#7 M step: estimate reduced skillspace
		if ( skillspace == "loglinear" ){
			res <- .gdm.est.skillspace(Ngroup, pi.k , Z=delta.designmatrix ,
					G , delta  )
			pi.k <- res$pi.k
			delta <- res$delta
			covdelta <- res$covdelta
					}
		if ( skillspace == "normal" ){
			pi.k <- .gdm.est.normalskills( pi.k , theta.k , irtmodel,
						G,D , mean.constraint , Sigma.constraint ,
						standardized.latent	)
					}

# cat("calc skillspace") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1							

			
		#*****
		#8 calculate likelihood
		# n.ik [ TP , I , K+1 , G ]
		# N.ik [ TP , I , G ]
		# probs [I , K+1 , TP ]
		ll <- 0
		for (gg in 1:G){ 
			ind.gg <- which(group==gg)
			ll <- ll + sum( weights[ind.gg] * log( rowSums( p.xi.aj[ind.gg,] * 
						outer( rep(1,nrow(p.xi.aj[ind.gg,])) , pi.k[,gg] ) ) ) )
							}
		dev <- -2*ll	

# cat("calc LL") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1				

		#~~~~~~~~~~~~~~~~~~~~~~~~~~
		# display progress		
		gg0 <- abs( a - a0 )
		gg1 <- abs( b - b0 )
		pardiff <- max( gg1 , gg0 )
		deltadiff <- abs( pi.k - pi.k0 )	
		conv1 <- max( c(pardiff,deltadiff))
		globconv1 <- abs( dev - dev0) 
		iter <- iter +1
    	cat(disp)	
	    cat("Iteration" , iter , "   " , paste( Sys.time() ) , "\n" )		
		cat( paste( "   Deviance = "  , 
                   round( dev , 4 ) , 
 				  if (iter > 1 ){ " | Deviance change = " } else {""} ,
						 if( iter>1){round( - dev + dev0 , 6 )} else { ""}	,sep=""))
		if ( dev > dev0 & (iter>1 ) ){ cat( "  Deviance increases!") } ; cat("\n")
		cat( paste( "    Maximum item intercept parameter change = " , 
                             round( max( gg1) , 6 ) ,  " \n"   )  )  
		cat( paste( "    Maximum item slope parameter change = " , 
                             round( max( gg0 ) , 6 ) ,  " \n"   )  )  							 
		cat( paste( "    Maximum distribution parameter change = " , 
                             round( max( deltadiff ) , 6 ) ,  " \n"   )  )  
        flush.console()							 
								}
		############################################
		# END MML Algorithm
		############################################
		
		# collect item parameters
		res <- .gdm.collect.itempars( data , K , D , b , a , TD , thetaDes ,
					irtmodel , se.b , se.a )
		.attach.environment( res , envir=e1 )					
		
		# calculate distribution properties
		res <- .gdm.calc.distributionmoments( D , G , pi.k , theta.k )
		.attach.environment( res , envir=e1 )	
		
		# Information criteria
		ic <- .gdm.calc.ic( dev , dat , G , skillspace , irtmodel , 
				K,D,TD,I,b.constraint,a.constraint,mean.constraint ,
			    Sigma.constraint, delta.designmatrix , standardized.latent)

	#########################################
	# item fit [ items , theta , categories ] 
	# # n.ik [ 1:TP , 1:I , 1:(K+1) , 1:G ]
	probs <- aperm( probs , c(3,1,2) )
	itemfit.rmsea <- itemfit.rmsea( n.ik , pi.k , probs )		
    item$itemfit.rmsea <- itemfit.rmsea

	# person parameters
	res <- .gdm.person.parameters( data , D , theta.k , p.xi.aj , p.aj.xi , weights )	
	.attach.environment( res , envir=e1 )
	
	#*************************
	# collect output	
	s2 <- Sys.time()
	res <- list( "item" = item , "person" = person , "EAP.rel" = EAP.rel , 
				"deviance"=dev , "ic" = ic , "b" = b , "se.b" = se.b , 
				"a" = a ,  "se.a" = se.a , 
				"itemfit.rmsea" = itemfit.rmsea , 
				"mean.rmsea" = mean(itemfit.rmsea) , 				
				"Qmatrix"=Qmatrix , "pi.k"=pi.k , 	
				"mean.trait"=mean.trait , "sd.trait" = sd.trait , 
				"skewness.trait" = skewness.trait , "correlation.trait" = correlation.trait , 
				"pjk" = probs , "n.ik" = n.ik ,  
				"G"=G , "D"=D , "I" = ncol(data) , "N" = nrow(data) , 
				"delta" = delta , "covdelta"=covdelta , "data" = data )
	res$p.xi.aj <- p.xi.aj ; res$posterior <- p.aj.xi 
	res$skill.levels <- skill.levels
	res$K.item <- K.item
	res$theta.k <- theta.k
	res$thetaDes <- thetaDes
	res$time <- list( "s1"=s1,"s2"=s2 , "timediff"=s2-s1)
	res$skillspace <- skillspace
	res$iter <- iter
	
	# some further values for modelfit.gdm
	res$AIC <- res$ic$AIC
	res$BIC <- res$ic$BIC	
	res$Npars <- res$ic$np	
	res$loglike <- - res$deviance / 2
	
                cat("----------------------------------- \n")
                cat("Start:" , paste( s1) , "\n")
                cat("End:" , paste(s2) , "\n")
                cat("Difference:" , print(s2 -s1), "\n")
                cat("----------------------------------- \n")
	class(res) <- "gdm"
	return(res)
				
		}
