#' @name BSTP
#' @title Fit a Bayesian Spatiotemporal Poisson model
#'
#' @description
#' Generate posterior samples for the parameters in a Bayesian Spatiotemporal Poisson Model
#'
#' @usage BSTP(y,X,A,nt,
#'             nchain=3,nsim=100,nburn=20,nthin=1)
#'
#' @param y vector of counts, must be non-negative
#' @param X matrix of covariates, numeric
#' @param A adjacency matrix, numeric
#' @param nt positive integer, number of time points
#' @param nchain positive integer, number of MCMC chains to be run
#' @param nsim positive integer, number of iterations in each chain
#' @param nburn non-negative integer, number of iterations to be discarded as burn-in samples
#' @param nthin positive integer, thinning interval
#'
#' @importFrom stats cov
#' @importFrom stats dnbinom
#' @importFrom stats rbinom
#' @importFrom stats reorder
#' @importFrom stats rgamma
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom stats spline
#' @importFrom stats var
#' @import dplyr
#' @import BayesLogit
#' @import spam
#' @import MCMCpack
#' @importFrom matrixcalc is.positive.definite
#'
#' @return list of posterior samples of the parameters of the model
#' @noRd
BSTP <- function(y, X, A, nt, nchain=3, nsim=100, nburn=20, nthin=1){

  ## Run the necessary checks
  if(!is.vector(y)){stop("y must be a vector")}
  if(!is.matrix(X)){stop("X must be a matrix")}
  if(!is.matrix(A)){stop("A must be a matrix")}
  if(nt <= 0){stop("nt must be positive integer")}
  if(nchain < 1){stop("nchain must be a positive integer")}
  if(nsim < 1){stop("nsim must be a positive integer")}
  if(nburn < 0){stop("nburn must be a non-negative integer")}
  if(nthin < 1){stop("nthin must be a positive integer")}
  y <- as.numeric(y)
  if(min(y,na.rm=T)<0){stop("y must be non-negative")}
  if(!is.numeric(X)){stop("X must be numeric")}
  if(!is.numeric(A)){stop("A must be numeric")}


  n <- nrow(A)			    # Number of spatial units
  nis <- rep(nt,n) 		# Number of individuals per county; here it's balanced -- 50 per county per year
  # Note: may need to lower proposal variance, s, below as n_i increases
  sid <- rep(1:n,nis)
  tid <- rep(1:nis[1],n)
  t <- tid / max(tid)
  N <- length(sid) 		  # Total number of observations
  p <- ncol(X)

  ##########
  # Priors #
  ##########
  beta0 <- rep(0,p)
  T0a <- diag(.01,p)
  T0b <- diag(.01,p)         # Uniform or Gamma(0.01,0.01) prior for r depending on MH or Gibbs
  s <- 0.0003                # Proposal variance  -- NOTE: may need to lower this as n_i increases
  kappa <- 0.999999
  Q <- as.spam(diag(apply(A,1,sum)))-kappa*as.spam(A)

  ############
  # Num Sims #
  ############
  lastit <- (nsim-nburn)/nthin	# Last stored value

  ############
  # Store #
  ############
  Beta <- array(0,c(lastit,p,nchain))
  colnames(Beta) <- colnames(X)
  R <- R2 <- matrix(0,lastit,nchain)
  Sigphi <- array(0,c(lastit,4,nchain))
  PHI3 <- PHI4 <- array(0,c(lastit,n,nchain))
  Eta <- array(0,c(lastit,N,nchain))

  for(chain in 1:3){

    #########
    # Inits #
    #########
    set.seed(2023+chain); beta <- rnorm(p)
    phi_init <- spam::rmvnorm(1,sigma=diag(.1,2*n))	  # Random effects
    phi_init <- matrix(phi_init,ncol=2,byrow=T)  # n x 3 matrix of spatial effects
    phi3 <- phi_init[,1]
    phi4 <- phi_init[,2]
    phi3 <- phi3-mean(phi3)
    phi4 <- phi4-mean(phi4)
    Phi3 <- rep(phi3,nis)
    Phi4 <- rep(phi4,nis)

    phimat <- cbind(phi3,phi4)
    Sigmaphi <- cov(phimat)
    r <- 1
    Acc <- 0
    N0 <- length(y[y==0])      # Number of observed 0's
    q <- rep(.5,N)             # 1-p=1/(1+exp(X%*%alpha)), used for updating y1


    ########
    # MCMC #
    ########

    for (i in 1:nsim){

      # Fix r
      r <- 1
      r2 <- 1

      # Update beta
      eta <- X%*%beta+Phi3+Phi4*t
      w <- rpg(N,y+r,eta)                               # Polya weights
      z <- (y-r)/(2*w)
      v <- solve(crossprod(X*sqrt(w))+T0b)
      m <- v%*%(T0b%*%beta0+t(sqrt(w)*X)%*%(sqrt(w)*(z-Phi3-Phi4*t)))
      beta <- c(spam::rmvnorm(1,m,v))

      # Update phi3
      priorprec <- as.numeric(1/(Sigmaphi[1,1]-Sigmaphi[1,-1]%*%solve(Sigmaphi[-1,-1])%*%Sigmaphi[-1,1]))*Q # Prior Prec of phi3|phi1,phi2,phi4
      priormean <- diag(n)%x%(Sigmaphi[1,-1]%*%solve(Sigmaphi[-1,-1]))%*%c(t(phimat[,-1]))      # Prior mean of phi3|phi1,phi2,phi4
      prec <- priorprec+as.spam(diag(tapply(w,sid,sum),n,n))
      tmp <- tapply(w*(z-X%*%beta-Phi4*t),sid,sum)
      m <- c(priorprec%*%priormean)+tmp
      if(is.positive.definite(prec%>%as.matrix)) phi3 <- spam::rmvnorm.canonical(1, m, prec)[1,]

      # Center
      phi3 <- phi3-mean(phi3)
      Phi3 <- rep(phi3,nis)

      # Update phi4
      priorprec <- as.numeric(1/(Sigmaphi[2,2]-Sigmaphi[2,-2]%*%solve(Sigmaphi[-2,-2])%*%Sigmaphi[-2,2]))*Q # Prior Prec of phi4|phi1,phi2,phi3
      priormean <- diag(n)%x%(Sigmaphi[2,-2]%*%solve(Sigmaphi[-2,-2]))%*%c(t(phimat[,-2]))      # Prior mean of phi4|phi1,phi2,phi3

      prec <- priorprec+as.spam(diag(tapply(w*t^2,sid,sum),n,n))
      tmp <- tapply(t*w*(z-X%*%beta-Phi3),sid,sum)
      m <- c(priorprec%*%priormean)+tmp
      if(is.positive.definite(prec%>%as.matrix)) phi4 <- spam::rmvnorm.canonical(1, m, prec)[1,]

      # Center
      phi4 <- phi4-mean(phi4)
      Phi4 <- rep(phi4,nis)

      # Update Sigma.phi
      phimat <- cbind(phi3,phi4)
      try({
        Sigmaphi <- riwish(2+n-1,diag(2)+t(phimat)%*%Q%*%phimat)
      })

      # Store
      if (i > nburn & i%%nthin==0) {
        j <- (i-nburn)/nthin
        Beta[j,,chain] <- beta
        R[j,chain] <- r
        R2[j,chain] <- r2
        Sigphi[j,,chain] <- c(Sigmaphi)
        PHI3[j,,chain] <- phi3
        PHI4[j,,chain] <- phi4
        Eta[j,,chain] <- eta
      }

      # if (i%%10==0) print(paste(chain, "/", nchain,"chain | ",round(i/nsim*100,2),"% completed"))
      if (i%%10==0) print(paste(chain, "/", nchain,"chain | ",round(i/nsim*100,2),"% completed |","Test:",conv.test(R[,chain])))
    }
  }

  list_params <- list(Beta=Beta, Sigphi=Sigphi, R=R, PHI3=PHI3, PHI4=PHI4, Eta1=Eta)
  return(list_params)
}
