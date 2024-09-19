#' @name BSTNB
#' @title Fit a Bayesian Spatiotemporal Negative Binomial model
#'
#' @description
#' Generate posterior samples for the parameters in a Bayesian Spatiotemporal Negative Binomial Model
#'
#' @usage BSTNB(y,Xtilde,A, oind=NULL,
#'             nchain=3,niter=1000,nburn=500,nthin=1)
#'
#' @param y vector of counts, must be non-negative
#' @param Xtilde matrix of covariates, numeric
#' @param A adjacency matrix, numeric
#' @param oind indices of offset
#' @param nchain positive integer, number of MCMC chains to be run
#' @param niter positive integer, number of iterations in each chain
#' @param nburn non-negative integer, number of iterations to be discarded as burn-in samples
#' @param nthin positive integer, thinning interval
#'
#' @import dplyr
#' @import mvtnorm
#' @import BayesLogit
#' @import spam
#' @import MCMCpack
#' @import msm
#' @import matrixcalc
#'
#' @return list of posterior samples of the parameters of the model
#' @export
BSTNB = function(y, Xtilde, A, oind = NULL, nchain=3, niter=1000, nburn=500, nthin=1){

  N = length(y)
  n <- nrow(A)			    # Number of spatial units
  nt <- N/n
  nis<-rep(nt,n) 		# Number of individuals per county; here it's balanced -- 50 per county per year
  # Note: may need to lower proposal variance, s, below as n_i increases
  sid<-rep(1:n,nis)
  tid<-rep(1:nis[1],n)
  t <- tid / max(tid)
  N<-length(sid) 		  # Total number of observations
  if(is.null(oind)){
    X = Xtilde
    x0 = rep(0,nrow(Xtilde))
  }else{
    X = matrix(Xtilde[,-c(oind)],nrow(Xtilde)) # Covariates
    x0 = Xtilde[,c(oind)]# Offset variable
  }
  p = ncol(X)

  ##########
  # Priors #
  ##########
  beta0<-rep(0,p)
  T0a<-diag(.01,p)
  T0b<-diag(.01,p)         # Uniform or Gamma(0.01,0.01) prior for r depending on MH or Gibbs
  s<-0.0003                # Proposal variance  -- NOTE: may need to lower this as n_i increases
  kappa<-.999999
  Q<-as.spam(diag(apply(A,1,sum)))-kappa*as.spam(A)

  ############
  # Num Sims #
  ############
  lastit<-(niter-nburn)/nthin	# Last stored value

  ############
  # Store #
  ############
  Beta<-array(NA,c(lastit,p,nchain))
  colnames(Beta) <- colnames(X)
  R<-R2<-matrix(0,lastit,nchain)
  Sigphi<-array(0,c(lastit,4,nchain))
  PHI3<-PHI4<-array(0,c(lastit,n,nchain))
  Eta<-array(0,c(lastit,N,nchain))

  for(chain in 1:3){

    #########
    # Inits #
    #########
    set.seed(2023+chain); beta<-rnorm(p)
    phi_init<-rmvnorm(1,sigma=diag(.1,2*n))	  # Random effects
    phi_init<-matrix(phi_init,ncol=2,byrow=T)  # n x 3 matrix of spatial effects
    phi3<-phi_init[,1]
    phi4<-phi_init[,2]
    phi3<-phi3-mean(phi3)
    phi4<-phi4-mean(phi4)
    Phi3<-rep(phi3,nis)
    Phi4<-rep(phi4,nis)

    phimat<-cbind(phi3,phi4)
    Sigmaphi<-cov(phimat)
    r<-1
    Acc<-0
    N0<-length(y[y==0])      # Number of observed 0's
    q<-rep(.5,N)             # 1-p=1/(1+exp(X%*%alpha)), used for updating y1


    ########
    # MCMC #
    ########

    for (i in 1:niter){

      # Update r
      rnew<-rtnorm(1,r,sqrt(s),lower=0)       # Treat r as continuous
      ratio<-sum(dnbinom(y,rnew,q,log=T))-sum(dnbinom(y,r,q,log=T))+
        dtnorm(r,rnew,sqrt(s),0,log=T) - dtnorm(rnew,r,sqrt(s),0,log=T)   # Uniform Prior for R
      # Proposal not symmetric
      if (log(runif(1))<ratio) {
        r<-rnew
        Acc<-Acc+1
      }

      # Update r2 using Gibbs as in Dadaneh et al and Zhou and Carin
      # Update latent counts, l
      l<-rep(0,N)
      ytmp<-y
      for(j in 1:N) l[j]<-sum(rbinom(ytmp[j],1,round(r/(r+1:ytmp[j]-1),6))) # Could try to avoid loop; rounding avoids numerical stability

      # Update r from conjugate gamma distribution given l and psi
      psi<-exp(eta2)/(1+exp(eta2))
      r2<-rgamma(1,0.01+sum(l),0.01-sum(log(1-psi)))

      # Update beta
      eta<-X%*%beta+x0+Phi3+Phi4*t
      w<-rpg(N,y+r,eta)                               # Polya weights
      z<-(y-r)/(2*w)
      v<-solve(crossprod(X*sqrt(w))+T0b)
      m<-v%*%(T0b%*%beta0+t(sqrt(w)*X)%*%(sqrt(w)*(z-x0-Phi3-Phi4*t)))
      beta<-c(rmvnorm(1,m,v))

      # Update phi3
      priorprec<-as.numeric(1/(Sigmaphi[1,1]-Sigmaphi[1,-1]%*%solve(Sigmaphi[-1,-1])%*%Sigmaphi[-1,1]))*Q # Prior Prec of phi3|phi1,phi2,phi4
      priormean<-diag(n)%x%(Sigmaphi[1,-1]%*%solve(Sigmaphi[-1,-1]))%*%c(t(phimat[,-1]))      # Prior mean of phi3|phi1,phi2,phi4
      prec<-priorprec+as.spam(diag(tapply(w,sid,sum),n,n))
      tmp<-tapply(w*(z-x0-Phi4*t),sid,sum)
      m<-c(priorprec%*%priormean)+tmp
      if(is.positive.definite(prec%>%as.matrix)) phi3<-rmvnorm.canonical(1, m, prec)[1,]

      # Center
      phi3<-phi3-mean(phi3)
      Phi3<-rep(phi3,nis)

      # Update phi4
      priorprec<-as.numeric(1/(Sigmaphi[2,2]-Sigmaphi[2,-2]%*%solve(Sigmaphi[-2,-2])%*%Sigmaphi[-2,2]))*Q # Prior Prec of phi4|phi1,phi2,phi3
      priormean<-diag(n)%x%(Sigmaphi[2,-2]%*%solve(Sigmaphi[-2,-2]))%*%c(t(phimat[,-2]))      # Prior mean of phi4|phi1,phi2,phi3

      prec<-priorprec+as.spam(diag(tapply(w*t^2,sid,sum),n,n))
      tmp<-tapply(t*w*(z-x0-Phi3),sid,sum)
      m<-c(priorprec%*%priormean)+tmp
      if(is.positive.definite(prec%>%as.matrix)) phi4<-rmvnorm.canonical(1, m, prec)[1,]

      # Center
      phi4<-phi4-mean(phi4)
      Phi4<-rep(phi4,nis)

      # Update Sigma.phi
      phimat<-cbind(phi3,phi4)
      try({
        Sigmaphi<-riwish(2+n-1,diag(2)+t(phimat)%*%Q%*%phimat)
      })

      # Store
      if (i> nburn & i%%nthin==0) {
        j<-(i-nburn)/nthin
        Beta[j,,chain]<-beta
        R[j,chain]<-r
        R2[j,chain]<-r2
        Sigphi[j,,chain]<-c(Sigmaphi)
        PHI3[j,,chain]<-phi3
        PHI4[j,,chain]<-phi4
        Eta[j,,chain]<-eta
      }

      # if (i%%10==0) print(paste(chain, "/", nchain,"chain | ",round(i/niter*100,2),"% completed"))
      if (i%%10==0) print(paste(chain, "/", nchain,"chain | ",round(i/niter*100,2),"% completed |","Test:",conv.test(R[,chain])))

    }
  }

  list.params = list(Alpha=NULL, Beta=Beta, R=R, R2=R2, Sigphi=Sigphi, PHI3=PHI3, PHI4=PHI4, Eta1=Eta)
  return(list.params)
}
