% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BSTNB.R
\name{BSTNB}
\alias{BSTNB}
\title{Fit a Bayesian Spatiotemporal Negative Binomial model}
\usage{
BSTNB(y,X,A,
            nchain=3,niter=100,nburn=20,nthin=1)
}
\arguments{
\item{y}{vector of counts, must be non-negative}

\item{X}{matrix of covariates, numeric}

\item{A}{adjacency matrix, numeric}

\item{nchain}{positive integer, number of MCMC chains to be run}

\item{niter}{positive integer, number of iterations in each chain}

\item{nburn}{non-negative integer, number of iterations to be discarded as burn-in samples}

\item{nthin}{positive integer, thinning interval}
}
\value{
list of posterior samples of the parameters of the model
}
\description{
Generate posterior samples for the parameters in a Bayesian Spatiotemporal Negative Binomial Model
}
\examples{
data(simdat)
y <- simdat$y
X <- cbind(simdat$V1,simdat$x)
data(county.adjacency)
data(USAcities)
IAcities <- subset(USAcities,state_id=="IA")
countyname <- unique(IAcities$county_name)
A <- get_adj_mat(county.adjacency,countyname,c("IA"))
\donttest{
res2 <- BSTNB(y, X, A, nchain=2, niter=100, nburn=20, nthin=1)
}

}
