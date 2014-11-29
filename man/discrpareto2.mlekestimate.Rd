\name{discrpareto2.mlekestimate}
\alias{discrpareto2.mlekestimate}
\title{Estimation of shape parameter for the Discretized Pareto-Type II distribution (MLE)}
\usage{
  discrpareto2.mlekestimate(x, s, kmin = 1e-04, kmax = 100)
}
\arguments{
  \item{x}{a non-negative numeric vector.}

  \item{s}{scale parameter, \eqn{s>0}.}

  \item{kmin}{lower bound for the shape parameter.}

  \item{kmax}{upper bound for the shape parameter.}
}
\value{
  A single numeric value is returned, the ML estimator of
  \eqn{k}.
}
\description{
  Finds the maximum likelihood estimator of the Discretized
  type II Pareto distribution's shape parameter \eqn{k} for
  given scale parameter \eqn{s}.
}
\details{
  If X has the Pareto-Type II distribution \eqn{P2(k,s)}
  then \code{Y=floor(X)} has the discretized Pareto-Type II
  distribution DP2(k,s).
}
\seealso{
  \code{\link{ppareto2}},
  \code{\link{discrpareto2.mleksestimate}},\cr
  \code{\link{discrpareto2.goftest}}
}

