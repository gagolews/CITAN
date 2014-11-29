\name{discrpareto2.mleksestimate}
\alias{discrpareto2.mleksestimate}
\title{Estimation of shape and scale parameters for the Discretized Pareto-Type II distribution (MLE)}
\usage{
  discrpareto2.mleksestimate(x, kmin = 1e-04, kmax = 100,
    smin = 1e-04, smax = 100)
}
\arguments{
  \item{x}{a non-negative numeric vector.}

  \item{kmin}{lower bound for the shape parameter.}

  \item{kmax}{upper bound for the shape parameter.}

  \item{smin}{lower bound for the scale parameter.}

  \item{smax}{upper bound for the scale parameter.}
}
\value{
  The list with the following components is passed as a
  result: \tabular{ll}{ \code{k} \tab the estimated
  parameter of shape.\cr \code{s} \tab the estimated
  parameter of scale.\cr }
}
\description{
  Finds the maximum likelihood estimator of the Discretized
  type II Pareto distribution's shape parameter \eqn{k} and
  scale parameter \eqn{s}.
}
\details{
  If X has the Pareto-Type II distribution \eqn{P2(k,s)}
  then \code{Y=floor(X)} has the Discretized Pareto-Type II
  distribution DP2(k,s).
}
\seealso{
  \code{\link{ppareto2}},
  \code{\link{discrpareto2.mlekestimate}},\cr
  \code{\link{discrpareto2.goftest}}
}

