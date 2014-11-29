\name{ppareto2}
\alias{ppareto2}
\title{Pareto distribution of the second kind - c.d.f.}
\usage{
  ppareto2(q, k = 1, s = 1)
}
\arguments{
  \item{q}{vector of quantiles.}

  \item{k}{vector of shape parameters, \eqn{k>0}.}

  \item{s}{vector of scale parameters, \eqn{s>0}.}
}
\value{
  The function gives the c.d.f..
}
\description{
  Cumulative distribution function for the Pareto Type-II
  (Lomax) distribution with shape parameter equal to
  \eqn{k>0} and scale parameter equal to \eqn{s>0}.
}
\details{
  The c.d.f. at \eqn{x\ge 0} is given by
  \code{F(x)=1-s^k/(s+x)^k}.
}
\seealso{
  \code{\link{dpareto2}}, \code{\link{qpareto2}},
  \code{\link{rpareto2}}, \code{\link{pareto2.zsestimate}},
  \code{\link{pareto2.mlekestimate}},
  \code{\link{pareto2.goftest}},
  \code{\link{pareto2.ftest}}
}

