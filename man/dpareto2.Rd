\name{dpareto2}
\alias{dpareto2}
\title{Pareto distribution of the second kind - density}
\usage{
  dpareto2(x, k = 1, s = 1)
}
\arguments{
  \item{x}{vector of quantiles.}

  \item{k}{vector of shape parameters, \eqn{k>0}.}

  \item{s}{vector of scale parameters, \eqn{s>0}.}
}
\value{
  The function gives the density.
}
\description{
  Density function for the Pareto Type-II (Lomax)
  distribution with shape parameter equal to \eqn{k>0} and
  scale parameter equal to \eqn{s>0}.
}
\details{
  The p.d.f. at \eqn{x\ge 0} is given by
  \code{f(x)=k*s^k/(s+x)^(k+1)}.
}
\seealso{
  \code{\link{ppareto2}}, \code{\link{qpareto2}},
  \code{\link{rpareto2}}, \code{\link{pareto2.zsestimate}},
  \code{\link{pareto2.mlekestimate}},
  \code{\link{pareto2.goftest}},
  \code{\link{pareto2.ftest}}
}

