\name{qpareto2}
\alias{qpareto2}
\title{Pareto distribution of the second kind - quantiles}
\usage{
  qpareto2(p, k = 1, s = 1)
}
\arguments{
  \item{p}{vector of probabilities, \eqn{p\in(0,1)}.}

  \item{k}{vector of shape parameters, \eqn{k>0}.}

  \item{s}{vector of scale parameters, \eqn{s>0}.}
}
\value{
  The function gives the theoretical quantiles.
}
\description{
  Quantile function for the Pareto Type-II (Lomax)
  distribution with shape parameter equal to \eqn{k>0} and
  scale parameter equal to \eqn{s>0}.
}
\seealso{
  \code{\link{dpareto2}}, \code{\link{ppareto2}},
  \code{\link{rpareto2}}, \code{\link{pareto2.zsestimate}},
  \code{\link{pareto2.mlekestimate}},
  \code{\link{pareto2.goftest}},
  \code{\link{pareto2.ftest}}
}

