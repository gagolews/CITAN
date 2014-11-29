\name{rpareto2}
\alias{rpareto2}
\title{Pareto distribution of the second kind - random deviates}
\usage{
  rpareto2(n, k = 1, s = 1)
}
\arguments{
  \item{n}{number of observations.}

  \item{k}{vector of shape parameters, \eqn{k>0}.}

  \item{s}{vector of scale parameters, \eqn{s>0}.}
}
\value{
  The function returns generated pseudorandom deviates.
}
\description{
  Random deviates generation for the Pareto Type-II (Lomax)
  distribution with shape parameter equal to \eqn{k>0} and
  scale parameter equal to \eqn{s>0}.
}
\seealso{
  \code{\link{dpareto2}}, \code{\link{ppareto2}},
  \code{\link{qpareto2}}, \code{\link{pareto2.zsestimate}},
  \code{\link{pareto2.mlekestimate}},
  \code{\link{pareto2.goftest}},
  \code{\link{pareto2.ftest}}
}

