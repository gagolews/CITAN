\name{pareto2.mlekestimate}
\alias{pareto2.mlekestimate}
\title{Estimation of shape parameter for the Pareto-II distribution (MLE)}
\usage{
  pareto2.mlekestimate(x, s)
}
\arguments{
  \item{x}{a non-negative numeric vector.}

  \item{s}{scale parameter, \eqn{s>0}.}
}
\value{
  A single numeric value is returned, the unbiased ML
  estimator of \eqn{k}.
}
\description{
  Finds the unbiased maximum likelihood estimator of the
  type II Pareto distribution's shape parameter \eqn{k} for
  known scale parameter \eqn{s}.
}
\seealso{
  \code{\link{dpareto2}}, \code{\link{pareto2.goftest}},
  \code{\link{pareto2.zsestimate}},
  \code{\link{pareto2.mleksestimate}}
}

