\name{pareto2.mleksestimate}
\alias{pareto2.mleksestimate}
\title{Estimation of shape and scale parameters for the Pareto-II distribution (MLE)}
\usage{
  pareto2.mleksestimate(x, tol = 1e-20, smin = 1e-04,
    smax = 20)
}
\arguments{
  \item{x}{a non-negative numeric vector.}

  \item{tol}{the desired accuracy (convergence tolerance).}

  \item{smin}{lower bound for the scale parameter.}

  \item{smax}{upper bound for the scale parameter.}
}
\value{
  The list with the following components is passed as a
  result: \tabular{ll}{ \code{k} \tab the estimated
  parameter of shape.\cr \code{s} \tab the estimated
  parameter of scale.\cr } or \code{NA} if the maximum of
  the likelihood function does not exist.
}
\description{
  Finds the maximum likelihood estimator of the type II
  Pareto distribution's shape parameter \eqn{k} and scale
  parameter \eqn{s}.
}
\details{
  Note that the maximum of the likelihood function may not
  exist for some input vectors.
}
\seealso{
  \code{\link{dpareto2}}, \code{\link{pareto2.goftest}},
  \code{\link{pareto2.zsestimate}}
}

