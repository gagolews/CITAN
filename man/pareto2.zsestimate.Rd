\name{pareto2.zsestimate}
\alias{pareto2.zsestimate}
\title{Estimation of parameters for the Pareto-II distribution (MMSE)}
\usage{
  pareto2.zsestimate(x)
}
\arguments{
  \item{x}{a non-negative numeric vector.}
}
\value{
  The list with the following components is passed as a
  result: \tabular{ll}{ \code{k} \tab the estimated
  parameter of shape.\cr \code{s} \tab the estimated
  parameter of scale.\cr }
}
\description{
  Finds the MMS estimator of the type II Pareto
  distribution parameters using the Bayesian method (and
  the R code) developed by Zhang and Stevens (2009).
}
\references{
  Zhang J., Stevens M.A., A New and Efficient Estimation
  Method for the Generalized Pareto Distribution,
  Technometrics 51(3), 2009, 316-325.\cr
}
\seealso{
  \code{\link{dpareto2}}, \code{\link{pareto2.goftest}},
  \code{\link{pareto2.mlekestimate}},
  \code{\link{pareto2.mleksestimate}}
}

