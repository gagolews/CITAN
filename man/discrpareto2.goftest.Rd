\name{discrpareto2.goftest}
\alias{discrpareto2.goftest}
\title{Goodness-of-fit test for the Discretized Pareto-II distribution}
\usage{
  discrpareto2.goftest(x, k = NULL, s = NULL, kmin = 1e-04,
    kmax = 100, smin = 1e-04, smax = 100)
}
\arguments{
  \item{x}{a non-negative numeric vector of data values of
  length >= 9.}

  \item{k}{scale parameter, \eqn{k>0} or \code{NULL}.}

  \item{s}{shape parameter, \eqn{s>0} or \code{NULL}.}

  \item{kmin}{lower bound for the shape parameter.}

  \item{kmax}{upper bound for the shape parameter.}

  \item{smin}{lower bound for the scale parameter.}

  \item{smax}{upper bound for the scale parameter.}
}
\value{
  The list of class \code{htest} with the following
  components is passed as a result: \tabular{ll}{
  \code{statistic} \tab the value of the test statistic.\cr
  \code{p.value} \tab the p-value of the test.\cr
  \code{alternative} \tab a character string describing the
  alternative hypothesis.\cr \code{method} \tab a character
  string indicating what type of test was performed.\cr
  \code{data.name} \tab a character string giving the
  name(s) of the data.\cr }
}
\description{
  Performs goodness-of-fit test for the Discretized
  Pareto-II distribution basing on MLE estimates and the
  chi-square test.
}
\details{
  If X has the Pareto-Type II distribution \eqn{P2(k,s)}
  then \code{Y=floor(X)} has the Discretized Pareto-Type II
  distribution DP2(k,s).

  It is known that the test has low power.
}
\seealso{
  \code{\link{discrpareto2.mleksestimate}},
  \code{\link{discrpareto2.mlekestimate}},
  \code{\link{chisq.test}}
}

