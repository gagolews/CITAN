\name{pareto2.goftest}
\alias{pareto2.goftest}
\title{Goodness-of-fit test for the Pareto-II distribution}
\usage{
  pareto2.goftest(x, k = NULL, s = NULL,
    method = c("anderson-darling", "kolmogorov"))
}
\arguments{
  \item{x}{a non-negative numeric vector of data values.}

  \item{k}{scale parameter, \eqn{k>0} or \code{NULL}.}

  \item{s}{shape parameter, \eqn{s>0} or \code{NULL}.}

  \item{method}{either "anderson-darling" or "kolmogorov".}
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
  Performs goodness-of-fit test for the Pareto-II
  distribution basing on MLE or MMSE estimates (Zhang,
  Stevens, 2009) and the Anderson-Darling or Kolmogorov
  test.
}
\details{
  This method, proposed e.g. by Zhang and Stevens (2009),
  uses either the function \code{\link[ADGofTest]{ad.test}}
  from package \pkg{ADGofTest} or \code{\link{ks.test}} to
  compute the selected test.

  It is known that the tests have low powers.

  If \code{k} and \code{s} are NULL, it bases on bayesian
  MMS estimators, see \code{\link{pareto2.zsestimate}}. If
  \code{s} is not NULL, then the unbiased maximum
  likelihood estimator is used to determine the scale
  parameter (see \code{\link{pareto2.mlekestimate}}) iff it
  is not given.
}
\references{
  Zhang J., Stevens M.A., A New and Efficient Estimation
  Method for the Generalized Pareto Distribution,
  Technometrics 51(3), 2009, 316-325.\cr
}
\seealso{
  \code{\link{dpareto2}}, \code{\link{pareto2.zsestimate}},
  \code{\link{pareto2.mlekestimate}},
  \code{\link{ks.test}}, \code{\link[ADGofTest]{ad.test}}
  from package \code{ADGofTest}
}

