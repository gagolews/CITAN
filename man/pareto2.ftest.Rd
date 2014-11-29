\name{pareto2.ftest}
\alias{pareto2.ftest}
\title{Two-sample F-test for equality of shape parameters for Type II-Pareto distributions with known common scale parameter}
\usage{
  pareto2.ftest(x, y, s,
    alternative = c("two.sided", "less", "greater"),
    significance = NULL)
}
\arguments{
  \item{x}{a non-negative numeric vector of data values.}

  \item{y}{a non-negative numeric vector of data values.}

  \item{s}{scale parameter, \eqn{s>0}.}

  \item{alternative}{indicates the alternative hypothesis
  and must be one of "two.sided" (default), "less", or
  "greater".}

  \item{significance}{significance level,
  \eqn{0<}\code{significance}\eqn{<1} or \code{NULL}. See
  Value for details.}
}
\value{
  If \code{significance} is not \code{NULL}, then the list
  of class \code{power.htest} with the following components
  is passed as a result: \tabular{ll}{ \code{statistic}
  \tab the value of the test statistic.\cr \code{result}
  \tab either FALSE (accept null hypothesis) or TRUE
  (reject).\cr \code{alternative} \tab a character string
  describing the alternative hypothesis.\cr \code{method}
  \tab a character string indicating what type of test was
  performed.\cr \code{data.name} \tab a character string
  giving the name(s) of the data.\cr } Otherwise, the list
  of class \code{htest} with the following components is
  passed as a result: \tabular{ll}{ \code{statistic} \tab
  the value of the test statistic.\cr \code{p.value} \tab
  the p-value of the test.\cr \code{alternative} \tab a
  character string describing the alternative
  hypothesis.\cr \code{method} \tab a character string
  indicating what type of test was performed.\cr
  \code{data.name} \tab a character string giving the
  name(s) of the data.\cr }
}
\description{
  Performs F-test for equality of shape parameters of two
  samples from the Pareto type-II distributions with known
  and equal scale parameters, \eqn{s>0}.
}
\details{
  Given two samples \eqn{(X_1,...,X_n)} i.i.d.
  \eqn{P2(k_x,s)} and \eqn{(Y_1,...,Y_m)} i.i.d.
  \eqn{P2(k_y,s)} this test verifies the null hypothesis
  \eqn{H_0: k_x=k_y} against two-sided or one-sided
  alternatives, depending on the value of
  \code{alternative}. It bases on test statistic
  \code{T=n/m*sum(log(1+Y/m))/sum(log(1+X/n))} which, under
  \eqn{H_0}, has the Snedecor's F distribution with
  \eqn{(2m, 2n)} degrees of freedom.

  Note that for \eqn{k_x < k_y}, then \eqn{X} dominates
  \eqn{Y} stochastically.
}
\seealso{
  \code{\link{dpareto2}}, \code{\link{pareto2.goftest}},
  \code{\link{pareto2.htest}},
  \code{\link{pareto2.htest.approx}}
}

