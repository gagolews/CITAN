\name{pareto2.htest.approx}
\alias{pareto2.htest.approx}
\title{Two-sample asymptotic h-test for equality of shape parameters for Type II-Pareto distributions with known common scale parameter}
\usage{
  pareto2.htest.approx(x, y, s,
    alternative = c("two.sided", "less", "greater"),
    significance = NULL)
}
\arguments{
  \item{x}{an n-element non-negative numeric vector of data
  values.}

  \item{y}{an n-element non-negative numeric vector of data
  values.}

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
  Performs asymptotic (approximate) \eqn{h}-test for
  equality of shape parameters of two samples from the
  Pareto type-II distributions with known and equal scale
  parameters, \eqn{s>0}.
}
\details{
  Given two equal-sized samples \eqn{X=(X_1,...,X_n)}
  i.i.d. \eqn{P2(k_x,s)} and \eqn{Y=(Y_1,...,Y_m)} i.i.d.
  \eqn{P2(k_y,s)} this test verifies the null hypothesis
  \eqn{H_0: k_x=k_y} against two-sided or one-sided
  alternatives, depending on the value of
  \code{alternative}. It bases on a test statistic that is
  a function of {H(Y)-H(X)}, where \eqn{H} denotes Hirsch's
  \eqn{h}-index (see \code{\link{index.h}}). This statistic
  approximately has asymptotically standardized normal
  distribution under \eqn{H_0}.

  Note that for \eqn{k_x < k_y}, then \eqn{X} dominates
  \eqn{Y} stochastically.
}
\references{
  Gagolewski M., Grzegorzewski P., S-Statistics and Their
  Basic Properties, In: Borgelt C. et al (Eds.), Combining
  Soft Computing and Statistical Methods in Data Analysis,
  Springer-Verlag, 2010, 281-288.\cr
}
\seealso{
  \code{\link{dpareto2}}, \code{\link{pareto2.goftest}},
  \code{\link{pareto2.ftest}}, \code{\link{pareto2.htest}},
  \code{\link{index.h}}
}

