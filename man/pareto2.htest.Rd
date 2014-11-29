\name{pareto2.htest}
\alias{pareto2.htest}
\title{Two-sample h-test for equality of shape parameters for Type II-Pareto distributions with known common scale parameter}
\usage{
  pareto2.htest(x, y, s,
    alternative = c("two.sided", "less", "greater"),
    significance = 0.05, wyg = NULL, verbose = TRUE,
    drho = 0.005, K = NULL, improve = TRUE)
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

  \item{significance}{significance level. See Value for
  details.}

  \item{wyg}{precomputed h-dependent acceptation region or
  \code{NULL}. See Value for details.}

  \item{verbose}{logical; if \code{TRUE} then the
  computation progress will be printed out.}

  \item{drho}{power calculation accuracy, a single number
  in [0.001, 0.1]. The smaller the value the slower
  computation, but more precise. This is used to determine
  \code{K} iff \code{K} is not given.}

  \item{K}{numeric vector; shape parameters for which to
  calculate the power function or \code{NULL}.}

  \item{improve}{logical; if \code{TRUE} then the greedy
  heuristic algorithm for improving the acceptation region
  will be run.}
}
\value{
  The list of class \code{power.htest} with the following
  components is passed as a result: \tabular{ll}{
  \code{statistic} \tab the value of the test statistic.\cr
  \code{result} \tab either FALSE (accept null hypothesis)
  or TRUE (reject).\cr \code{alternative} \tab a character
  string describing the alternative hypothesis.\cr
  \code{method} \tab a character string indicating what
  type of test was performed.\cr \code{data.name} \tab a
  character string giving the name(s) of the data.\cr
  \code{wyg} \tab a numeric vector giving the h-dependent
  acceptation region used.\cr \code{size} \tab size of the
  test corresponding to \code{wyg}.\cr \code{qual} \tab
  quality of the test corresponding to \code{wyg}, the
  closer to \code{significance}, the better.\cr } Currently
  no method for determining the p-value of this test is
  implemented.
}
\description{
  Performs \eqn{h}-test for equality of shape parameters of
  two samples from the Pareto type-II distributions with
  known and equal scale parameters, \eqn{s>0}.
}
\details{
  Given two equal-sized samples \eqn{X=(X_1,...,X_n)}
  i.i.d. \eqn{P2(k_x,s)} and \eqn{Y=(Y_1,...,Y_m)} i.i.d.
  \eqn{P2(k_y,s)} this test verifies the null hypothesis
  \eqn{H_0: k_x=k_y} against two-sided or one-sided
  alternatives, depending on the value of
  \code{alternative}. It bases on test statistic
  \code{T=H(Y)-H(X)} where \eqn{H} denotes Hirsch's
  \eqn{h}-index (see \code{\link{index.h}}).

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
  \code{\link{pareto2.ftest}},
  \code{\link{pareto2.htest.approx}}, \code{\link{index.h}}
}

