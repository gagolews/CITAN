\name{pareto2.confint.h.upper}
\alias{pareto2.confint.h.upper}
\title{Right-sided exact confidence interval for the theoretical h-index}
\usage{
  pareto2.confint.h.upper(h, s, n, conf.level = 0.95,
    tol = 1e-12)
}
\arguments{
  \item{h}{observed value of the \eqn{h}-index}

  \item{s}{scale parameter, \eqn{s>0}.}

  \item{n}{sample size.}

  \item{conf.level}{confidence level; defaults 0.95.}

  \item{tol}{the desired accuracy (convergence tolerance).}
}
\value{
  Upper bound of the confidence interval.
}
\description{
  Computes the exact right-sided confidence interval for
  the theoretical \eqn{h}-index of a probability
  distribution in an \eqn{(X_1,\dots,X_n)} i.i.d.
  Pareto-type II model with known scale parameter
  \eqn{s>0}.
}
\details{
  See \code{\link{pareto2.confint.h}} for details.
}
\references{
  Gagolewski M., Grzegorzewski P., S-Statistics and Their
  Basic Properties, In: Borgelt C. et al (Eds.), Combining
  Soft Computing and Statistical Methods in Data Analysis,
  Springer-Verlag, 2010, 281-288.\cr
}
\seealso{
  \code{\link{index.h}}, \code{\link{ppareto2}},
  \code{\link{rho.get}},
  \code{\link{pareto2.confint.rho}},\cr
  \code{\link{pareto2.confint.h}},
  \code{\link{pareto2.confint.h.lower}}
}

