\name{pareto2.confint.h}
\alias{pareto2.confint.h}
\title{Two-sided exact confidence interval for the theoretical h-index}
\usage{
  pareto2.confint.h(h, s, n, conf.level = 0.95,
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
  Vector of length 2 with the computed bounds of the
  confidence interval.
}
\description{
  Computes the exact two-sided confidence interval for the
  theoretical \eqn{h}-index of a probability distribution
  in an \eqn{(X_1,\dots,X_n)} i.i.d. Pareto-type II model
  with known scale parameter \eqn{s>0}.
}
\details{
  The \dfn{Theoretical \eqn{h}-index} for a sequence of
  \eqn{n} i.i.d. random variables with common increasing
  and continuous c.d.f. \eqn{F} defined on \eqn{[0,\infty)}
  is equal to \eqn{n\varrho_\kappa}{n*\rho_\kappa}, where
  \eqn{\rho_\kappa} is the \eqn{\rho}-index of \eqn{F} for
  \eqn{\kappa(x)=nx}, see \code{\link{rho.get}} for
  details.
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
  \code{\link{pareto2.confint.h.upper}},
  \code{\link{pareto2.confint.h.upper}}
}

