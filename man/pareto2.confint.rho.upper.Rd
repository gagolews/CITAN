\name{pareto2.confint.rho.upper}
\alias{pareto2.confint.rho.upper}
\title{Right-sided exact confidence interval for the rho-index}
\usage{
  pareto2.confint.rho.upper(v, kappa, s, n,
    conf.level = 0.95, tol = 1e-12)
}
\arguments{
  \item{v}{observed value of the S-statistic w.r.t.
  \eqn{\kappa}.}

  \item{kappa}{an increasing function, \eqn{\kappa}, a
  so-called control function.}

  \item{s}{scale parameter, \eqn{s>0}.}

  \item{n}{sample size.}

  \item{conf.level}{confidence level; defaults 0.95.}

  \item{tol}{the desired accuracy (convergence tolerance).}
}
\value{
  Upper bound of the confidence interval.
}
\description{
  #/internal/ #/internal/ Computes the exact right-sided
  confidence interval for the \eqn{\rho}-index of a
  probability distribution in an \eqn{(X_1,\dots,X_n)}
  i.i.d. Pareto-type II model with known scale parameter
  \eqn{s>0}. The confidence interval bases on the observed
  value of S-statistic w.r.t. to the given control function
  \eqn{\kappa}.
}
\details{
  For more information see man page on
  \code{\link{rho.get}}, \code{\link{Sstat}} and the paper
  (Gagolewski, Grzegorzewski, 2010).
}
\references{
  Gagolewski M., Grzegorzewski P., S-Statistics and Their
  Basic Properties, In: Borgelt C. et al (Eds.), Combining
  Soft Computing and Statistical Methods in Data Analysis,
  Springer-Verlag, 2010, 281-288.\cr
}
\seealso{
  \code{\link{ppareto2}},
  \code{\link{pareto2.confint.rho.approx}},
  \code{\link{Sstat}},
  \code{\link{pareto2.confint.rho.lower}},
  \code{\link{pareto2.confint.rho}}, \code{\link{rho.get}}
}

