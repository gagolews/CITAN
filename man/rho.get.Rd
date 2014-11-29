\name{rho.get}
\alias{rho.get}
\title{Rho-index of a continuous probability distribution}
\usage{
  rho.get(cdf, kappa, ..., tol = 1e-20)
}
\arguments{
  \item{cdf}{a cumulative distribution function, e.g.
  \code{\link{ppareto2}}.}

  \item{kappa}{an increasing function, \eqn{\kappa} (see
  Details), a so-called control function.}

  \item{...}{optional arguments to \code{cdf}.}

  \item{tol}{the desired accuracy (convergence tolerance).}
}
\value{
  The function returns a single number.
}
\description{
  Numerically computes the rho-index of a given continuous
  cumulative distribution function
}
\details{
  Let \eqn{F} be a continuous c.d.f that is strictly
  increasing on \eqn{[a,b]}, where \eqn{a=\inf\{x:
  F(x)>0\}}{a=inf{x: F(x)>0}} and \eqn{b=\sup\{x:
  F(x)<1\}}{b=sup{x: F(x)<1}}.

  A \dfn{control function} is any function
  \eqn{\kappa:[0,1]\to[c,d]\subseteq[a,b]}{\kappa:[0,1]->[c,d]c[a,b]}
  that is continuous and strictly increasing and which
  fulfills \eqn{\kappa(0)=c} and \eqn{\kappa(1)=d}.

  The \dfn{\eqn{\rho}-index} of the distribution \eqn{F}
  (Gagolewski, Grzegorzewski, 2010) is a number
  \eqn{\rho_\kappa\in(0,1)}{0<\rho_\kappa<1} such that
  \deqn{\rho_\kappa=1-F(\kappa(\rho_\kappa)).}

  It turns out that under certain conditions in a model of
  i.i.d. random variables the S-statistic associated with
  \eqn{\kappa} is an asymptotically unbiased, normal and
  strongly consistent estimator of \eqn{\rho_\kappa}.
}
\examples{
kappa <- function(x) { pmax(0,pmin(x,1)) } # identity function on [0,1]
rho.get(ppareto2, kappa, 1, 1)             # golden ratio
}
\references{
  Gagolewski M., Grzegorzewski P., S-Statistics and Their
  Basic Properties, In: Borgelt C. et al (Eds.), Combining
  Soft Computing and Statistical Methods in Data Analysis,
  Springer-Verlag, 2010, 281-288.\cr
}
\seealso{
  \code{\link{phirsch}}, \code{\link{dhirsch}},
  \code{\link{psstat}}, \code{\link{dsstat}},
  \code{\link{Sstat}}, \code{\link{Sstat2}}
}

