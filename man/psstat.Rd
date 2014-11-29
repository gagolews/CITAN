\name{psstat}
\alias{psstat}
\title{Distribution of S-statistics - c.d.f.}
\usage{
  psstat(x, n, cdf, kappa, ...)
}
\arguments{
  \item{x}{numeric vector.}

  \item{n}{sample size.}

  \item{cdf}{a cumulative distribution function \eqn{F},
  e.g. \code{\link{ppareto2}}.}

  \item{kappa}{an increasing function, \eqn{\kappa} (see
  Details), a so-called control function.}

  \item{...}{optional arguments to \code{cdf} and
  \code{pdf}.}
}
\value{
  The value of the c.d.f. at \code{x}.
}
\description{
  Computes the cumulative distribution function of the
  S-statistic w.r.t. to a control function in an i.i.d.
  model with common increasing and continuous c.d.f.
  \eqn{F} defined on \eqn{[0,\infty)}.
}
\details{
  Let \eqn{F} (parameter \code{cdf}) be a continuous c.d.f
  that is strictly increasing on \eqn{[a,b]}, where
  \eqn{a=\inf\{x: F(x)>0\}}{a=inf{x: F(x)>0}} and
  \eqn{b=\sup\{x: F(x)<1\}}{b=sup{x: F(x)<1}}.

  Moreover, let
  \eqn{\kappa:[0,1]\to[c,d]\subseteq[a,b]}{\kappa:[0,1]->[c,d]c[a,b]}
  be any control function (parameter \code{kappa}), i.e. a
  function that is continuous and strictly increasing and
  which fulfills \eqn{\kappa(0)=c} and \eqn{\kappa(1)=d}.

  The function computes the value of the c.d.f. of an
  S-statistic w.r.t. to the control function for sample of
  size \code{n}. This result was given in (Gagolewski,
  Grzegorzewski, 2010).

  Note that under certain conditions the distribution of an
  S-statistic is asymptotically normal with expectation
  \eqn{\rho_\kappa} (see \code{\link{rho.get}}) and
  variance \eqn{\rho_\kappa
  (1-\rho_\kappa)/n/(1+F(\kappa(\rho_\kappa)))^2}.
}
\references{
  Gagolewski M., Grzegorzewski P., S-Statistics and Their
  Basic Properties, In: Borgelt C. et al (Eds.), Combining
  Soft Computing and Statistical Methods in Data Analysis,
  Springer-Verlag, 2010, 281-288.\cr
}
\seealso{
  \code{\link{Sstat}}, \code{\link{Sstat2}},
  \code{\link{dsstat}}, \code{\link{rho.get}}
}

