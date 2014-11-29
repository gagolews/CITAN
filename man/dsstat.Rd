\name{dsstat}
\alias{dsstat}
\title{Distribution of S-statistics - density}
\usage{
  dsstat(x, n, cdf, pdf, kappa, kappaInvDer, ...)
}
\arguments{
  \item{x}{numeric vector.}

  \item{n}{sample size.}

  \item{cdf}{a cumulative distribution function \eqn{F},
  e.g. \code{\link{ppareto2}}.}

  \item{pdf}{a density function \eqn{F'}, e.g.
  \code{\link{dpareto2}}.}

  \item{kappa}{an increasing function, \eqn{\kappa}, a
  so-called control function.}

  \item{kappaInvDer}{the derivative of the inverse of
  \eqn{\kappa}.}

  \item{...}{optional arguments to \code{cdf} and
  \code{pdf}.}
}
\value{
  The value of the density at \code{x}.
}
\description{
  Computes the density function of the S-statistic w.r.t.
  to a control function in an i.i.d. model with common
  increasing and continuous c.d.f. \eqn{F} defined on
  \eqn{[0,\infty)}.
}
\details{
  The function computes the value of the p.d.f. of an
  S-statistic w.r.t. to the control function \code{kappa}
  for sample of size \code{n}. Note that the result is
  valid only at continuity points of \eqn{F'}.

  For more information see man page on \code{\link{psstat}}
  and the paper (Gagolewski, Grzegorzewski, 2010).
}
\references{
  Gagolewski M., Grzegorzewski P., S-Statistics and Their
  Basic Properties, In: Borgelt C. et al (Eds.), Combining
  Soft Computing and Statistical Methods in Data Analysis,
  Springer-Verlag, 2010, 281-288.\cr
}
\seealso{
  \code{\link{Sstat}}, \code{\link{Sstat2}},
  \code{\link{psstat}}, \code{\link{rho.get}}
}

