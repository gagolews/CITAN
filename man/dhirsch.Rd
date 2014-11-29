\name{dhirsch}
\alias{dhirsch}
\title{Distribution of the h-index - p.m.f.}
\usage{
  dhirsch(x, n, cdf, ...)
}
\arguments{
  \item{x}{numeric vector.}

  \item{n}{sample size.}

  \item{cdf}{a continuous cumulative distribution function
  \eqn{F}, e.g. \code{\link{ppareto2}}.}

  \item{...}{optional arguments to \code{cdf}.}
}
\value{
  The value of the p.m.f. at \code{x}.
}
\description{
  The probability mass function of Hirsch's \eqn{h}-index
  for sample of size \code{n} in an i.i.d. model with
  common increasing and continuous c.d.f. \eqn{F} defined
  on \eqn{[0,\infty)}.
}
\references{
  Gagolewski M., Grzegorzewski P., S-Statistics and Their
  Basic Properties, In: Borgelt C. et al (Eds.), Combining
  Soft Computing and Statistical Methods in Data Analysis,
  Springer-Verlag, 2010, 281-288.\cr
}
\seealso{
  \code{\link{index.h}}, \code{\link{rho.get}},
  \code{\link{phirsch}}, \code{\link{pareto2.dhirsch}}
}

