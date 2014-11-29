\name{pareto2.dhirsch}
\alias{pareto2.dhirsch}
\title{Distribution of the h-index - p.m.f.}
\usage{
  pareto2.dhirsch(x, n, k, s)
}
\arguments{
  \item{x}{numeric vector.}

  \item{n}{sample size.}

  \item{k}{shape parameter, \eqn{k>0}.}

  \item{s}{scale parameter, \eqn{s>0}.}
}
\value{
  The value of the p.m.f. at \code{x}.
}
\description{
  The probability mass function of Hirsch's \eqn{h}-index
  for sample of size \code{n} in an i.i.d. Pareto-Type II
  model.
}
\details{
  This is a significantly faster (written in C) version of
  the more general function \code{\link{dhirsch}}.
}
\references{
  Gagolewski M., Grzegorzewski P., S-Statistics and Their
  Basic Properties, In: Borgelt C. et al (Eds.), Combining
  Soft Computing and Statistical Methods in Data Analysis,
  Springer-Verlag, 2010, 281-288.\cr
}
\seealso{
  \code{\link{index.h}}, \code{\link{rho.get}},
  \code{\link{dhirsch}}, \code{\link{phirsch}},
  \code{\link{pareto2.phirsch}}
}

