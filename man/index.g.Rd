\name{index.g}
\alias{index.g}
\title{Egghe's g-index}
\usage{
  index.g(x, sorted.dec = FALSE, disable.check = FALSE)
}
\arguments{
  \item{x}{a non-negative numeric vector.}

  \item{sorted.dec}{logical; \code{TRUE} if the vector has
  already been sorted non-increasingly; defaults
  \code{FALSE}.}

  \item{disable.check}{logical; \code{TRUE} to disable some
  validity checks on the input vector; defaults
  \code{FALSE}.}
}
\value{
  The function returns a single number or NA if improper
  input has been given.
}
\description{
  Computes the "Classical" \eqn{g}-index of a numeric
  vector.
}
\details{
  Given a sequence of \eqn{n} non-negative numbers
  \eqn{x=(x_1,\dots,x_n)}, where \eqn{x_i \ge x_j} for
  \eqn{i \le j}, the \dfn{\eqn{g}-index} (Egghe, 2006) for
  \eqn{x} is defined as \deqn{G(x)=\max\{i=1,\dots,n:
  \sum_{j=1}^i x_i \ge i^2\},}{G(x)=max{i=1,\dots,n:
  x_1+\dots+x_i \ge i^2}} if \eqn{n \ge 1} and \eqn{x_1 \ge
  1}, or \eqn{G(x)=0} otherwise.

  If \code{disable.check} is set to \code{FALSE}, then
  eventual \code{NA} values are removed from the input
  vector.

  If a non-increasingly sorted vector is given as input
  (set \code{sorted.dec} to \code{TRUE}) the value is
  calculated in linear time.
}
\references{
  Egghe L., Theory and practise of the g-index,
  Scientometrics 69(1), 131-152, 2006.\cr
}
\seealso{
  \code{\link{index.h}}, \code{\link{index.rp}},
  \code{\link{index.lp}}, \code{\link{Sstat}},
  \code{\link{Sstat2}}
}

