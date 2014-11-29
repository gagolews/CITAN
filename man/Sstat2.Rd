\name{Sstat2}
\alias{Sstat2}
\title{S-statistic2}
\usage{
  Sstat2(x, kappaInv, sorted.dec = FALSE,
    disable.check = FALSE)
}
\arguments{
  \item{x}{a vector of real numbers.}

  \item{kappaInv}{a nondecreasing function ranging on
  [0,1], \eqn{\gamma} (see Details), the inverse of a
  so-called control function.}

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
  Computes the S-statistic2 w.r.t. to the identity function
  for data transformed by the inverse of a control
  function.
}
\details{
  Given a sequence of \eqn{n} non-negative numbers
  \eqn{x=(x_1,\dots,x_n)}, where \eqn{x_i\ge x_j} for
  \eqn{i\le j}, and a nondecreasing function \eqn{\gamma:
  R\to[0,1]}{\gamma: R->[0,1]}, the \dfn{S-statistic2}
  (Gagolewski, Grzegorzewski, 2010) for \eqn{x} is defined
  as \deqn{V_n(x)=\max_{i=1,\dots,n}\{\min\{\gamma(x_i),
  i/n \}\}}{V_n(x)=max{ min{i/n, \gamma(x_i)} } for
  i=1,\dots,n}

  If \code{disable.check} is set to \code{FALSE}, then
  eventual \code{NA} values are removed from the input
  vector.

  If a non-increasingly sorted vector is given as input
  (set \code{sorted.dec} to \code{TRUE}) the result is
  computed in linear time.
}
\examples{
x <- rpareto2(25, 1.05, 1);
kappaInv <- function(x) { pmax(0,pmin(1,x/25)); }
Sstat2(x, kappaInv, FALSE, TRUE);
}
\references{
  Gagolewski M., Grzegorzewski P., S-Statistics and Their
  Basic Properties, In: Borgelt C. et al (Eds.), Combining
  Soft Computing and Statistical Methods in Data Analysis,
  (SMPS 2010), Springer-Verlag, 2010, 281-288.
}
\seealso{
  \code{\link{index.h}}, \code{\link{index.g}},
  \code{\link{index.rp}}, \code{\link{index.lp}},
  \code{\link{Sstat}}, \code{\link{psstat}},
  \code{\link{dsstat}}
}

