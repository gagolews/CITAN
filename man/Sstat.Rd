\name{Sstat}
\alias{Sstat}
\title{S-statistic}
\usage{
  Sstat(x, kappa, sorted.dec = FALSE,
    disable.check = FALSE)
}
\arguments{
  \item{x}{a vector of real numbers.}

  \item{kappa}{an increasing function, \eqn{\kappa} (see
  Details), a so-called control function.}

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
  Computes the S-statistic w.r.t. to a control function.
}
\details{
  Given a sequence of \eqn{n} non-negative numbers
  \eqn{x=(x_1,\dots,x_n)}, where \eqn{x_i\ge x_j} for
  \eqn{i\le j}, and an increasing function \eqn{\kappa:
  [0,1]\to[a,b]}{\kappa: [0,1]->[a,b]} for some \eqn{a,b},
  the \dfn{S-statistic} (Gagolewski, Grzegorzewski, 2010)
  w.r.t. \eqn{\kappa} for \eqn{x} is defined as
  \deqn{V_n(x)=\max_{i=1,\dots,n}\{\min\{x_i, \kappa(i/n)
  \}\}}{V_n(x)=max{ min{\kappa(i/n), x_i} } for
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
kappa <- function(x) { pmax(0,pmin(1,x))*25; }
Sstat(x, kappa, FALSE, TRUE);
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
  \code{\link{Sstat2}}, \code{\link{psstat}},
  \code{\link{dsstat}}
}

