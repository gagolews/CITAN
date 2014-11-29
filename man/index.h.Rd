\name{index.h}
\alias{index.h}
\title{Hirsch's h-index}
\usage{
  index.h(x, sorted.dec = FALSE, disable.check = FALSE,
    algorithm = c("log-time", "linear-time"))
}
\arguments{
  \item{x}{a non-negative numeric vector.}

  \item{sorted.dec}{logical; \code{TRUE} if the vector has
  already been sorted non-increasingly; defaults
  \code{FALSE}.}

  \item{disable.check}{logical; \code{TRUE} to disable some
  validity checks on the input vector; defaults
  \code{FALSE}.}

  \item{algorithm}{type of algorithm, "linear-time" or
  "log-time" (default).}
}
\value{
  The function returns a single number or NA if improper
  input has been given.
}
\description{
  Computes the "Classical" \eqn{h}-index of a numeric
  vector.
}
\details{
  Given a sequence of \eqn{n} non-negative numbers
  \eqn{x=(x_1,\dots,x_n)}, where \eqn{x_i \ge x_j} for
  \eqn{i \le j}, the \dfn{\eqn{h}-index} (Hirsch, 2005) for
  \eqn{x} is defined as \deqn{H(x)=\max\{i=1,\dots,n: x_i
  \ge i\}}{H(x)=max{i=1,\dots,n: x_i \ge i}} if \eqn{n \ge
  1} and \eqn{x_1 \ge 1}, or \eqn{H(x)=0} otherwise.

  If \code{disable.check} is set to \code{FALSE}, then
  eventual \code{NA} values are removed from the input
  vector.

  If a non-increasingly sorted vector is given as input
  (set \code{sorted.dec} to \code{TRUE}) the value is
  calculated in linear or log-time, depending on the value
  of the \code{algorithm} parameter.
}
\examples{
authors <- list(  # a list of numeric sequences
                  # (e.g. citation counts of the articles
                  # written by some authors)
    "A" =c(23,21,4,2,1,0,0),
    "B" =c(11,5,4,4,3,2,2,2,2,2,1,1,1,0,0,0,0),
    "C" =c(53,43,32,23,14,13,12,8,4,3,2,1,0)
 );
index.h(authors$A);
lapply(authors, index.h);
}
\references{
  Hirsch J.E., An index to quantify individual's scientific
  research output, Proceedings of the National Academy of
  Sciences 102(46), 16569-16572, 2005.\cr
}
\seealso{
  \code{\link{pareto2.confint.h}},
  \code{\link{pareto2.htest}},
  \code{\link{pareto2.htest.approx}},\cr
  \code{\link{phirsch}}, \code{\link{dhirsch}},
  \code{\link{index.g}}, \code{\link{index.rp}},
  \code{\link{index.lp}}, \code{\link{Sstat}},
  \code{\link{Sstat2}}
}

