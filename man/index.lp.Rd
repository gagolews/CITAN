\name{index.lp}
\alias{index.lp}
\title{The l_p-index}
\usage{
  index.lp(x, p = Inf, sorted.dec = FALSE,
    disable.check = FALSE)
}
\arguments{
  \item{x}{a non-negative numeric vector.}

  \item{p}{index order, \eqn{p \in [1,\infty]}{p in
  [1,\infty]}; defaults \eqn{\infty} (\code{Inf}).}

  \item{sorted.dec}{logical; \code{TRUE} if the vector has
  already been sorted non-increasingly; defaults
  \code{FALSE}.}

  \item{disable.check}{logical; \code{TRUE} to disable some
  validity checks on the input vector; defaults
  \code{FALSE}.}
}
\value{
  The function returns a numeric vector of length 2 equal
  to \eqn{(i,x_i)} or NA if improper input has been given.
}
\description{
  Computes the \eqn{l_p}-index of a numeric vector for
  given \eqn{p}.
}
\details{
  Given a sequence of \eqn{n} non-negative numbers
  \eqn{x=(x_1,\dots,x_n)}, where \eqn{x_i \ge x_j} for
  \eqn{i \le j}, the \dfn{\eqn{l_p}-index} for
  \eqn{p=\infty} equals to \deqn{l_p(x)=\arg\max_{(i,x_i),
  i=1,\dots,n} \{ i x_i \}}{l_p(x) = arg max_(i,x_i) {
  i*x_i } for i=1,\dots,n} if \eqn{n \ge 1}, or
  \eqn{l_\infty(x)=0} otherwise. Note that if
  \eqn{(i,x_i)=l_\infty(x)}, then \deqn{MAXPROD(x) = i
  x_i,}{MAXPROD(x) = i*x_i,} where \eqn{MAXPROD} is the
  index proposed in (Kosmulski, 2007).

  For the definition of the \eqn{l_p}-index for \eqn{p <
  \infty} we refer to (Gagolewski, Grzegorzewski, 2009a).

  If \code{disable.check} is set to \code{FALSE}, then
  eventual \code{NA} values are removed from the input
  vector.

  If a non-increasingly sorted vector is given as input
  (set \code{sorted.dec} to \code{TRUE}) the result is
  computed in linear time (see Gagolewski, Debski,
  Nowakiewicz, 2009b).
}
\examples{
x <- runif(100, 0, 100);
index.lp(x);                # two-dimensional value, can not be used
                            # directly in the analysis
prod(index.lp(x));          # the MAXPROD-index (one-dimensional)
mean(index.lp(x,1));        # some other one-dimensional impact index
}
\references{
  Gagolewski M., Grzegorzewski P., A geometric approach to
  the construction of scientific impact indices,
  Scientometrics, 81(3), 2009a, pp. 617-634.\cr Gagolewski
  M., Debski M., Nowakiewicz M., Efficient algorithms for
  computing ''geometric'' scientific impact indices,
  Research Report of Systems Research Institute, Polish
  Academy of Sciences RB/1/2009, 2009b.\cr Kosmulski M.,
  MAXPROD - A new index for assessment of the scientific
  output of an individual, and a comparison with the
  h-index, Cybermetrics, 11(1), 2007.\cr
}
\seealso{
  \code{\link{index.h}}, \code{\link{index.g}},
  \code{\link{index.rp}}, \code{\link{Sstat}},
  \code{\link{Sstat2}}
}

