\name{index.rp}
\alias{index.rp}
\title{The r_p-index}
\usage{
  index.rp(x, p = Inf, sorted.dec = FALSE,
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
  The function returns a single number or NA if improper
  input has been given.
}
\description{
  Computes the \eqn{r_p}-index of a numeric vector for
  given \eqn{p}.
}
\details{
  Given a sequence of \eqn{n} non-negative numbers
  \eqn{x=(x_1,\dots,x_n)}, where \eqn{x_i \ge x_j} for
  \eqn{i \le j}, the \dfn{\eqn{r_p}-index} for
  \eqn{p=\infty} equals to \deqn{r_p(x)=\max_{i=1,\dots,n}
  \{ \min\{i,x_i\} \}}{r_p(x) = max{ min{i, x_i} } for
  i=1,\dots,n} if \eqn{n \ge 1}, or \eqn{r_\infty(x)=0}
  otherwise. For the definition of the \eqn{r_p}-index for
  \eqn{p < \infty} we refer to (Gagolewski, Grzegorzewski,
  2009).

  Note that if \eqn{x_1,\dots,x_n} are integers, then
  \deqn{r_\infty(x)=H(x),} where \eqn{H} is the
  \eqn{h}-index (Hirsch, 2005) and \deqn{r_1(x)=W(x),}
  where \eqn{W} is the \eqn{w}-index (Woeginger, 2008).

  If \code{disable.check} is set to \code{FALSE}, then
  eventual \code{NA} values are removed from the input
  vector.

  If a non-increasingly sorted vector is given as input
  (set \code{sorted.dec} to \code{TRUE}) the value of
  \eqn{r_\infty} is calculated in log time (note that it
  may be determined in linear time using \code{max(pmin(x,
  1:length(x)))}). Otherwise, linear time is needed.
}
\examples{
x <- runif(100, 0, 100);
index.rp(x);            # the r_oo-index
floor(index.rp(x));     # the h-index
index.rp(floor(x), 1);  # the w-index
}
\references{
  Gagolewski M., Grzegorzewski P., A geometric approach to
  the construction of scientific impact indices,
  Scientometrics, 81(3), 2009, pp. 617-634.\cr Hirsch J.E.,
  An index to quantify individual's scientific research
  output, Proceedings of the National Academy of Sciences
  102(46), 16569-16572, 2005.\cr Woeginger G.J., An
  axiomatic characterization of the Hirsch-index,
  Mathematical Social Sciences, 56(2), 224-232, 2008.\cr
}
\seealso{
  \code{\link{index.h}}, \code{\link{index.g}},
  \code{\link{index.lp}}, \code{\link{Sstat}},
  \code{\link{Sstat2}}
}

