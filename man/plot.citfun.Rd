\name{plot.citfun}
\alias{plot.citfun}
\title{Draw the citation function of a given vector}
\usage{
  plot.citfun(x, ..., xmarg = 10, add = FALSE, xlab = "",
    ylab = "", main = "", ylim = c(0, max(x)),
    xlim = c(0, length(x) + xmarg))
}
\arguments{
  \item{x}{non-negative numeric vector.}

  \item{xmarg}{x-margin on the right.}

  \item{add}{logical; indicates whether to start a new
  plot, \code{FALSE} by default.}

  \item{ylim,xlim,xlab,ylab,main,...}{additional graphical
  parameters.}
}
\description{
  The citation function of a vector \eqn{x=(x_1,\dots,x_n)}
  is a mapping \deqn{\pi(y)=x_{(n-\lfloor
  y+1\rfloor+1)}}{\pi(y)=x_{(n-floor(y+1)+1)}} defined for
  \eqn{0\le y<n}, where \eqn{x_{(i)}} denotes the
  \eqn{i}-th smallest value of \eqn{x}.
}
\examples{
john_s <- c(11,5,4,4,3,2,2,2,2,2,1,1,1,0,0,0,0);
plot.citfun(john_s, main="Smith, John", col="red");
}
\seealso{
  \code{\link{curve.add.rp}}, \code{\link{curve.add.lp}},
  \code{\link{plot.default}}
}

