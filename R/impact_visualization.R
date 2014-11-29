## This file is part of the CITAN library.
##
## Copyright 2011 Marek Gagolewski
##
##
## CITAN is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## CITAN is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
## GNU Lesser General Public License for more details.
##
## You should have received a copy of the GNU Lesser General Public License
## along with CITAN. If not, see <http://www.gnu.org/licenses/>.


#' The citation function of a vector \eqn{x=(x_1,\dots,x_n)}
#' is a mapping \deqn{\pi(y)=x_{(n-\lfloor y+1\rfloor+1)}}{\pi(y)=x_{(n-floor(y+1)+1)}}
#' defined for \eqn{0\le y<n}, where \eqn{x_{(i)}} denotes the
#' \eqn{i}-th smallest value of \eqn{x}.
#'
#' @title Draw the citation function of a given vector
#' @param x non-negative numeric vector.
#' @param xmarg x-margin on the right.
#' @param add logical; indicates whether to start a new plot, \code{FALSE} by default.
#' @param ylim,xlim,xlab,ylab,main,... additional graphical parameters.
#' @seealso \code{\link{curve.add.rp}}, \code{\link{curve.add.lp}}, \code{\link{plot.default}}
#' @examples
#' john_s <- c(11,5,4,4,3,2,2,2,2,2,1,1,1,0,0,0,0);
#' plot.citfun(john_s, main="Smith, John", col="red");
#' @export
plot.citfun <- function(x, ..., xmarg=10, add=FALSE, xlab="", ylab="", main="", ylim=c(0, max(x)), xlim=c(0,length(x)+xmarg))
{
   n <- length(x);
   if (n == 0) stop("Zero-length vector given.");
   if (mode(x) != "numeric") stop("Non-numeric vector given.");
   if (any(x < 0)) warning("The function should be used with non-negative vectors.");

   x <- sort(x, decreasing=T);

   px <- n;
   py <- x[n];

   for (i in (n-1):1)
   {
      if (x[i]!=x[i+1])
      {
         px <- c(i,    px);
         py <- c(x[i], py);
      }
   }



   if (!add)
   {
      plot(px, py, ylim=ylim, xlim=xlim, type='p', xlab=xlab, ylab=ylab, main=main, ...);
   } else
   {
      points(px, py, ...);
   }

# 	segments(px, py, px, c(py[-1],0), ...);
   segments(c(0,px[-length(px)]), py, px, py, ...);
}




#' The \eqn{r_p}-curve appears in the definition of the \eqn{r_p}-index
#' (see Gagolewski, Grzegorzewski, 2009) and the \code{\link{index.rp}} function.
#'
#' @title Draw the r_p-curve of given radius
#' @param r radius of the \eqn{r_p}-curve; \eqn{r>0}.
#' @param p index order, \eqn{p \in [1,\infty]}{p in [1,\infty]}; defaults \eqn{\infty} (\code{Inf}).
#' @param n integer; the maximal number of values at which to evaluate the underlying function.
#' @param ... additional graphical parameters.
#' @seealso \code{\link{index.rp}}, \code{\link{plot.citfun}}, \code{\link{plot.default}}
#' @examples
#' john_s <- c(11,5,4,4,3,2,2,2,2,2,1,1,1,0,0,0,0);
#' plot.citfun(john_s, main="Smith, John");
#' curve.add.rp(index.rp(john_s), col="green");
#' curve.add.rp(index.rp(john_s,1), 1, col="blue");
#' curve.add.rp(index.rp(john_s,2), 2, col="red");
#' @references
#' Gagolewski M., Grzegorzewski P., A geometric approach to the construction of scientific impact indices, Scientometrics, 81(3), 2009a, 617-634.\cr
#' @export
curve.add.rp <- function(r, p=Inf, n=101, ...)
{
   if (length(p) != 1 || mode(p) != "numeric") stop("p must be a single numeric value");
   if (p < 1) stop("p must be >= 1");

   if (length(r) != 1 || mode(r) != "numeric") stop("r must be a single numeric value");
   if (r <= 0) stop("r must be > 0");

   if (is.finite(p))
   {
      px <- seq(0,r,length=n);
      py <- (r^p-px^p)^(1.0/p);
      lines(px, py, ...);
   } else
   {
      lines(c(0,r), c(r,r), ...);
      points(r, r, ...);
   }
}




#' The \eqn{l_p}-curve appears in the definition of the \eqn{l_p}-index
#' (see Gagolewski, Grzegorzewski, 2009) and the \code{\link{index.lp}} function.
#'
#' @title Draw the l_p-curve of given size
#' @param ab size of the \eqn{l_p}-curve; positive numeric vector of length 2.
#' @param p index order, \eqn{p \in [1,\infty]}{p in [1,\infty]}; defaults \eqn{\infty} (\code{Inf}).
#' @param n integer; the maximal number of values at which to evaluate the underlying function.
#' @param ... additional graphical parameters.
#' @seealso \code{\link{index.lp}}, \code{\link{plot.citfun}}, \code{\link{plot.default}}
#' @examples
#' john_s <- c(11,5,4,4,3,2,2,2,2,2,1,1,1,0,0,0,0);
#' plot.citfun(john_s, main="Smith, John");
#' curve.add.lp(index.lp(john_s), col="green");
#' curve.add.lp(index.lp(john_s,1), 1, col="blue");
#' curve.add.lp(index.lp(john_s,2), 2, col="red");
#' @references
#' Gagolewski M., Grzegorzewski P., A geometric approach to the construction of scientific impact indices, Scientometrics, 81(3), 2009a, 617-634.\cr
#' @export
curve.add.lp <- function(ab, p=Inf, n=101, ...)
{
   if (length(p) != 1 || mode(p) != "numeric") stop("p must be a single numeric value");
   if (p < 1) stop("p must be >= 1");

   if (length(ab) != 2 || mode(ab) != "numeric") stop("ab must be a numeric vector of length 2");
   if (any(ab <= 0)) stop("ab must be > 0");

   a <- ab[1];
   b <- ab[2];

   if (is.finite(p))
   {
      px <- seq(0,a,length=n);
      py <- (b^p-(b/a*px)^p)^(1.0/p);
      lines(px, py, ...);
   } else
   {
      lines(c(0,a), c(b,b), ...);
      points(a, b, ...);
   }
}
