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


#' The probability mass function of Hirsch's \eqn{h}-index
#' for sample of size \code{n} in an i.i.d. model with common increasing and continuous c.d.f. \eqn{F} defined on \eqn{[0,\infty)}.
#'
#' @references
#' Gagolewski M., Grzegorzewski P., S-Statistics and Their Basic Properties, In: Borgelt C. et al (Eds.), Combining Soft Computing and Statistical Methods in Data Analysis, Springer-Verlag, 2010, 281-288.\cr
#'
#' @title Distribution of the h-index - p.m.f.
#' @param x numeric vector.
#' @param n sample size.
#' @param cdf a continuous cumulative distribution function \eqn{F}, e.g. \code{\link{ppareto2}}.
#' @param ... optional arguments to \code{cdf}.
#' @return The value of the p.m.f. at \code{x}.
#' @export
#' @seealso \code{\link{index.h}}, \code{\link{rho.get}}, \code{\link{phirsch}}, \code{\link{pareto2.dhirsch}}
dhirsch <- function(x, n, cdf, ...)
{
	phirsch(x+1e-9, n, cdf, ...)-phirsch(x-1e-9, n, cdf, ...)
}


#' The cumulative distribution function of Hirsch's \eqn{h}-index
#' for sample of size \code{n} in an i.i.d. model with common increasing and continuous c.d.f. \eqn{F} defined on \eqn{[0,\infty)}.
#'
#' @references
#' Gagolewski M., Grzegorzewski P., S-Statistics and Their Basic Properties, In: Borgelt C. et al (Eds.), Combining Soft Computing and Statistical Methods in Data Analysis, Springer-Verlag, 2010, 281-288.\cr
#'
#' @title Distribution of the h-index - c.d.f.
#' @param x numeric vector.
#' @param n sample size.
#' @param cdf a continuous cumulative distribution function \eqn{F}, e.g. \code{\link{ppareto2}}.
#' @param ... optional arguments to \code{cdf}.
#' @return The value of the c.d.f. at \code{x}.
#' @export
#' @seealso \code{\link{index.h}}, \code{\link{rho.get}}, \code{\link{dhirsch}}, \code{\link{pareto2.phirsch}}
phirsch <- function(x, n, cdf, ...)
{
	warn <- getOption("warn");
	options("warn"=-1);
	y <- ifelse(x>n-1e-9, 1.0,
	     ifelse(x<0, 0.0,
	     pbeta(cdf(floor(x+1), ...), n-floor(x),floor(x)+1)));
	options("warn"=warn);
	return(y);
}





#' The cumulative distribution function of Hirsch's \eqn{h}-index
#' for sample of size \code{n} in an i.i.d. Pareto-Type II model.
#'
#' This is a significantly faster (written in C) version of the more general function \code{\link{phirsch}}.
#'
#' @references
#' Gagolewski M., Grzegorzewski P., S-Statistics and Their Basic Properties, In: Borgelt C. et al (Eds.), Combining Soft Computing and Statistical Methods in Data Analysis, Springer-Verlag, 2010, 281-288.\cr
#'
#' @title Distribution of the h-index - c.d.f.
#' @param x numeric vector.
#' @param n sample size.
#' @param k shape parameter, \eqn{k>0}.
#' @param s scale parameter, \eqn{s>0}.
#' @return The value of the c.d.f. at \code{x}.
#' @export
#' @seealso \code{\link{index.h}}, \code{\link{rho.get}}, \code{\link{phirsch}}, \code{\link{dhirsch}}, \code{\link{pareto2.dhirsch}}
pareto2.phirsch <- function(x, n, k, s)
{
	m <- as.integer(length(x));
# 	.C("pareto2_phirsch", as.double(x), m, as.integer(n), as.double(k), as.double(s), out=as.double(m), DUP=FALSE, PACKAGE="CITAN")
# 	return(out);
	.C("pareto2_phirsch", out=as.double(x), m, as.double(n), as.double(k), as.double(s), DUP=FALSE, PACKAGE="CITAN")$out
}



#' The probability mass function of Hirsch's \eqn{h}-index
#' for sample of size \code{n} in an i.i.d. Pareto-Type II model.
#'
#' This is a significantly faster (written in C) version of the more general function \code{\link{dhirsch}}.
#'
#' @references
#' Gagolewski M., Grzegorzewski P., S-Statistics and Their Basic Properties, In: Borgelt C. et al (Eds.), Combining Soft Computing and Statistical Methods in Data Analysis, Springer-Verlag, 2010, 281-288.\cr
#'
#' @title Distribution of the h-index - p.m.f.
#' @param x numeric vector.
#' @param n sample size.
#' @param k shape parameter, \eqn{k>0}.
#' @param s scale parameter, \eqn{s>0}.
#' @return The value of the p.m.f. at \code{x}.
#' @export
#' @seealso \code{\link{index.h}}, \code{\link{rho.get}}, \code{\link{dhirsch}}, \code{\link{phirsch}}, \code{\link{pareto2.phirsch}}
pareto2.dhirsch <- function(x, n, k, s)
{
	m <- as.integer(length(x));
# 	.C("pareto2_phirsch", as.double(x), m, as.integer(n), as.double(k), as.double(s), out=as.double(m), DUP=FALSE, PACKAGE="CITAN")
# 	return(out);
	.C("pareto2_dhirsch", out=as.double(x), m, as.double(n), as.double(k), as.double(s), DUP=FALSE, PACKAGE="CITAN")$out
}
