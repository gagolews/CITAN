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


#' Computes the density function of the S-statistic w.r.t. to a control function in an i.i.d. model with common increasing and continuous c.d.f. \eqn{F} defined on \eqn{[0,\infty)}.
#'
#' The function computes the value of the p.d.f. of an S-statistic
#' w.r.t. to the control function \code{kappa} for sample of size \code{n}.
#' Note that the result is valid
#' only at continuity points of \eqn{F'}.
#'
#' For more information see man page on  \code{\link{psstat}} and the paper (Gagolewski, Grzegorzewski, 2010).
#'
#' @references
#' Gagolewski M., Grzegorzewski P., S-Statistics and Their Basic Properties, In: Borgelt C. et al (Eds.), Combining Soft Computing and Statistical Methods in Data Analysis, Springer-Verlag, 2010, 281-288.\cr
#'
#' @title Distribution of S-statistics - density
#' @param x numeric vector.
#' @param n sample size.
#' @param cdf a cumulative distribution function \eqn{F}, e.g. \code{\link{ppareto2}}.
#' @param pdf a density function \eqn{F'}, e.g. \code{\link{dpareto2}}.
#' @param kappa an increasing function, \eqn{\kappa}, a so-called control function.
#' @param kappaInvDer the derivative of the inverse of \eqn{\kappa}.
#' @param ... optional arguments to \code{cdf} and \code{pdf}.
#' @return The value of the density at \code{x}.
#' @export
#' @seealso \code{\link{Sstat}}, \code{\link{Sstat2}}, \code{\link{psstat}}, \code{\link{rho.get}}
dsstat <- function(x, n, cdf, pdf, kappa, kappaInvDer, ...)
{
	(cdf(kappa(x), ...)^(n-floor(x*n)-1))*
	((1-cdf(kappa(x), ...))^floor(x*n))*
	(pdf(kappa(x), ...)/abs(kappaInvDer(kappa(x))))/
	(beta(n-floor(x*n),floor(x*n)+1))
}


#' Computes the cumulative distribution function of the S-statistic w.r.t. to a control function in an i.i.d. model with common increasing and continuous c.d.f. \eqn{F} defined on \eqn{[0,\infty)}.
#'
#' Let \eqn{F} (parameter \code{cdf}) be a continuous c.d.f that is strictly increasing on \eqn{[a,b]},
#' where \eqn{a=\inf\{x: F(x)>0\}}{a=inf{x: F(x)>0}} and
#' \eqn{b=\sup\{x: F(x)<1\}}{b=sup{x: F(x)<1}}.
#'
#' Moreover, let
#' \eqn{\kappa:[0,1]\to[c,d]\subseteq[a,b]}{\kappa:[0,1]->[c,d]c[a,b]}
#' be any control function (parameter \code{kappa}), i.e. a function that
#' is continuous and strictly increasing
#' and which fulfills \eqn{\kappa(0)=c} and \eqn{\kappa(1)=d}.
#'
#' The function computes the value of the c.d.f. of an S-statistic
#' w.r.t. to the control function for sample of size \code{n}.
#' This result was given in (Gagolewski, Grzegorzewski, 2010).
#'
#' Note that under certain conditions the distribution of an S-statistic
#' is asymptotically normal with expectation \eqn{\rho_\kappa} (see \code{\link{rho.get}})
#' and variance \eqn{\rho_\kappa (1-\rho_\kappa)/n/(1+F(\kappa(\rho_\kappa)))^2}.
#'
#' @references
#' Gagolewski M., Grzegorzewski P., S-Statistics and Their Basic Properties, In: Borgelt C. et al (Eds.), Combining Soft Computing and Statistical Methods in Data Analysis, Springer-Verlag, 2010, 281-288.\cr
#'
#' @title Distribution of S-statistics - c.d.f.
#' @param x numeric vector.
#' @param n sample size.
#' @param cdf a cumulative distribution function \eqn{F}, e.g. \code{\link{ppareto2}}.
#' @param kappa an increasing function, \eqn{\kappa} (see Details), a so-called control function.
#' @param ... optional arguments to \code{cdf} and \code{pdf}.
#' @return The value of the c.d.f. at \code{x}.
#' @export
#' @seealso \code{\link{Sstat}}, \code{\link{Sstat2}}, \code{\link{dsstat}}, \code{\link{rho.get}}
psstat <- function(x, n, cdf, kappa, ...)
{
	warn <- getOption("warn");
	options("warn"=-1);
	y <- ifelse(x>1-1e-16, 1.0,
	     ifelse(x<  1e-16, 0.0,
	                       pbeta(cdf(kappa(x), ...), n-floor(n*x),floor(n*x)+1)));
	options("warn"=warn);
	return(y);
}
