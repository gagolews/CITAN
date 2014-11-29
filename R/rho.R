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


#' Numerically computes the rho-index of a given continuous cumulative distribution function
#'
#' Let \eqn{F} be a continuous c.d.f that is strictly increasing on \eqn{[a,b]},
#' where \eqn{a=\inf\{x: F(x)>0\}}{a=inf{x: F(x)>0}} and
#' \eqn{b=\sup\{x: F(x)<1\}}{b=sup{x: F(x)<1}}.
#'
#' A \dfn{control function} is any function
#' \eqn{\kappa:[0,1]\to[c,d]\subseteq[a,b]}{\kappa:[0,1]->[c,d]c[a,b]} that
#' is continuous and strictly increasing
#' and which fulfills \eqn{\kappa(0)=c} and \eqn{\kappa(1)=d}.
#'
#' The \dfn{\eqn{\rho}-index} of the distribution \eqn{F} (Gagolewski, Grzegorzewski, 2010)
#' is a number \eqn{\rho_\kappa\in(0,1)}{0<\rho_\kappa<1}
#' such that
#' \deqn{\rho_\kappa=1-F(\kappa(\rho_\kappa)).}
#'
#' It turns out that under certain conditions in a model of i.i.d. random variables
#' the S-statistic associated with \eqn{\kappa} is an asymptotically
#' unbiased, normal and strongly consistent estimator of \eqn{\rho_\kappa}.
#'
#'
#' @references
#' Gagolewski M., Grzegorzewski P., S-Statistics and Their Basic Properties, In: Borgelt C. et al (Eds.), Combining Soft Computing and Statistical Methods in Data Analysis, Springer-Verlag, 2010, 281-288.\cr
#'
#' @title Rho-index of a continuous probability distribution
#' @param cdf a cumulative distribution function, e.g. \code{\link{ppareto2}}.
#' @param kappa an increasing function, \eqn{\kappa} (see Details), a so-called control function.
#' @param ... optional arguments to \code{cdf}.
#' @param tol the desired accuracy (convergence tolerance).
#' @return The function returns a single number.
#' @export
#' @seealso \code{\link{phirsch}}, \code{\link{dhirsch}}, \code{\link{psstat}}, \code{\link{dsstat}}, \code{\link{Sstat}}, \code{\link{Sstat2}}
#' @examples
#' kappa <- function(x) { pmax(0,pmin(x,1)) } # identity function on [0,1]
#' rho.get(ppareto2, kappa, 1, 1)             # golden ratio
rho.get <- function(cdf, kappa, ..., tol=1e-20)
{
	uniroot(function(x) { 1-cdf(kappa(x), ...)-x; },
		c(0,1), tol=tol)$root;
}
