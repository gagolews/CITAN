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


#' Computes the exact right-sided confidence interval for the theoretical \eqn{h}-index
#' of a probability distribution in an \eqn{(X_1,\dots,X_n)} i.i.d. Pareto-type II
#' model with known scale parameter \eqn{s>0}.
#'
#' See \code{\link{pareto2.confint.h}} for details.
#'
#' @references
#' Gagolewski M., Grzegorzewski P., S-Statistics and Their Basic Properties, In: Borgelt C. et al (Eds.),
#' Combining Soft Computing and Statistical Methods in Data Analysis, Springer-Verlag, 2010, 281-288.\cr
#'
#' @title Right-sided exact confidence interval for the theoretical h-index
#' @param h observed value of the \eqn{h}-index
#' @param s scale parameter, \eqn{s>0}.
#' @param n sample size.
#' @param conf.level confidence level; defaults 0.95.
#' @param tol the desired accuracy (convergence tolerance).
#' @return Upper bound of the confidence interval.
#' @export
#' @seealso \code{\link{index.h}}, \code{\link{ppareto2}}, \code{\link{rho.get}},
#' \code{\link{pareto2.confint.rho}},\cr
#' \code{\link{pareto2.confint.h}}, \code{\link{pareto2.confint.h.lower}}
pareto2.confint.h.upper <- function(h, s, n, conf.level=0.95, tol=1e-12)
{
	gamma <- 1-conf.level;

	if (mode(s) != "numeric" || length(s) != 1 || s <= 0) stop("'s' should be > 0");

	if (length(h) != 1 || h < 0 || h > n)
		stop("Incorrect h value!");
	h <- round(h);

	if (h > n-1e-9) return(n);
	if (gamma < 1e-9) return(n);


	xsol <- uniroot(function(x,h,s,n,gamma) {
		pareto2.phirsch(h+1e-9,n,x,s)-gamma;
	}, c(0,1e25),h,s,n,gamma, tol=tol, maxiter=1000)$root;

	n*rho.get(ppareto2, function(x) { pmin(1,pmax(0,x))*n }, xsol, s); # return value
}




#' Computes the exact left-sided confidence interval for the theoretical \eqn{h}-index
#' of a probability distribution in an \eqn{(X_1,\dots,X_n)} i.i.d. Pareto-type II
#' model with known scale parameter \eqn{s>0}.
#'
#' See \code{\link{pareto2.confint.h}} for details.
#'
#' @references
#' Gagolewski M., Grzegorzewski P., S-Statistics and Their Basic Properties, In: Borgelt C. et al (Eds.),
#' Combining Soft Computing and Statistical Methods in Data Analysis, Springer-Verlag, 2010, 281-288.\cr
#'
#' @title Left-sided exact confidence interval for the theoretical h-index
#' @param h observed value of the \eqn{h}-index
#' @param s scale parameter, \eqn{s>0}.
#' @param n sample size.
#' @param conf.level confidence level; defaults 0.95.
#' @param tol the desired accuracy (convergence tolerance).
#' @return Lower bound of the confidence interval.
#' @export
#' @seealso \code{\link{index.h}}, \code{\link{ppareto2}}, \code{\link{rho.get}},
#' \code{\link{pareto2.confint.rho}},\cr
#' \code{\link{pareto2.confint.h}}, \code{\link{pareto2.confint.h.lower}}
pareto2.confint.h.lower <- function(h, s, n, conf.level=0.95, tol=1e-12)
{
	gamma <- 1-conf.level;

	if (mode(s) != "numeric" || length(s) != 1 || s <= 0) stop("'s' should be > 0");

	if (length(h) != 1 || h < 0 || h > n)
		stop("Incorrect h value!");
	h <- round(h);

	if (h < 1e-9) return(0);
	if (gamma < 1e-9) return(0);


	xsol <- uniroot(function(x,h,s,n,gamma) {
		pareto2.phirsch(h-1+1e-9,n,x,s)-1+gamma;
	}, c(0,1e25),h,s,n,gamma, tol=tol, maxiter=1000)$root;

	n*rho.get(ppareto2, function(x) { pmin(1,pmax(0,x))*n }, xsol, s); # return value
}



#' Computes the exact two-sided confidence interval for the theoretical \eqn{h}-index
#' of a probability distribution in an \eqn{(X_1,\dots,X_n)} i.i.d. Pareto-type II
#' model with known scale parameter \eqn{s>0}.
#'
#' The \dfn{Theoretical \eqn{h}-index} for a sequence of \eqn{n} i.i.d. random variables
#' with common increasing and continuous c.d.f. \eqn{F} defined on \eqn{[0,\infty)}
#' is equal to \eqn{n\varrho_\kappa}{n*\rho_\kappa}, where \eqn{\rho_\kappa}
#' is the \eqn{\rho}-index of \eqn{F} for \eqn{\kappa(x)=nx}, see \code{\link{rho.get}} for details.
#'
#' @references
#' Gagolewski M., Grzegorzewski P., S-Statistics and Their Basic Properties, In: Borgelt C. et al (Eds.),
#' Combining Soft Computing and Statistical Methods in Data Analysis, Springer-Verlag, 2010, 281-288.\cr
#'
#' @title Two-sided exact confidence interval for the theoretical h-index
#' @param h observed value of the \eqn{h}-index
#' @param s scale parameter, \eqn{s>0}.
#' @param n sample size.
#' @param conf.level confidence level; defaults 0.95.
#' @param tol the desired accuracy (convergence tolerance).
#' @return Vector of length 2 with the computed bounds of the confidence interval.
#' @export
#' @seealso \code{\link{index.h}}, \code{\link{ppareto2}}, \code{\link{rho.get}},
#' \code{\link{pareto2.confint.rho}},\cr
#' \code{\link{pareto2.confint.h.upper}}, \code{\link{pareto2.confint.h.upper}}
pareto2.confint.h <- function(h, s, n, conf.level=0.95, tol=1e-12)
{
	gamma <- 1-conf.level;
	return(c(
		pareto2.confint.h.lower(h,s,n,1-gamma*0.5,tol),
		pareto2.confint.h.upper(h,s,n,1-gamma*0.5,tol)
	));
}
