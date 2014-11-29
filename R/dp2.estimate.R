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


#' Finds the maximum likelihood estimator of the Discretized type II Pareto distribution's
#' shape parameter \eqn{k} and scale parameter \eqn{s}.
#'
#' If X has the Pareto-Type II distribution \eqn{P2(k,s)}
#' then \code{Y=floor(X)} has the Discretized Pareto-Type II distribution DP2(k,s).
#'
#' @title Estimation of shape and scale parameters for the Discretized Pareto-Type II distribution (MLE)
#' @param x a non-negative numeric vector.
#' @param kmin lower bound for the shape parameter.
#' @param kmax upper bound for the shape parameter.
#' @param smin lower bound for the scale parameter.
#' @param smax upper bound for the scale parameter.
#' @return
#' The list  with the following components is passed as a result:
#' \tabular{ll}{
#' \code{k} \tab	the estimated parameter of shape.\cr
#' \code{s} \tab	the estimated parameter of scale.\cr
#' }
#' @export
#' @seealso \code{\link{ppareto2}}, \code{\link{discrpareto2.mlekestimate}},\cr
#' \code{\link{discrpareto2.goftest}}
discrpareto2.mleksestimate <- function(x, kmin=1e-4, kmax=100, smin=1e-4, smax=100)
{
	if (mode(x) != "numeric") stop("'x' should be numeric");
	x <- x[!is.na(x)];
	n <- length(x);

	if (n < 2) stop("'x' should be of length at least 2");
	if (any(x < 0.0)) stop("'x' should be non-negative");

	res <- optim(c((kmin+kmax)*0.5,(smin+smax)*0.5), function(p,x,n) {
		-n*p[1]*log(p[2])-sum(log((p[2]+x)^(-p[1])-(1+p[2]+x)^(-p[1])))
	}, gr=NULL, x, n, method="L-BFGS-B", lower=c(kmin,smin), upper=c(kmax,smax))$par

	list(k=res[1], s=res[2]);
}


#' Finds the maximum likelihood estimator of the Discretized type II Pareto distribution's
#' shape parameter \eqn{k} for given scale parameter \eqn{s}.
#'
#' If X has the Pareto-Type II distribution \eqn{P2(k,s)}
#' then \code{Y=floor(X)} has the discretized Pareto-Type II distribution DP2(k,s).
#'
#' @title Estimation of shape parameter for the Discretized Pareto-Type II distribution (MLE)
#' @param x a non-negative numeric vector.
#' @param s scale parameter, \eqn{s>0}.
#' @param kmin lower bound for the shape parameter.
#' @param kmax upper bound for the shape parameter.
#' @return
#' A single numeric value is returned, the ML estimator of \eqn{k}.
#' @export
#' @seealso \code{\link{ppareto2}}, \code{\link{discrpareto2.mleksestimate}},\cr
#' \code{\link{discrpareto2.goftest}}
discrpareto2.mlekestimate <- function(x, s, kmin=1e-4, kmax=100)
{
	if (mode(s) != "numeric" || length(s) != 1 || s <= 0) stop("'s' should be > 0");

	if (mode(x) != "numeric") stop("'x' should be numeric");
	x <- x[!is.na(x)];
	n <- length(x);

	if (n < 2) stop("'x' should be of length at least 2");
	if (any(x < 0.0)) stop("'x' should be non-negative");

	res <- optim((kmin+kmax)*0.5, function(k,x,n,s) {
		-n*k*log(s)-sum(log((s+x)^(-k)-(1+s+x)^(-k)))
	}, gr=NULL, x, n, s, method="L-BFGS-B", lower=c(kmin), upper=c(kmax))$par

	res[1];
}
