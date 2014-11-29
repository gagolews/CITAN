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


#' Finds the MMS estimator of the type II Pareto distribution parameters
#' using the Bayesian method (and the R code) developed by
#' Zhang and Stevens (2009).
#'
#' @title Estimation of parameters for the Pareto-II distribution (MMSE)
#' @param x a non-negative numeric vector.
#' @return
#' The list  with the following components is passed as a result:
#' \tabular{ll}{
#' \code{k} \tab	the estimated parameter of shape.\cr
#' \code{s} \tab	the estimated parameter of scale.\cr
#' }
#' @export
#' @seealso \code{\link{dpareto2}}, \code{\link{pareto2.goftest}}, \code{\link{pareto2.mlekestimate}}, \code{\link{pareto2.mleksestimate}}
#' @references
#' Zhang J., Stevens M.A., A New and Efficient Estimation Method for the Generalized Pareto Distribution, Technometrics 51(3), 2009, 316-325.\cr
pareto2.zsestimate <- function(x)
{
	if (mode(x) != "numeric") stop("'x' should be numeric");
	x <- x[!is.na(x)];
	n <- length(x);
	x <- sort(x);

	if (n < 2) stop("'x' should be of length at least 2");
	if (x[n] <= 0.0) stop("'x' should be non-negative");


	lx <- function(b, x)
	{
		k <- -mean(log(1-b*x));
		log(b/k)+k-1;
	}

	m <- 20+floor(sqrt(n))

	b <- w <- L <- 1/x[n]+(1-sqrt(m/((1:m)-0.5)))/3/x[floor(n/4+0.5)]

	for (i in 1:m) L[i] <- n*lx(b[i],x)

	for (i in 1:m) w[i]<- 1/sum(exp(L-L[i]))

	b <- sum(b*w);

	k <- 1/mean(log(1-b*x));
	s <- -1/b;

	if (k<=0)  warning("estimated shape parameter <= 0");
	if (s<=0)  warning("estimated scale parameter <= 0");

	list(k=k, s=s);
}


#' Finds the unbiased maximum likelihood estimator of the type II Pareto distribution's
#' shape parameter \eqn{k} for known scale parameter \eqn{s}.
#'
#' @title Estimation of shape parameter for the Pareto-II distribution (MLE)
#' @param x a non-negative numeric vector.
#' @param s scale parameter, \eqn{s>0}.
#' @return
#' A single numeric value is returned, the unbiased ML estimator of \eqn{k}.
#' @export
#' @seealso \code{\link{dpareto2}}, \code{\link{pareto2.goftest}}, \code{\link{pareto2.zsestimate}}, \code{\link{pareto2.mleksestimate}}
pareto2.mlekestimate <- function(x, s)
{
	if (mode(s) != "numeric" || length(s) != 1 || s <= 0) stop("'s' should be > 0");

	if (mode(x) != "numeric") stop("'x' should be numeric");
	x <- x[!is.na(x)];
	n <- length(x);

	if (n < 2) stop("'x' should be of length at least 2");
	if (any(x < 0.0)) stop("'x' should be non-negative");

	return((n-1)/sum(log(1+x/s)));
}



#' Finds the maximum likelihood estimator of the type II Pareto distribution's
#' shape parameter \eqn{k} and scale parameter \eqn{s}.
#'
#' Note that the maximum of the likelihood function may not exist
#' for some input vectors.
#'
#' @title Estimation of shape and scale parameters for the Pareto-II distribution (MLE)
#' @param x a non-negative numeric vector.
#' @param tol the desired accuracy (convergence tolerance).
#' @param smin lower bound for the scale parameter.
#' @param smax upper bound for the scale parameter.
#' @return
#' The list  with the following components is passed as a result:
#' \tabular{ll}{
#' \code{k} \tab	the estimated parameter of shape.\cr
#' \code{s} \tab	the estimated parameter of scale.\cr
#' }
#' or \code{NA} if the maximum of the likelihood function does not exist.
#' @export
#' @seealso \code{\link{dpareto2}}, \code{\link{pareto2.goftest}}, \code{\link{pareto2.zsestimate}}
pareto2.mleksestimate <- function(x, tol=1e-20, smin=1e-4, smax=20)
{
	if (mode(x) != "numeric") stop("'x' should be numeric");
	x <- x[!is.na(x)];
	n <- length(x);

	if (n < 2) stop("'x' should be of length at least 2");
	if (any(x < 0.0)) stop("'x' should be non-negative");


	flow <- 1+sum(log(1+x/smin))/n-n/sum(1/(1+x/smin));
	fupp <- 1+sum(log(1+x/smax))/n-n/sum(1/(1+x/smax));

	if (flow*fupp >= 0)
	{
		warning("Maximum of the likelihood function does not exist.");
		return(list(k=NA, s=NA));
	}

	s <- uniroot( function(s,x) {
		dx <- 1+x/s;
		1+sum(log(dx))/n-n/sum(1/dx);
	}, c(smin, smax), x, f.lower=flow, f.upper=fupp, tol=tol)$root;
	k <- (n/sum(log(1+x/s)));

	list(k=k, s=s);
}

