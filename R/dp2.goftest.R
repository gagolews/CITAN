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


#' Performs goodness-of-fit test for the Discretized Pareto-II distribution
#' basing on MLE estimates and the chi-square test.
#'
#' If X has the Pareto-Type II distribution \eqn{P2(k,s)}
#' then \code{Y=floor(X)} has the Discretized Pareto-Type II distribution DP2(k,s).
#'
#' It is known that the test has low power.
#'
#' @title Goodness-of-fit test for the Discretized Pareto-II distribution
#' @param x a non-negative numeric vector of data values of length >= 9.
#' @param k scale parameter, \eqn{k>0} or \code{NULL}.
#' @param s shape parameter, \eqn{s>0} or \code{NULL}.
#' @param kmin lower bound for the shape parameter.
#' @param kmax upper bound for the shape parameter.
#' @param smin lower bound for the scale parameter.
#' @param smax upper bound for the scale parameter.
#' @return
#' The list of class \code{htest} with the following components is passed as a result:
#' \tabular{ll}{
#' \code{statistic} \tab	the value of the test statistic.\cr
#' \code{p.value} \tab	the p-value of the test.\cr
#' \code{alternative} \tab	a character string describing the alternative hypothesis.\cr
#' \code{method} \tab	a character string indicating what type of test was performed.\cr
#' \code{data.name} \tab	a character string giving the name(s) of the data.\cr
#' }
#' @export
#' @seealso \code{\link{discrpareto2.mleksestimate}}, \code{\link{discrpareto2.mlekestimate}},
#' \code{\link{chisq.test}}
discrpareto2.goftest <- function(x, k=NULL, s=NULL, kmin=1e-4, kmax=100, smin=1e-4, smax=100)
{
	DNAME <- deparse(substitute(x));

	x <- x[!is.na(x)];
	nx <- length(x);
	if (nx < 2L || any(x<0)) stop("incorrect 'x' data");

	if (!is.null(s) && (mode(s) != "numeric" || length(s) != 1 || s <= 0)) stop("'s' should be > 0");
	if (!is.null(k) && (mode(k) != "numeric" || length(k) != 1 || k <= 0)) stop("'k' should be > 0");

	if (is.null(s)) {
		if (!is.null(k)) warning("'k' given but 's' not given. ignoring");
		par <- discrpareto2.mleksestimate(x, kmin=kmin, kmax=kmax, smin=smin, smax=smax);
	} else {
		if (is.null(k)) k <- discrpareto2.mlekestimate(x, s, kmin=kmin, kmax=kmax);
		par <- list(k=k, s=s);
	}

	stopifnot(par$k > 0 && is.finite(par$k));
	stopifnot(par$s > 0 && is.finite(par$s));

	## -------------------------------------------------------------------

	sn <- min(floor(sqrt(nx)), ceiling(1/ppareto2(1, par$k, par$s)));

	if (sn < 3) stop("could not create at least 3 classes.");
	
	q <- qpareto2(seq(0, 1, length=sn+1), par$k, par$s);
	q <- ceiling(q);
	
	if (length(unique(q))!=length(q))
		stop("could not create distinct classes for given 'x'.");
	
	p <- ppareto2(q[-1], par$k, par$s)-ppareto2(q[-(sn+1)], par$k, par$s);
	
	x <- table(cut(x, breaks=q, right=F));
	
	RVAL <- chisq.test(x, p=p, rescale.p=TRUE);

	RVAL$method <-  sprintf("Chi-square goodness-of-fit test for the discretized Pareto-II distribution DP2(%g, %g)", par$k, par$s);
	RVAL$data.name <- DNAME;
	
	return(RVAL);
}






