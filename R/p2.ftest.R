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


#' Performs F-test for equality of shape parameters
#' of two samples from the Pareto type-II distributions with known
#' and equal scale parameters, \eqn{s>0}.
#'
#' Given two samples \eqn{(X_1,...,X_n)} i.i.d. \eqn{P2(k_x,s)}
#' and \eqn{(Y_1,...,Y_m)} i.i.d. \eqn{P2(k_y,s)}
#' this test verifies the null hypothesis
#' \eqn{H_0: k_x=k_y}
#' against two-sided or one-sided alternatives, depending
#' on the value of \code{alternative}.
#' It bases on test statistic
#' \code{T=n/m*sum(log(1+Y/m))/sum(log(1+X/n))}
#' which, under \eqn{H_0}, has the Snedecor's F distribution with \eqn{(2m, 2n)}
#' degrees of freedom.
#'
#' Note that for \eqn{k_x < k_y}, then \eqn{X} dominates \eqn{Y} stochastically.
#'
#' @title Two-sample F-test for equality of shape parameters for Type II-Pareto distributions with known common scale parameter
#' @param x a non-negative numeric vector of data values.
#' @param y a non-negative numeric vector of data values.
#' @param s scale parameter, \eqn{s>0}.
#' @param alternative indicates the alternative hypothesis and must be one of "two.sided" (default), "less", or "greater".
#' @param significance significance level, \eqn{0<}\code{significance}\eqn{<1} or \code{NULL}. See Value for details.
#' @return
#' If \code{significance} is not \code{NULL}, then
#' the list of class \code{power.htest} with the following components is passed as a result:
#' \tabular{ll}{
#' \code{statistic} \tab	the value of the test statistic.\cr
#' \code{result} \tab	either FALSE (accept null hypothesis) or TRUE (reject).\cr
#' \code{alternative} \tab	a character string describing the alternative hypothesis.\cr
#' \code{method} \tab	a character string indicating what type of test was performed.\cr
#' \code{data.name} \tab	a character string giving the name(s) of the data.\cr
#' }
#' Otherwise, the list of class \code{htest} with the following components is passed as a result:
#' \tabular{ll}{
#' \code{statistic} \tab	the value of the test statistic.\cr
#' \code{p.value} \tab	the p-value of the test.\cr
#' \code{alternative} \tab	a character string describing the alternative hypothesis.\cr
#' \code{method} \tab	a character string indicating what type of test was performed.\cr
#' \code{data.name} \tab	a character string giving the name(s) of the data.\cr
#' }
#' @export
#' @seealso \code{\link{dpareto2}}, \code{\link{pareto2.goftest}}, \code{\link{pareto2.htest}}, \code{\link{pareto2.htest.approx}}
pareto2.ftest <- function(x, y, s, alternative = c("two.sided", "less", "greater"), significance=NULL)
{
	alternative <- match.arg(alternative);
	DNAME <- deparse(substitute(x));
	DNAME <- paste(DNAME, "and", deparse(substitute(y)));

	x <- x[!is.na(x)];
	nx <- length(x);
	if (nx < 1L || any(x<0)) stop("incorrect 'x' data");

	y <- y[!is.na(y)];
	ny <- length(y);
	if (ny < 1L || any(y<0)) stop("incorrect 'y' data");


	if (mode(s) != "numeric" || length(s) != 1 || s <= 0) stop("'s' should be > 0");




	STATISTIC <- nx/ny*sum(log(1+y/s))/sum(log(1+x/s));
	names(STATISTIC) <- "F";

	METHOD <- "Two-sample F-test for equality of shape parameters for Type II-Pareto distributions with known common scale parameter";
	nm_alternative <- switch(alternative, two.sided = "two-sided",
		less = "kx < ky",
		greater = "kx > ky");



	if (!is.null(significance))
	{
		if (length(significance) != 1 || significance <= 0 || significance >=1)
			stop("incorrect 'significance'");

		if (significance > 0.2) warning("'significance' is possibly incorrect");


		RESULT <- ifelse(alternative == "two.sided", (STATISTIC<qf(significance*0.5,2*ny,2*nx) || STATISTIC>qf(1-significance*0.5,2*ny,2*nx)),
		          ifelse(alternative == "greater",    STATISTIC>qf(1-significance,2*ny,2*nx),
		                                              STATISTIC<qf(significance,2*ny,2*nx)));

		RVAL <- list(statistic = STATISTIC, result = RESULT, alternative = nm_alternative,
			method = METHOD, data.name = DNAME);
		class(RVAL) <- "power.htest";
		return(RVAL);
	} else {

		PVAL <- ifelse(alternative == "two.sided", (0.5-abs(pf(STATISTIC, 2*ny, 2*nx)-0.5))*2,
		        ifelse(alternative == "greater",   1-pf(STATISTIC, 2*ny, 2*nx),
		                                             pf(STATISTIC, 2*ny, 2*nx)));


		RVAL <- list(statistic = STATISTIC, p.value = PVAL, alternative = nm_alternative,
			method = METHOD, data.name = DNAME);
		class(RVAL) <- "htest";
		return(RVAL);
	}
}

