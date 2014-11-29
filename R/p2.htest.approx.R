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


#' Performs asymptotic (approximate) \eqn{h}-test for equality of shape parameters
#' of two samples from the Pareto type-II distributions with known
#' and equal scale parameters, \eqn{s>0}.
#'
#' Given two equal-sized samples \eqn{X=(X_1,...,X_n)} i.i.d. \eqn{P2(k_x,s)}
#' and \eqn{Y=(Y_1,...,Y_m)} i.i.d. \eqn{P2(k_y,s)}
#' this test verifies the null hypothesis
#' \eqn{H_0: k_x=k_y}
#' against two-sided or one-sided alternatives, depending
#' on the value of \code{alternative}.
#' It bases on a test statistic that is a function of {H(Y)-H(X)},
#' where \eqn{H} denotes Hirsch's \eqn{h}-index (see \code{\link{index.h}}).
#' This statistic approximately has asymptotically standardized normal distribution under \eqn{H_0}.
#'
#' Note that for \eqn{k_x < k_y}, then \eqn{X} dominates \eqn{Y} stochastically.
#'
#' @title Two-sample asymptotic h-test for equality of shape parameters for Type II-Pareto distributions with known common scale parameter
#' @param x an n-element non-negative numeric vector of data values.
#' @param y an n-element non-negative numeric vector of data values.
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
#' @seealso \code{\link{dpareto2}}, \code{\link{pareto2.goftest}}, \code{\link{pareto2.ftest}}, \code{\link{pareto2.htest}}, \code{\link{index.h}}
#' @references
#' Gagolewski M., Grzegorzewski P., S-Statistics and Their Basic Properties, In: Borgelt C. et al (Eds.),
#' Combining Soft Computing and Statistical Methods in Data Analysis, Springer-Verlag, 2010, 281-288.\cr
pareto2.htest.approx <- function(x, y, s, alternative = c("two.sided", "less", "greater"), significance=NULL)
{
	alternative <- match.arg(alternative);
	DNAME <- deparse(substitute(x));
	DNAME <- paste(DNAME, "and", deparse(substitute(y)));

	if (mode(s) != "numeric" || length(s) != 1 || s <= 0) stop("'s' should be > 0");

	if (mode(x) != "numeric" || mode(y) != "numeric") stop("non-numeric data given");

	x <- x[!is.na(x)];
	n <- length(x);
	if (n < 1L || any(x<0)) stop("incorrect 'x' data");

	y <- y[!is.na(y)];
	if (length(y) < 1L || any(y<0)) stop("incorrect 'y' data");

	if (length(y) != n) stop("non-equal-sized vectors given on input");



	METHOD <- "Two-sample asymptotic h-test for equality of shape parameters for Type II-Pareto distributions with known common scale parameter";

	nm_alternative <- switch(alternative, two.sided = "two-sided",
			less = "kx < ky",
			greater = "kx > ky");


	# -----------------------------------------------------------------------

	HYn <- index.h(y,disable.check=TRUE)/n;
	HXn <- index.h(x,disable.check=TRUE)/n;

	kappa    <- function(x) { pmax(0,pmin(1,x))*n; }

# 	if (HXn < 1e-9)
# 	{
# 		gprimex <- 0.0;
# 	} else {
		kx <- uniroot(function(k,s,targetrho,kappa)
				{
					1-ppareto2(kappa(targetrho),k,s)-targetrho;
				}, c(1e-15,1e15), s, HXn, kappa, tol=1e-20)$root;
		gprimex <- dpareto2(HXn*n,kx,s)*n;
# 	}

# 	if (HYn < 1e-9)
# 	{
# 		gprimey <- 0.0;
# 	} else {
		ky <- uniroot(function(k,s,targetrho,kappa)
				{
					1-ppareto2(kappa(targetrho),k,s)-targetrho;
				}, c(1e-15,1e15), s, HYn, kappa, tol=1e-20)$root;
		gprimey <- dpareto2(HYn*n,ky,s)*n;
# 	}

	sigmax2 <- HXn*(1-HXn)/n/(1+gprimex)^2;
	sigmay2 <- HYn*(1-HYn)/n/(1+gprimey)^2;


	STATISTIC <- (HYn-HXn)/sqrt(sigmax2+sigmay2);

	names(STATISTIC) <- "T";

	# -----------------------------------------------------------------------

	if (!is.null(significance))
	{
		if (length(significance) != 1 || significance <= 0 || significance >=1)
			stop("incorrect 'significance'");

		if (significance > 0.2) warning("'significance' is possibly incorrect");

		RESULT <- ifelse(alternative == "two.sided", (STATISTIC<qnorm(significance*0.5) || STATISTIC>qnorm(1-significance*0.5)),
		          ifelse(alternative == "greater",    STATISTIC>qnorm(1-significance),
		                                              STATISTIC<qnorm(significance)));

		RVAL <- list(statistic = STATISTIC, result = RESULT, alternative = nm_alternative,
			method = METHOD, data.name = DNAME);
		class(RVAL) <- "power.htest";
		return(RVAL);
	} else {

		PVAL <- ifelse(alternative == "two.sided", (0.5-abs(pnorm(STATISTIC)-0.5))*2,
		        ifelse(alternative == "greater",   1-pnorm(STATISTIC),
		                                             pnorm(STATISTIC)));


		RVAL <- list(statistic = STATISTIC, p.value = PVAL, alternative = nm_alternative,
			method = METHOD, data.name = DNAME);
		class(RVAL) <- "htest";
		return(RVAL);
	}
}
