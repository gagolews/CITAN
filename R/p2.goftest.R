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


#' Performs goodness-of-fit test for the Pareto-II distribution
#' basing on MLE or MMSE estimates (Zhang, Stevens, 2009) and the Anderson-Darling
#' or Kolmogorov test.
#'
#' This method, proposed e.g. by Zhang and Stevens (2009), uses either the function \code{\link[ADGofTest]{ad.test}} from package \pkg{ADGofTest}
#' or \code{\link{ks.test}} to compute the selected test.
#'
#' It is known that the tests have low powers.
#'
#' If \code{k} and \code{s} are NULL, it bases on bayesian MMS estimators, see \code{\link{pareto2.zsestimate}}.
#' If \code{s} is not NULL, then the unbiased maximum likelihood estimator
#' is used to determine the scale parameter (see \code{\link{pareto2.mlekestimate}}) iff it is not given.
#'
#' @title Goodness-of-fit test for the Pareto-II distribution
#' @param x a non-negative numeric vector of data values.
#' @param k scale parameter, \eqn{k>0} or \code{NULL}.
#' @param s shape parameter, \eqn{s>0} or \code{NULL}.
#' @param method either "anderson-darling" or "kolmogorov".
#'
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
#' @seealso \code{\link{dpareto2}}, \code{\link{pareto2.zsestimate}}, \code{\link{pareto2.mlekestimate}}, \code{\link{ks.test}}, \code{\link[ADGofTest]{ad.test}} from package \code{ADGofTest}
#' @references
#' Zhang J., Stevens M.A., A New and Efficient Estimation Method for the Generalized Pareto Distribution, Technometrics 51(3), 2009, 316-325.\cr
pareto2.goftest <- function(x, k=NULL, s=NULL, #alternative = c("two.sided", "less", "greater"),
	method = c("anderson-darling", "kolmogorov"))
{
# 	alternative <- match.arg(alternative);
	DNAME <- deparse(substitute(x));

	method <- match.arg(method);

	x <- x[!is.na(x)];
	nx <- length(x);
	if (nx < 2L || any(x<0)) stop("incorrect 'x' data");

	if (!is.null(s) && (mode(s) != "numeric" || length(s) != 1 || s <= 0)) stop("'s' should be > 0");
	if (!is.null(k) && (mode(k) != "numeric" || length(k) != 1 || k <= 0)) stop("'k' should be > 0");

	if (is.null(s)) {
		if (!is.null(k)) warning("'k' given but 's' not given. ignoring");
		params <- pareto2.zsestimate(x);
	} else {
		if (is.null(k)) k <- pareto2.mlekestimate(x, s);
		params <- list(k=k, s=s);
	}

	stopifnot(params$k > 0 && is.finite(params$k));
	stopifnot(params$s > 0 && is.finite(params$s));

	RVAL <- switch(method,
		"anderson-darling" = ad.test(x, ppareto2, params$k, params$s),
		"kolmogorov" = ks.test(x, "ppareto2", params$k, params$s)
	);

	RVAL$method = switch(method,
		"anderson-darling" = sprintf("Anderson-Darling Goodness-of-Fit test for the Pareto-II distribution P2(%g, %g)", params$k, params$s),
		"kolmogorov" = sprintf("Kolmogorov Goodness-of-Fit test for the Pareto-II distribution P2(%g, %g)", params$k, params$s),
	);
	
	RVAL$data.name <- DNAME;

	return(RVAL);
}






