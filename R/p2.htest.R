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


# #' /internal-generalized version/
# .htest.getpowerupper <- function(n, wyg, cdf, PARAM_X, PARAM_Y, ...)
# {
# 	x <- (0:n);
#
# 	stopifnot(length(wyg)==n+1);
#
# 	POW <- numeric(length(PARAM_X));
# 	for (i in 1:length(PARAM_X))
# 	{
# 		POW[i] <- 1-sum(dhirsch(x,n,cdf,PARAM_X[i], ...)*phirsch(x+wyg+0.25,n,cdf,PARAM_Y[i], ...));
#
# 	}
#
# 	return(POW);
# }

# #' /internal-generalized version/
# .htest.acceptreg.improve <- function(n, cdf, wyg, j0, j1, alpha, verbose, PARAM, ...)
# {
# 	epsbound <- 1e-7;
# 	wn <- n+1;
#
# 	for (i in j0:j1)
# 	{
# 		lbound <- ifelse(i>1 && i < wn, min(wyg[i-1], wyg[i+1]), 0)
# 		if (wyg[i] > lbound)
# 		{
# 			repeat
# 			{
# 				wyg[i] <- wyg[i]-1L;
# 				POW <- .htest.getpowerupper(n, wyg, cdf, PARAM, PARAM, ...);
# 				if (max(POW)>alpha)
# 				{
# 					wyg[i] <- wyg[i]+1L;
# 					break;
# 				}
# 				if (wyg[i] == lbound) break;
# 			}
#
# 			if (verbose)
# 			{
# 				POW <- .htest.getpowerupper(n, wyg, cdf, PARAM, PARAM, ...);
# 				cat(sprintf("%3.0f%% complete in current iteration, q=%.4f.\r", abs((i-j0)/(j1-j0))*100, mean(POW)));
# 			}
# 		}
# 	}
#
# 	if (verbose)
# 	{
# 		POW <- .htest.getpowerupper(n, wyg, cdf, PARAM, PARAM, ...);
# 		cat(sprintf("%3.0f%% complete in current iteration, q=%.4f.\r", 100, mean(POW)));
# 	}
#
# 	return(wyg);
# }



#' /internal/
.pareto2.htest.getpowerupper <- function(n, wyg, K_X, K_Y, s)
{
	x <- (0:n);

	stopifnot(length(wyg)==n+1);

	POW <- numeric(length(K_X));
	for (i in 1:length(K_X))
	{
		POW[i] <- 1-sum(pareto2.dhirsch(x,n,K_X[i],s)*pareto2.phirsch(x+wyg+0.25,n,K_Y[i],s));
	}

	return(POW);
}






#' /internal/
.pareto2.htest.getsize_optimized <- function(n, x, wyg, K, s, dhirp2mat)
{
	POW <- numeric(length(K));
	for (i in 1:length(K))
	{
		POW[i] <- 1-sum(dhirp2mat[i,]*pareto2.phirsch(x+wyg+0.25,n,K[i], s));
	}

	return(POW);
}


#' /internal/
.pareto2.htest.acceptreg.improve2 <- function(n, wyg, alpha, verbose, K, s)
{
	epsbound <- 1e-7;
	wn <- n+1;


	# prepare input data for .pareto2.htest.getsize_optimized
	# this significanlty speeds up computations
	xhir0n <- (0:n);
	dhirp2mat <- matrix(nrow=length(K), ncol=length(xhir0n));
	for (i in 1:length(K))
		dhirp2mat[i,] <- pareto2.dhirsch(xhir0n,n,K[i],s);
	# -------------------------------------------------------

	stopifnot(length(wyg)==wn);

	v <- max(wyg);
	wyg <- rep(v,wn);
	stopifnot(max(.pareto2.htest.getsize_optimized(n, xhir0n, wyg, K, s, dhirp2mat)) <= alpha);
	stopifnot(max(.pareto2.htest.getsize_optimized(n, xhir0n, wyg-1, K, s, dhirp2mat)) > alpha);

	k1 <- wn-1;
	while (max(.pareto2.htest.getsize_optimized(n, xhir0n, c(rep(v-1,k1), rep(v,wn-k1)), K, s, dhirp2mat))>=alpha)
		k1 <- k1-1;

	k1 <- k1+1;
	k0 <- k1-1;
	maxdiff <- 2;
# 	while (max(.pareto2.htest.getsize_optimized(n, xhir0n, c(rep(v,k0), rep(v-1,wn-k0)), K, s, dhirp2mat))>=alpha)
# 		k0 <- k0+1;

	cat(sprintf("wn=%g; k0=%g; k1=%g; maxdiff=%g\n", wn, k0, k1, maxdiff));
	stopifnot(k0<k1);


	sizemax <- 0;
	qualmax <- 10;
	wygmax <- wyg;

	imprrec <- function(wyg, k0, k1)
	{
		if (k0 <= 0 && k1 > wn) return();

# 		if (k0 <= 0) k0 <- 1;
# 		if (k1 > wn) k1 <- wn;
		v0 <- ifelse(k0<=0, wyg[1], wyg[k0]);
		v1 <- ifelse(k1>wn, wyg[wn], wyg[k1]);

		j0 <- 0;
		while(T)
		{
			if (k0 >= 1) wyg [1:k0] <- v0-j0;

			j1 <- 0;
			while (T)
			{
				if (k1 <= wn) wyg[k1:wn] <- v1-j1;

				POW <- .pareto2.htest.getsize_optimized(n, xhir0n, wyg, K, s, dhirp2mat);
				mp <- max(POW);

				if (mp <= alpha)
				{
					mep <- mean((POW-alpha)^2);

					if (mep<=qualmax)
					{
						wygmax <<- wyg;
						qualmax <<- mep;
						sizemax <<- mp;
						cat(sprintf("size=%.5f and q=%.5f. ", sizemax, qualmax));
						print(wyg);
					}

					imprrec(wyg, k0-1, k1+1);
				} else {
					j1 <- maxdiff;
				}

				j1 <- j1 + 1;
				if (v1-j1 < 0 || j1 > maxdiff || k1 > wn) break;
			}

			j0 <- j0 + 1;
			if (v0-j0 < 0 || j0 > maxdiff || k0 < 1) break;
		}
	}

	imprrec(wygmax, k0, k1);
# 	imprrec(wygmax, 0, k1);
# 	imprrec(wygmax, k0, wn+1);

	wyg <- wygmax;

# 	for (i in j0:j1)
# 	{
# 		lbound <- ifelse(i>1 && i < wn, min(wyg[i-1], wyg[i+1]), 0)
# 		if (wyg[i] > lbound)
# 		{
# 			repeat
# 			{
# 				wyg[i] <- wyg[i]-1L;
# 				POW <- .pareto2.htest.getsize_optimized(n, xhir0n, wyg, K, s, dhirp2mat);
# 				if (max(POW)>alpha)
# 				{
# 					wyg[i] <- wyg[i]+1L;
# 					break;
# 				}
# 				if (wyg[i] == lbound) break;
# 			}
#
# 			if (verbose)
# 			{
# 				POW <- .pareto2.htest.getsize_optimized(n, xhir0n, wyg, K, s, dhirp2mat);
# 				cat(sprintf("%3.0f%% complete in current iteration, q=%.4f.\r", abs((i-j0)/(j1-j0))*100, mean(POW)));
# 			}
# 		}
# 	}

	if (verbose)
	{
		POW <- .pareto2.htest.getsize_optimized(n, xhir0n, wyg, K, s, dhirp2mat);
		cat(sprintf("%3.0f%% complete in current iteration, q=%.4f.\r", 100, mean(POW)));
	}

	return(wyg);
}



#' /internal/
.pareto2.htest.acceptreg.improve <- function(n, wyg, j0, j1, alpha, verbose, K, s)
{
	epsbound <- 1e-7;
	wn <- n+1;


	# prepare input data for .pareto2.htest.getsize_optimized
	# this significanlty speeds up computations
	xhir0n <- (0:n);
	dhirp2mat <- matrix(nrow=length(K), ncol=length(xhir0n));
	for (i in 1:length(K))
		dhirp2mat[i,] <- pareto2.dhirsch(xhir0n,n,K[i],s);
	# -------------------------------------------------------



	for (i in j0:j1)
	{
		lbound <- ifelse(i>1 && i < wn, min(wyg[i-1], wyg[i+1]), 0)
		if (wyg[i] > lbound)
		{
			repeat
			{
				wyg[i] <- wyg[i]-1L;
				POW <- .pareto2.htest.getsize_optimized(n, xhir0n, wyg, K, s, dhirp2mat);
				if (max(POW)>alpha)
				{
					wyg[i] <- wyg[i]+1L;
					break;
				}
				if (wyg[i] == lbound) break;
			}

			if (verbose)
			{
				POW <- .pareto2.htest.getsize_optimized(n, xhir0n, wyg, K, s, dhirp2mat);
				cat(sprintf("%3.0f%% complete in current iteration, q=%.4f.\r", abs((i-j0)/(j1-j0))*100, mean(POW)));
			}
		}
	}

	if (verbose)
	{
		POW <- .pareto2.htest.getsize_optimized(n, xhir0n, wyg, K, s, dhirp2mat);
		cat(sprintf("%3.0f%% complete in current iteration, q=%.4f.\r", 100, mean(POW)));
	}

	return(wyg);
}


#' /internal/
.pareto2.htest.getK <- function(n, drho, s)
{
	# control function:
	kappa    <- function(x) { pmax(0,pmin(1,x))*n; }

	# kappa-indices (rho) for which to determine the power function
	RHO <- c(1e-3,seq(1e-3+drho, 1-drho-1e-3, drho), 1-1e-3);
	nrho <- length(RHO);

	stopifnot(RHO[1]<RHO[length(RHO)] && RHO[1]<RHO[2]);

	K   <- numeric(nrho);
	for (i in 1:nrho)
	{
		K[i] <- uniroot(function(k,s,targetrho,kappa)
			{
				1-ppareto2(kappa(targetrho),k,s)-targetrho;
			}, c(1e-15,ifelse(i==1,1e15,K[i-1])), # note that we assume RHO is sorted increasingly-that's much faster
			   s, RHO[i], kappa, tol=1e-20)$root;
	}

	return(K);
}



#' Performs \eqn{h}-test for equality of shape parameters
#' of two samples from the Pareto type-II distributions with known
#' and equal scale parameters, \eqn{s>0}.
#'
#' Given two equal-sized samples \eqn{X=(X_1,...,X_n)} i.i.d. \eqn{P2(k_x,s)}
#' and \eqn{Y=(Y_1,...,Y_m)} i.i.d. \eqn{P2(k_y,s)}
#' this test verifies the null hypothesis
#' \eqn{H_0: k_x=k_y}
#' against two-sided or one-sided alternatives, depending
#' on the value of \code{alternative}.
#' It bases on test statistic
#' \code{T=H(Y)-H(X)}
#' where \eqn{H} denotes Hirsch's \eqn{h}-index (see \code{\link{index.h}}).
#'
#' Note that for \eqn{k_x < k_y}, then \eqn{X} dominates \eqn{Y} stochastically.
#'
#' @title Two-sample h-test for equality of shape parameters for Type II-Pareto distributions with known common scale parameter
#' @param x an n-element non-negative numeric vector of data values.
#' @param y an n-element non-negative numeric vector of data values.
#' @param s scale parameter, \eqn{s>0}.
#' @param alternative indicates the alternative hypothesis and must be one of "two.sided" (default), "less", or "greater".
#' @param significance significance level. See Value for details.
#' @param wyg precomputed h-dependent acceptation region or \code{NULL}. See Value for details.
#' @param verbose logical; if \code{TRUE} then the computation progress will be printed out.
#' @param drho power calculation accuracy, a single number in [0.001, 0.1]. The smaller the value the slower computation, but more precise. This is used to determine \code{K} iff \code{K} is not given.
#' @param K numeric vector; shape parameters for which to calculate the power function or \code{NULL}.
#' @param improve logical; if \code{TRUE} then the greedy heuristic algorithm for improving the acceptation region will be run.
#' @return
#' The list of class \code{power.htest} with the following components is passed as a result:
#' \tabular{ll}{
#' \code{statistic} \tab	the value of the test statistic.\cr
#' \code{result} \tab	either FALSE (accept null hypothesis) or TRUE (reject).\cr
#' \code{alternative} \tab	a character string describing the alternative hypothesis.\cr
#' \code{method} \tab	a character string indicating what type of test was performed.\cr
#' \code{data.name} \tab	a character string giving the name(s) of the data.\cr
#' \code{wyg} \tab	a numeric vector giving the h-dependent acceptation region used.\cr
#' \code{size} \tab	size of the test corresponding to \code{wyg}.\cr
#' \code{qual} \tab	quality of the test corresponding to \code{wyg}, the closer to \code{significance}, the better.\cr
#' }
#' Currently no method for determining the p-value of this test is implemented.
#' @export
#' @seealso \code{\link{dpareto2}}, \code{\link{pareto2.goftest}}, \code{\link{pareto2.ftest}}, \code{\link{pareto2.htest.approx}}, \code{\link{index.h}}
#' @references
#' Gagolewski M., Grzegorzewski P., S-Statistics and Their Basic Properties, In: Borgelt C. et al (Eds.),
#' Combining Soft Computing and Statistical Methods in Data Analysis, Springer-Verlag, 2010, 281-288.\cr
pareto2.htest <- function(x, y, s, alternative = c("two.sided", "less", "greater"), significance=0.05, wyg=NULL, verbose=TRUE, drho=0.005, K=NULL, improve=TRUE)
{
	if (length(significance) != 1 || significance <= 0 || significance >= 1) stop("incorrect significance level");

	if (significance > 0.2) warning("'significance' is possibly incorrect");

	alternative <- match.arg(alternative);
	DNAME <- deparse(substitute(x));
	DNAME <- paste(DNAME, "and", deparse(substitute(y)));

	if (mode(s) != "numeric" || length(s) != 1 || s <= 0) stop("'s' should be > 0");

	if (length(drho) != 1 || drho < 0.001 || drho > 0.1)
		stop("drho should be a single numeric value in [0.000001, 0.1]");

	if (!is.null(K) && (mode(K) != "numeric" || length(K) < 10 || any(K<=0) || any(is.infinite(K))))
		stop("incorrect 'K'");

	if (mode(x) != "numeric" || mode(y) != "numeric") stop("non-numeric data given");

	x <- x[!is.na(x)];
	n <- length(x);
	if (n < 1L || any(x<0)) stop("incorrect 'x' data");

	y <- y[!is.na(y)];
	if (length(y) < 1L || any(y<0)) stop("incorrect 'y' data");

	if (length(y) != n) stop("non-equal-sized vectors given on input");

# 	if (n > 50 && is.null(wyg)) warning("n is large - Do you know what you're doing? It's damn slow! :-)");


	HY <- index.h(y,disable.check=TRUE);
	HX <- index.h(x,disable.check=TRUE);
	STATISTIC <- HY-HX;
	names(STATISTIC) <- "H";

	METHOD <- "Two-sample h-test for equality of shape parameters for Type II-Pareto distributions with known common scale parameter";

	nm_alternative <- switch(alternative, two.sided = "two-sided",
			less = "kx < ky",
			greater = "kx > ky");

	if (alternative == "two.sided") {
		powerscale <- 2;
	} else {
		powerscale <- 1;
	}
	alpha <- significance/powerscale;
	size <- NA;
	qual <- NA;

	wn <- n+1;     # number of possible h-index values

	# -----------------------------------------------------------------------

	if (is.null(wyg))
	{
		if (is.null(K))
		{
			# find scale parameters corresponding to kappa-indices
			if (verbose) cat(sprintf("Determining 'K' for 'drho'=%g...\n", drho));
			K <- .pareto2.htest.getK(n, drho, s)
		}


		# determine bound of the acceptation region that is independent on the value
		# of the h-indices (constant 'wyg')
		if (verbose) cat(sprintf("Determining h-independent bound of the acceptation region...\n"));

		v <- 0L; # initial solution
		wyg <- rep(v,wn) # a constant function of observed h-index
		POW <- .pareto2.htest.getpowerupper(n, wyg, K, K, s);

		while (max(POW) > alpha)
		{
			v <- v+1L;
			if (v > n) stop("h-independent bound could not be found");
			wyg <- rep(v,wn) # a constant function of observed h-index
			POW <- .pareto2.htest.getpowerupper(n, wyg, K, K, s);
		}

		if (verbose)
		{
			size <- max(POW)*powerscale;
			qual <- mean(POW)*powerscale;
			cat(sprintf("OK, v=%d for n=%d. This gives test size=%f and qual=%f.\n",
				v, n, size, qual));
		}
	} else { # wyg is not null

		if (length(wyg) != n+1 || mode(wyg) != "numeric" || any(wyg<0) || any(wyg>n))
			stop("incorrect 'wyg' for this test case");

		if (verbose || improve) # check whether given 'wyg' guarantees desired significance level
		{
			if (is.null(K))
			{
				# find scale parameters corresponding to kappa-indices
				if (verbose) cat(sprintf("Determining 'K' for 'drho'=%g...\n", drho));
				K <- .pareto2.htest.getK(n, drho, s);
			}

			POW <- .pareto2.htest.getpowerupper(n, wyg, K, K, s);
			if (verbose)
			{
				size <- max(POW)*powerscale;
				qual <- mean(POW)*powerscale;
				cat(sprintf("Given test size=%f and qual=%f for n=%d.\n",
					size, qual, n));
			}

			if (max(POW) > alpha) stop("Given 'wyg' does not guarantee desired significance level");
		}
	}


	# improve the acceptation region
	if (improve)
	{
		if (verbose) cat(sprintf("Improving h-dependent bounds of the acceptation region...\n"));

# 				wyg <- .htest.acceptregimprove(n, ppareto2, wyg, 1, wn,  alpha, verbose, K, s);
# 				wyg <- .htest.acceptregimprove(n, ppareto2, wyg, wn, 1,  alpha, verbose, K, s);
# 		wyg <- .pareto2.htest.acceptreg.improve(n, wyg, wn, 1, alpha, verbose, K, s)
# 		wyg <- .pareto2.htest.acceptreg.improve(n, wyg, 1, wn, alpha, verbose, K, s)
		wyg <- .pareto2.htest.acceptreg.improve2(n, wyg, alpha, verbose, K, s);

		if (verbose)
		{
			POW <- .pareto2.htest.getpowerupper(n, wyg, K, K, s);
			size <- max(POW)*powerscale;
			qual <- mean(POW)*powerscale;
			cat(sprintf("This now gives test size=%f and qual=%f.\n",
				size, qual));
			stopifnot(max(POW) <= alpha);
		}
	}


	# -----------------------------------------------------------------------

	if (alternative == "two.sided") {
		if (HX<HY)
		{
			RESULT <- (abs(STATISTIC)>wyg[HX+1]);
		} else {
			RESULT <- (abs(STATISTIC)>wyg[HY+1]);
		}
	} else if (alternative == "less") {
		if (HX<HY)
		{
			RESULT <- FALSE;
		} else {
			RESULT <- (abs(STATISTIC)>wyg[HY+1]);
		}
	} else {
		if (HX>HY)
		{
			RESULT <- FALSE;
		} else {
			RESULT <- (abs(STATISTIC)>wyg[HX+1]);
		}
	}

	RVAL <- list(statistic = STATISTIC, result = RESULT, alternative = nm_alternative,
		method = METHOD, data.name = DNAME, wyg = wyg, size = size, qual = qual);
	class(RVAL) <- "power.htest";
	return(RVAL);
}
