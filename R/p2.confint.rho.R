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


# #' /internal/
# .argmaxgreater <- function(f, interval, ..., lower=min(interval), upper=max(interval),
# 	f.lower = f(lower, ...), f.upper = f(upper, ...),
# 	tol = .Machine$double.eps^0.25, maxiter = 1000)
# {
# 	if (!missing(interval) && length(interval) != 2L)
# 		stop("'interval' must be a vector of length 2")
# 	if (!is.numeric(lower) || !is.numeric(upper) || lower >= upper)
# 		stop("lower < upper  is not fulfilled")
# 	if (is.na(f.lower))
# 		stop("f.lower = f(lower) is NA")
# 	if (is.na(f.upper))
# 		stop("f.upper = f(upper) is NA")
#
#
# 	if (!(f.lower < 0 && f.upper > 0))
# 		stop("f() values at end points not of opposite sign or a trivial solution")
#
# 	iter <- 1;
#
#
# 	while (iter < maxiter && (upper-lower) > tol)
# 	{
# 		mid <- (upper+lower)*0.5;
#
# 		if (f(mid, ...) > -1e-9) {
# 			upper <- mid;
# 		} else {
# 			lower <- mid;
# 		}
#
# 		iter <- iter + 1;
# 	}
#
# 	if (iter == maxiter)
# 		warning("_NOT_ converged in ", maxiter, " iterations");
#
# 	list(root = upper, f.root = f(upper, ...), iter = iter,
#         estim.prec = upper-lower)
# }
#
#
# #' /internal/
# .argmaxle <- function(f, interval, ..., lower=min(interval), upper=max(interval),
# 	f.lower = f(lower, ...), f.upper = f(upper, ...),
# 	tol = .Machine$double.eps^0.25, maxiter = 1000)
# {
# 	if (!missing(interval) && length(interval) != 2L)
# 		stop("'interval' must be a vector of length 2")
# 	if (!is.numeric(lower) || !is.numeric(upper) || lower >= upper)
# 		stop("lower < upper  is not fulfilled")
# 	if (is.na(f.lower))
# 		stop("f.lower = f(lower) is NA")
# 	if (is.na(f.upper))
# 		stop("f.upper = f(upper) is NA")
#
#
# 	if (!(f.lower < 0 && f.upper > 0))
# 		stop("f() values at end points not of opposite sign or a trivial solution")
#
# 	iter <- 1;
#
# 	while (iter < maxiter && (upper-lower) > tol)
# 	{
# 		mid <- (upper+lower)*0.5;
#
# 		if (f(mid, ...) <= 0.0) {
# 			lower <- mid;
# 		} else {
# 			upper <- mid;
# 		}
#
# 		iter <- iter + 1;
# 	}
#
# 	if (iter == maxiter)
# 		warning("_NOT_ converged in ", maxiter, " iterations");
#
# 	list(root = lower, f.root = f(lower, ...), iter = iter,
#         estim.prec = upper-lower)
# }






#' Computes the exact right-sided confidence interval for the \eqn{\rho}-index of
#' a probability distribution in an \eqn{(X_1,\dots,X_n)} i.i.d. Pareto-type II
#' model with known scale parameter \eqn{s>0}.
#' The confidence interval bases on the observed value
#' of S-statistic w.r.t. to the given control function \eqn{\kappa}.
#'
#'
#' For more information see man page on  \code{\link{rho.get}}, \code{\link{Sstat}} and the paper (Gagolewski, Grzegorzewski, 2010).
#'
#' @references
#' Gagolewski M., Grzegorzewski P., S-Statistics and Their Basic Properties, In: Borgelt C. et al (Eds.),
#' Combining Soft Computing and Statistical Methods in Data Analysis, Springer-Verlag, 2010, 281-288.\cr
#'
#' @title Right-sided exact confidence interval for the rho-index
#' @param v observed value of the S-statistic w.r.t. \eqn{\kappa}.
#' @param kappa an increasing function, \eqn{\kappa}, a so-called control function.
#' @param s scale parameter, \eqn{s>0}.
#' @param n sample size.
#' @param conf.level confidence level; defaults 0.95.
#' @param tol the desired accuracy (convergence tolerance).
#' @return Upper bound of the confidence interval.
#' @seealso \code{\link{ppareto2}}, \code{\link{pareto2.confint.rho.approx}}, \code{\link{Sstat}}, \code{\link{pareto2.confint.rho.lower}},
#' \code{\link{pareto2.confint.rho}}, \code{\link{rho.get}}
#' @export
pareto2.confint.rho.upper <- function(v, kappa, s, n, conf.level=0.95, tol=1e-12)
{
	# rho_lower = min{rho: D_n,rho (v-) <= gamma}
	# will be found by:
	# klower = max{k: psstat(v-,n,ppareto2,kappa,k,s)<=gamma}

	gamma <- 1-conf.level;

	if (mode(s) != "numeric" || length(s) != 1 || s <= 0) stop("'s' should be > 0");

	if (!is.numeric(v) || length(v) != 1)
		stop("v must be a single numeric value")

# 	if (v < 1e-13) return(0.0);
	if (v > 1-1e-6) return(1.0);
	if (gamma < 1e-9) return(1.0);


	xsol <- uniroot(function(x,v,kappa,s,n,gamma) {
		psstat(v,n,ppareto2,kappa,x,s)-gamma;
	}, c(0,1e25),v,kappa,s,n,gamma, tol=tol, maxiter=1000)$root;

# 	xsol <- .argmaxle(function(x,v,kappa,s,n,gamma) {
# 		psstat(v,n,ppareto2,kappa,x,s)-gamma;
# 	}, c(max(0,xsol_initial-1e-6),xsol_initial+1e-6),v,kappa,s,n,gamma, tol=tol, maxiter=10000)$root;

	rho.get(ppareto2, kappa, xsol, s,tol=tol); # return value
}





#' Computes the exact left-sided confidence interval for the \eqn{\rho}-index of
#' a probability distribution in an \eqn{(X_1,\dots,X_n)} i.i.d. Pareto-type II
#' model with known scale parameter \eqn{s>0}.
#' The confidence interval bases on the observed value
#' of S-statistic w.r.t. to the given control function \eqn{\kappa}.
#'
#'
#' For more information see man page on  \code{\link{rho.get}}, \code{\link{Sstat}} and the paper (Gagolewski, Grzegorzewski, 2010).
#'
#' @references
#' Gagolewski M., Grzegorzewski P., S-Statistics and Their Basic Properties, In: Borgelt C. et al (Eds.),
#' Combining Soft Computing and Statistical Methods in Data Analysis, Springer-Verlag, 2010, 281-288.\cr
#'
#' @title Left-sided exact confidence interval for the rho-index
#' @param v observed value of the S-statistic w.r.t. \eqn{\kappa}.
#' @param kappa an increasing function, \eqn{\kappa}, a so-called control function.
#' @param s scale parameter, \eqn{s>0}.
#' @param n sample size.
#' @param conf.level confidence level; defaults 0.95.
#' @param tol the desired accuracy (convergence tolerance).
#' @return Lower bound of the confidence interval.
#' @seealso \code{\link{ppareto2}}, \code{\link{pareto2.confint.rho.approx}}, \code{\link{Sstat}},
#' \code{\link{pareto2.confint.rho.upper}},
#' \code{\link{pareto2.confint.rho}}, \code{\link{rho.get}}
#' @export
pareto2.confint.rho.lower <- function(v, kappa, s, n, conf.level=0.95, tol=1e-12)
{
	# rho_lower = max{rho: D_n,rho (v) >= 1-gamma}
	# will be found by:
	# klower = min{k: psstat(v,n,ppareto2,kappa,k,s)>=1-gamma}

	gamma <- 1-conf.level;

	if (mode(s) != "numeric" || length(s) != 1 || s <= 0) stop("'s' should be > 0");

	if (!is.numeric(v) || length(v) != 1)
		stop("v must be a single numeric value")

	if (v < 1e-6) return(0.0);
# 	if (v > 1-1e-13) return (1.0);
	if (gamma < 1e-9) return(0.0);

	v <- max(0,v-1e-12);

	xsol <- uniroot(function(x,v,kappa,s,n,gamma) {
		psstat(v,n,ppareto2,kappa,x,s)-1+gamma;
	}, c(0,1e25),v,kappa,s,n,gamma, tol=tol, maxiter=1000)$root;

# 	xsol <- .argmaxgreater(function(x,v,kappa,s,n,gamma) {
# 		psstat(v,n,ppareto2,kappa,x,s)-1+gamma;
# 	}, c(max(0,xsol_initial-1e-6),xsol_initial+1e-6),v,kappa,s,n,gamma, tol=tol, maxiter=10000)$root;

	rho.get(ppareto2, kappa, xsol, s,tol=tol); # return value
}






#' Computes the exact two-sided confidence interval for the \eqn{\rho}-index of
#' a probability distribution in an \eqn{(X_1,\dots,X_n)} i.i.d. Pareto-type II
#' model with known scale parameter \eqn{s>0}.
#' The confidence interval bases on the observed value
#' of S-statistic w.r.t. to the given control function \eqn{\kappa}.
#'
#'
#' For more information see man page on  \code{\link{rho.get}}, \code{\link{Sstat}} and the paper (Gagolewski, Grzegorzewski, 2010).
#'
#' @references
#' Gagolewski M., Grzegorzewski P., S-Statistics and Their Basic Properties, In: Borgelt C. et al (Eds.),
#' Combining Soft Computing and Statistical Methods in Data Analysis, Springer-Verlag, 2010, 281-288.\cr
#'
#' @title Two-sided exact confidence interval for the rho-index
#' @param v observed value of the S-statistic w.r.t. \eqn{\kappa}.
#' @param kappa an increasing function, \eqn{\kappa}, a so-called control function.
#' @param s scale parameter, \eqn{s>0}.
#' @param n sample size.
#' @param conf.level confidence level; defaults 0.95.
#' @param tol the desired accuracy (convergence tolerance).
#' @return Vector of length 2 with the computed bounds of the confidence interval.
#' @seealso \code{\link{ppareto2}}, \code{\link{pareto2.confint.rho.approx}}, \code{\link{Sstat}},
#' \code{\link{pareto2.confint.rho.lower}},  \code{\link{pareto2.confint.rho.upper}},
#' \code{\link{pareto2.confint.rho}}, \code{\link{rho.get}}
#' @export
pareto2.confint.rho <- function(v, kappa, s, n, conf.level=0.95, tol=1e-12)
{
	gamma <- 1-conf.level;
	return(c(
		pareto2.confint.rho.lower(v,kappa,s,n,1-gamma*0.5,tol),
		pareto2.confint.rho.upper(v,kappa,s,n,1-gamma*0.5,tol)
	));
}



