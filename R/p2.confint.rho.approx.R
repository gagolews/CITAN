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


#' Computes the approximate (asymptotic) left-sided confidence interval for the \eqn{\rho}-index of
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
#' @title Left-sided approximate confidence interval for the rho-index
#' @param v observed value of the S-statistic w.r.t. \eqn{\kappa}.
#' @param kappa an increasing function, \eqn{\kappa}, a so-called control function.
#' @param kappaInvDer the derivative of the inverse of \eqn{\kappa}.
#' @param s scale parameter, \eqn{s>0}.
#' @param n sample size.
#' @param conf.level confidence level; defaults 0.95.
#' @param tol the desired accuracy (convergence tolerance).
#' @return Lower bound of the confidence interval.
#' @seealso \code{\link{ppareto2}}, \code{\link{pareto2.confint.rho}}, \code{\link{Sstat}},
#' \code{\link{pareto2.confint.rho.approx.upper}},
#' \code{\link{pareto2.confint.rho.approx}}, \code{\link{rho.get}}
#' @export
pareto2.confint.rho.approx.lower <- function(v, kappa, kappaInvDer, s, n, conf.level=0.95, tol=1e-20)
{
	gamma <- 1-conf.level;

	if (!is.numeric(v) || length(v) != 1)
		stop("v must be a single numeric value");

	if (mode(s) != "numeric" || length(s) != 1 || s <= 0) stop("'s' should be > 0");

	if (v < 1e-6) return(0.0);
# 	if (v > 1-1e-13) return (1.0);
	if (gamma < 1e-5) return(0.0);


	bord <- uniroot(function(rho, v, kappa, kappaInvDer, s, n, gamma)
	{
		k <- uniroot(function(k, s, kappa, rho) {
				1-ppareto2(kappa(rho), k, s)-rho
			}, c(1e-15,1e10), s, kappa, rho, tol=tol)$root;

# 		k <- log(rho)/(log(s/(s+rho))); # round-off errors :-(
# 		k <- log1p(rho-1)/(log1p(s/(s+rho)-1)); # round-off errors :-(

		gprimerho <- dpareto2(kappa(rho), k, s)/abs(kappaInvDer(kappa(rho)));

		qnorm(1-gamma, rho, sqrt(rho*(1-rho)/n)/(1+gprimerho))-v;
	}, c(1e-7,1-1e-7), v, kappa, kappaInvDer, s, n, gamma, tol=tol)$root;

	return(bord);

# 	return(qnorm(gamma/2, v, sqrt(v*(1-v)/n)/(1+gprime)));
}






#' Computes the approximate (asymptotic) right-sided confidence interval for the \eqn{\rho}-index of
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
#' @title Right-sided approximate confidence interval for the rho-index
#' @param v observed value of the S-statistic w.r.t. \eqn{\kappa}.
#' @param kappa an increasing function, \eqn{\kappa}, a so-called control function.
#' @param kappaInvDer the derivative of the inverse of \eqn{\kappa}.
#' @param s scale parameter, \eqn{s>0}.
#' @param n sample size.
#' @param conf.level confidence level; defaults 0.95.
#' @param tol the desired accuracy (convergence tolerance).
#' @return Upper bound of the confidence interval.
#' @seealso \code{\link{ppareto2}}, \code{\link{pareto2.confint.rho}}, \code{\link{Sstat}},
#' \code{\link{pareto2.confint.rho.approx.lower}},
#' \code{\link{pareto2.confint.rho.approx}}, \code{\link{rho.get}}
#' @export
pareto2.confint.rho.approx.upper <- function(v, kappa, kappaInvDer, s, n, conf.level=0.95, tol=1e-20)
{
	gamma <- 1-conf.level;

	if (!is.numeric(v) || length(v) != 1)
		stop("v must be a single numeric value");

	if (mode(s) != "numeric" || length(s) != 1 || s <= 0) stop("'s' should be > 0");

	if (v > 1-1e-6) return(1.0);
# 	if (v < 1e-13) return (0.0);
	if (gamma < 1e-5) return(1.0);


	bord <- uniroot(function(rho, v, kappa, kappaInvDer, s, n, gamma)
	{
		k <- uniroot(function(k, s, kappa, rho) {
				1-ppareto2(kappa(rho), k, s)-rho
			}, c(1e-15,1e10), s, kappa, rho, tol=tol)$root;

# 		k <- log(rho)/(log(s/(s+rho))); # round-off errors :-(
# 		k <- log1p(rho-1)/(log1p(s/(s+rho)-1)); # round-off errors :-(

		gprimerho <- dpareto2(kappa(rho), k, s)/abs(kappaInvDer(kappa(rho)));

		qnorm(gamma, rho, sqrt(rho*(1-rho)/n)/(1+gprimerho))-v;
	}, c(1e-7,1-1e-7), v, kappa, kappaInvDer, s, n, gamma, tol=tol)$root;

	return(bord);

# 	return(qnorm(gamma/2, v, sqrt(v*(1-v)/n)/(1+gprime)));
}







#' Computes the approximate (asymptotic) two-sided confidence interval for the \eqn{\rho}-index of
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
#' @title Two-sided approximate confidence interval for the rho-index
#' @param v observed value of the S-statistic w.r.t. \eqn{\kappa}.
#' @param kappa an increasing function, \eqn{\kappa}, a so-called control function.
#' @param kappaInvDer the derivative of the inverse of \eqn{\kappa}.
#' @param s scale parameter, \eqn{s>0}.
#' @param n sample size.
#' @param conf.level confidence level; defaults 0.95.
#' @param tol the desired accuracy (convergence tolerance).
#' @return Vector of length 2 with the computed bounds of the confidence interval.
#' @seealso \code{\link{ppareto2}}, \code{\link{pareto2.confint.rho}}, \code{\link{Sstat}},
#' \code{\link{pareto2.confint.rho.approx.lower}},
#' \code{\link{pareto2.confint.rho.approx.upper}}, \code{\link{rho.get}}
#' @export
pareto2.confint.rho.approx <- function(v, kappa, kappaInvDer, s, n, conf.level=0.95, tol=1e-20)
{
	gamma <- 1-conf.level;
	return(c(
		pareto2.confint.rho.approx.lower(v,kappa, kappaInvDer,s,n,1-gamma*0.5,tol),
		pareto2.confint.rho.approx.upper(v,kappa, kappaInvDer,s,n,1-gamma*0.5,tol)
	));
}
