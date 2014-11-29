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


#' Random deviates generation for the Pareto Type-II (Lomax)  distribution with shape
#' parameter equal to \eqn{k>0} and scale parameter equal to \eqn{s>0}.
#'
#' @title Pareto distribution of the second kind - random deviates
#' @param n number of observations.
#' @param k vector of shape parameters, \eqn{k>0}.
#' @param s vector of scale parameters, \eqn{s>0}.
#' @return The function returns generated pseudorandom deviates.
#' @export
#' @seealso \code{\link{dpareto2}}, \code{\link{ppareto2}}, \code{\link{qpareto2}}, \code{\link{pareto2.zsestimate}}, \code{\link{pareto2.mlekestimate}}, \code{\link{pareto2.goftest}}, \code{\link{pareto2.ftest}}
rpareto2 <- function(n, k=1, s=1)
{
	s*((runif(n)^(-1/k)) - 1); # s*(exp(-log(runif(n))/k)-1); - faster?
}


#' Cumulative distribution function
#' for the Pareto Type-II (Lomax)  distribution with shape parameter equal to \eqn{k>0}
#' and scale parameter equal to \eqn{s>0}.
#'
#' The c.d.f. at \eqn{x\ge 0} is given by \code{F(x)=1-s^k/(s+x)^k}.
#'
#' @title Pareto distribution of the second kind - c.d.f.
#' @param q vector of quantiles.
#' @param k vector of shape parameters, \eqn{k>0}.
#' @param s vector of scale parameters, \eqn{s>0}.
#' @return The function gives the c.d.f..
#' @export
#' @seealso \code{\link{dpareto2}}, \code{\link{qpareto2}}, \code{\link{rpareto2}}, \code{\link{pareto2.zsestimate}}, \code{\link{pareto2.mlekestimate}}, \code{\link{pareto2.goftest}}, \code{\link{pareto2.ftest}}
ppareto2 <- function(q, k=1, s=1)
{
	ifelse(q<0, 0, (1-(s/(s+q))^k));
}


#' Quantile function
#' for the Pareto Type-II (Lomax) distribution with shape parameter equal to \eqn{k>0}
#' and scale parameter equal to \eqn{s>0}.
#'
#' @title Pareto distribution of the second kind - quantiles
#' @param p vector of probabilities, \eqn{p\in(0,1)}.
#' @param k vector of shape parameters, \eqn{k>0}.
#' @param s vector of scale parameters, \eqn{s>0}.
#' @return The function gives the theoretical quantiles.
#' @export
#' @seealso \code{\link{dpareto2}}, \code{\link{ppareto2}}, \code{\link{rpareto2}}, \code{\link{pareto2.zsestimate}}, \code{\link{pareto2.mlekestimate}}, \code{\link{pareto2.goftest}}, \code{\link{pareto2.ftest}}
qpareto2 <- function(p, k=1, s=1)
{
	ifelse(p<0 | p>1, NA, s*((1-p)^(-1/k)-1));
}


#' Density function for the Pareto Type-II (Lomax)  distribution with shape parameter
#' equal to \eqn{k>0} and scale parameter equal to \eqn{s>0}.
#'
#' The p.d.f. at \eqn{x\ge 0} is given by \code{f(x)=k*s^k/(s+x)^(k+1)}.
#'
#' @title Pareto distribution of the second kind - density
#' @param x vector of quantiles.
#' @param k vector of shape parameters, \eqn{k>0}.
#' @param s vector of scale parameters, \eqn{s>0}.
#' @return The function gives the density.
#' @export
#' @seealso \code{\link{ppareto2}}, \code{\link{qpareto2}}, \code{\link{rpareto2}}, \code{\link{pareto2.zsestimate}}, \code{\link{pareto2.mlekestimate}}, \code{\link{pareto2.goftest}}, \code{\link{pareto2.ftest}}
dpareto2 <- function(x, k=1, s=1)
{
	ifelse(x<=0, 0, k/(s+x)*(s/(s+x))^k);
}
