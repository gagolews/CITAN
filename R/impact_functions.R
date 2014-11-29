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


#' Computes the "Classical" \eqn{h}-index of a numeric vector.
#'
#' Given a sequence of \eqn{n} non-negative numbers \eqn{x=(x_1,\dots,x_n)},
#' where \eqn{x_i \ge x_j} for \eqn{i \le j},
#' the \dfn{\eqn{h}-index} (Hirsch, 2005) for \eqn{x} is defined as
#' \deqn{H(x)=\max\{i=1,\dots,n: x_i \ge i\}}{H(x)=max{i=1,\dots,n: x_i \ge i}}
#' if \eqn{n \ge 1} and \eqn{x_1 \ge 1}, or \eqn{H(x)=0} otherwise.
#'
#' If \code{disable.check} is set to \code{FALSE}, then
#' eventual \code{NA} values are removed from the input vector.
#'
#' If a non-increasingly sorted vector is given as input (set \code{sorted.dec} to \code{TRUE})
#' the value is calculated in linear or log-time, depending on the value of the \code{algorithm} parameter.
#'
#' @references Hirsch J.E., An index to quantify individual's scientific research output, Proceedings of the National Academy of Sciences 102(46), 16569-16572, 2005.\cr
#'
#' @title Hirsch's h-index
#' @param x a non-negative numeric vector.
#' @param sorted.dec logical; \code{TRUE} if the vector has already been sorted non-increasingly; defaults \code{FALSE}.
#' @param disable.check logical; \code{TRUE} to disable some validity checks on the input vector; defaults \code{FALSE}.
#' @param algorithm type of algorithm, "linear-time" or "log-time" (default).
#' @return The function returns a single number or NA if improper input has been given.
#' @seealso
#' \code{\link{pareto2.confint.h}}, \code{\link{pareto2.htest}},
#' \code{\link{pareto2.htest.approx}},\cr
#' \code{\link{phirsch}}, \code{\link{dhirsch}}, 
#' \code{\link{index.g}}, \code{\link{index.rp}}, \code{\link{index.lp}},
#' \code{\link{Sstat}}, \code{\link{Sstat2}}
#'
#' @examples
#' authors <- list(  # a list of numeric sequences
#'                   # (e.g. citation counts of the articles
#'                   # written by some authors)
#'     "A" =c(23,21,4,2,1,0,0),
#'     "B" =c(11,5,4,4,3,2,2,2,2,2,1,1,1,0,0,0,0),
#'     "C" =c(53,43,32,23,14,13,12,8,4,3,2,1,0)
#'  );
#' index.h(authors$A);
#' lapply(authors, index.h);
#' @export
index.h <- function(x, sorted.dec=FALSE, disable.check=FALSE, algorithm=c("log-time", "linear-time"))
{
	algorithm <- match.arg(algorithm);

	if (!disable.check)
	{
		if (length(x) == 0) return(0);
		if (mode(x) != "numeric") return(NA);
		if (any(x < 0)) return(NA);
		x <- x[!is.na(x)];
	}

	if (!sorted.dec)
		x <- sort(x, decreasing=TRUE);

	if (algorithm == "linear-time") {
		.C("index_h", as.double(x), as.integer(length(x)), out=double(1), DUP=FALSE, PACKAGE="CITAN")$out;
	} else {
		# (may be slower):
		.C("index_h_log", as.double(x), as.integer(length(x)), out=double(1), DUP=FALSE, PACKAGE="CITAN")$out;
	}
}



#' Computes the "Classical" \eqn{g}-index of a numeric vector.
#'
#' Given a sequence of \eqn{n} non-negative numbers \eqn{x=(x_1,\dots,x_n)},
#' where \eqn{x_i \ge x_j} for \eqn{i \le j},
#' the \dfn{\eqn{g}-index} (Egghe, 2006) for \eqn{x} is defined as
#' \deqn{G(x)=\max\{i=1,\dots,n: \sum_{j=1}^i x_i \ge i^2\},}{G(x)=max{i=1,\dots,n: x_1+\dots+x_i \ge i^2}}
#' if \eqn{n \ge 1} and \eqn{x_1 \ge 1}, or \eqn{G(x)=0} otherwise.
#'
#' If \code{disable.check} is set to \code{FALSE}, then
#' eventual \code{NA} values are removed from the input vector.
#'
#' If a non-increasingly sorted vector is given as input (set \code{sorted.dec} to \code{TRUE})
#' the value is calculated in linear time.
#'
#' @references Egghe L., Theory and practise of the g-index, Scientometrics 69(1), 131-152, 2006.\cr
#'
#' @title Egghe's g-index
#' @param x a non-negative numeric vector.
#' @param sorted.dec logical; \code{TRUE} if the vector has already been sorted non-increasingly; defaults \code{FALSE}.
#' @param disable.check logical; \code{TRUE} to disable some validity checks on the input vector; defaults \code{FALSE}.
#' @return The function returns a single number or NA if improper input has been given.
#' @seealso \code{\link{index.h}}, \code{\link{index.rp}}, \code{\link{index.lp}}, \code{\link{Sstat}}, \code{\link{Sstat2}}
#' @export
index.g <- function(x, sorted.dec=FALSE, disable.check=FALSE)
{
	if (!disable.check)
	{
		if (length(x) == 0) return(0);
		if (mode(x) != "numeric") return(NA);
		if (any(x < 0)) return(NA);
		x <- x[!is.na(x)];
	}

	if (!sorted.dec)
		x <- sort(x, decreasing=TRUE);

	.C("index_g", as.double(x), as.integer(length(x)), out=double(1), DUP=FALSE, PACKAGE="CITAN")$out;
}



#' Computes the \eqn{r_p}-index of a numeric vector for given \eqn{p}.
#'
#' Given a sequence of \eqn{n} non-negative numbers \eqn{x=(x_1,\dots,x_n)},
#' where \eqn{x_i \ge x_j} for \eqn{i \le j},
#' the \dfn{\eqn{r_p}-index} for \eqn{p=\infty} equals to
#' \deqn{r_p(x)=\max_{i=1,\dots,n} \{ \min\{i,x_i\} \}}{r_p(x) = max{ min{i, x_i} } for i=1,\dots,n}
#' if \eqn{n \ge 1}, or \eqn{r_\infty(x)=0} otherwise.
#' For the definition of the \eqn{r_p}-index for \eqn{p < \infty} we refer
#' to (Gagolewski, Grzegorzewski, 2009).
#'
#' Note that if \eqn{x_1,\dots,x_n} are integers, then
#' \deqn{r_\infty(x)=H(x),} where \eqn{H} is the \eqn{h}-index (Hirsch, 2005) and
#' \deqn{r_1(x)=W(x),} where \eqn{W} is the \eqn{w}-index (Woeginger, 2008).
#'
#' If \code{disable.check} is set to \code{FALSE}, then
#' eventual \code{NA} values are removed from the input vector.
#'
#' If a non-increasingly sorted vector is given as input (set \code{sorted.dec} to \code{TRUE})
#' the value of \eqn{r_\infty} is calculated in log time
#' (note that it may be determined in linear time using \code{max(pmin(x, 1:length(x)))}).
#' Otherwise, linear time is needed.
#'
#' @references
#' Gagolewski M., Grzegorzewski P., A geometric approach to the construction of scientific impact indices, Scientometrics, 81(3), 2009, pp. 617-634.\cr
#' Hirsch J.E., An index to quantify individual's scientific research output, Proceedings of the National Academy of Sciences 102(46), 16569-16572, 2005.\cr
#' Woeginger G.J., An axiomatic characterization of the Hirsch-index, Mathematical Social Sciences, 56(2), 224-232, 2008.\cr
#'
#' @title The r_p-index
#' @param x a non-negative numeric vector.
#' @param p index order, \eqn{p \in [1,\infty]}{p in [1,\infty]}; defaults \eqn{\infty} (\code{Inf}).
#' @param sorted.dec logical; \code{TRUE} if the vector has already been sorted non-increasingly; defaults \code{FALSE}.
#' @param disable.check logical; \code{TRUE} to disable some validity checks on the input vector; defaults \code{FALSE}.
#' @return The function returns a single number or NA if improper input has been given.
#' @seealso \code{\link{index.h}}, \code{\link{index.g}}, \code{\link{index.lp}}, \code{\link{Sstat}}, \code{\link{Sstat2}}
#' @examples
#' x <- runif(100, 0, 100);
#' index.rp(x);            # the r_oo-index
#' floor(index.rp(x));     # the h-index
#' index.rp(floor(x), 1);  # the w-index
#' @export
index.rp <- function(x, p=Inf, sorted.dec=FALSE, disable.check=FALSE)
{
	if (!disable.check)
	{
		if (length(x) == 0) return(0);
		if (mode(x) != "numeric") return(NA);
		if (any(x < 0)) return(NA);
		x <- x[!is.na(x)];
	}

	if (mode(p) != "numeric" || length(p)!=1 || p < 1) stop("'p' should be a single numeric value >= 1");

	if (!sorted.dec)
		x <- sort(x, decreasing=TRUE);

	if (is.finite(p))
	{
		if (p > 50) warning("'p' is quite large. possible accuracy problems. maybe you should try 'p'==Inf?");
		.C("index_rp_finite", as.double(x), as.integer(length(x)), as.double(p), out=double(1), DUP=FALSE, PACKAGE="CITAN")$out;
	} else
	{
		.C("index_rp_infinite", as.double(x), as.integer(length(x)), out=double(1), DUP=FALSE, PACKAGE="CITAN")$out;
	}
}




#' Computes the \eqn{l_p}-index of a numeric vector for given \eqn{p}.
#'
#' Given a sequence of \eqn{n} non-negative numbers \eqn{x=(x_1,\dots,x_n)},
#' where \eqn{x_i \ge x_j} for \eqn{i \le j},
#' the \dfn{\eqn{l_p}-index} for \eqn{p=\infty} equals to
#' \deqn{l_p(x)=\arg\max_{(i,x_i), i=1,\dots,n} \{ i x_i \}}{l_p(x) = arg max_(i,x_i) { i*x_i } for i=1,\dots,n}
#' if \eqn{n \ge 1}, or \eqn{l_\infty(x)=0} otherwise.
#' Note that if \eqn{(i,x_i)=l_\infty(x)}, then
#' \deqn{MAXPROD(x) = i x_i,}{MAXPROD(x) = i*x_i,} where \eqn{MAXPROD} is the index proposed in (Kosmulski, 2007).
#'
#' For the definition of the \eqn{l_p}-index for \eqn{p < \infty} we refer
#' to (Gagolewski, Grzegorzewski, 2009a).
#'
#' If \code{disable.check} is set to \code{FALSE}, then
#' eventual \code{NA} values are removed from the input vector.
#'
#' If a non-increasingly sorted vector is given as input (set \code{sorted.dec} to \code{TRUE})
#' the result is computed in linear time (see Gagolewski, Debski, Nowakiewicz, 2009b).
#'
#' @references
#' Gagolewski M., Grzegorzewski P., A geometric approach to the construction of scientific impact indices, Scientometrics, 81(3), 2009a, pp. 617-634.\cr
#' Gagolewski M., Debski M., Nowakiewicz M., Efficient algorithms for computing ''geometric'' scientific impact indices, Research Report of Systems Research Institute, Polish Academy of Sciences RB/1/2009, 2009b.\cr
#' Kosmulski M., MAXPROD - A new index for assessment of the scientific output of an individual, and a comparison with the h-index, Cybermetrics, 11(1), 2007.\cr
#'
#' @title The l_p-index
#' @param x a non-negative numeric vector.
#' @param p index order, \eqn{p \in [1,\infty]}{p in [1,\infty]}; defaults \eqn{\infty} (\code{Inf}).
#' @param sorted.dec logical; \code{TRUE} if the vector has already been sorted non-increasingly; defaults \code{FALSE}.
#' @param disable.check logical; \code{TRUE} to disable some validity checks on the input vector; defaults \code{FALSE}.
#' @return The function returns a numeric vector of length 2 equal to \eqn{(i,x_i)} or NA if improper input has been given.
#' @seealso \code{\link{index.h}}, \code{\link{index.g}}, \code{\link{index.rp}}, \code{\link{Sstat}}, \code{\link{Sstat2}}
#' @examples
#' x <- runif(100, 0, 100);
#' index.lp(x);                # two-dimensional value, can not be used
#'                             # directly in the analysis
#' prod(index.lp(x));          # the MAXPROD-index (one-dimensional)
#' mean(index.lp(x,1));        # some other one-dimensional impact index
#' @export
index.lp <- function(x, p=Inf, sorted.dec=FALSE, disable.check=FALSE)
{
	if (!disable.check)
	{
		if (length(x) == 0) return(0);
		if (mode(x) != "numeric") return(NA);
		if (any(x < 0)) return(NA);
		x <- x[!is.na(x)];
	}

	if (mode(p) != "numeric" || length(p)!=1 || p < 1) stop("'p' should be a single numeric value >= 1");

	if (!sorted.dec)
		x <- sort(x, decreasing=TRUE);

	if (is.finite(p))
	{
		if (p > 50) warning("'p' is quite large. possible accuracy problems. maybe you should try 'p'==Inf?");
		.C("index_lp_finite", as.double(x), as.integer(length(x)), as.double(p), integer(length(x)+1), out=double(2), DUP=FALSE, PACKAGE="CITAN")$out;
	} else
	{
		.C("index_lp_infinite", as.double(x), as.integer(length(x)), out=double(2), DUP=FALSE, PACKAGE="CITAN")$out;
	}
}



#' Computes the S-statistic2 w.r.t. to the identity function for data transformed by the inverse of a control function.
#'
#' Given a sequence of \eqn{n} non-negative numbers \eqn{x=(x_1,\dots,x_n)},
#' where \eqn{x_i\ge x_j} for \eqn{i\le j},
#' and a nondecreasing function \eqn{\gamma: R\to[0,1]}{\gamma: R->[0,1]},
#' the \dfn{S-statistic2} (Gagolewski, Grzegorzewski, 2010) for \eqn{x} is defined as
#' \deqn{V_n(x)=\max_{i=1,\dots,n}\{\min\{\gamma(x_i), i/n \}\}}{V_n(x)=max{ min{i/n, \gamma(x_i)} } for i=1,\dots,n}
#'
#' If \code{disable.check} is set to \code{FALSE}, then
#' eventual \code{NA} values are removed from the input vector.
#'
#' If a non-increasingly sorted vector is given as input (set \code{sorted.dec} to \code{TRUE})
#' the result is computed in linear time.
#'
#' @references Gagolewski M., Grzegorzewski P., S-Statistics and Their Basic Properties,
#' In: Borgelt C. et al (Eds.), Combining Soft Computing and Statistical Methods in Data Analysis, (SMPS 2010), Springer-Verlag, 2010, 281-288.
#'
#' @title S-statistic2
#' @param x a vector of real numbers.
#' @param kappaInv a nondecreasing function ranging on [0,1], \eqn{\gamma} (see Details), the inverse of a so-called control function.
#' @param sorted.dec logical; \code{TRUE} if the vector has already been sorted non-increasingly; defaults \code{FALSE}.
#' @param disable.check logical; \code{TRUE} to disable some validity checks on the input vector; defaults \code{FALSE}.
#' @return The function returns a single number or NA if improper input has been given.
#' @seealso \code{\link{index.h}}, \code{\link{index.g}}, \code{\link{index.rp}}, \code{\link{index.lp}}, \code{\link{Sstat}}, \code{\link{psstat}}, \code{\link{dsstat}}
#' @examples
#' x <- rpareto2(25, 1.05, 1);
#' kappaInv <- function(x) { pmax(0,pmin(1,x/25)); }
#' Sstat2(x, kappaInv, FALSE, TRUE);
#' @export
Sstat2 <- function(x, kappaInv, sorted.dec=FALSE, disable.check=FALSE)
{
	if (!disable.check)
	{
		if (length(x) == 0) return(0);
		if (mode(x) != "numeric") return(NA);
		if (any(x < 0)) return(NA);
		x <- x[!is.na(x)];
	}

	if (!sorted.dec)
		x <- sort(x, decreasing=TRUE);

	# internal method needs O(log n) time, however kappaInv(x) is O(n)
	.C("Sstat2", as.double(kappaInv(x)), as.integer(length(x)), out=double(1), DUP=FALSE, PACKAGE="CITAN")$out;
}


#' Computes the S-statistic w.r.t. to a control function.
#'
#' Given a sequence of \eqn{n} non-negative numbers \eqn{x=(x_1,\dots,x_n)},
#' where \eqn{x_i\ge x_j} for \eqn{i\le j},
#' and an increasing function \eqn{\kappa: [0,1]\to[a,b]}{\kappa: [0,1]->[a,b]} for some \eqn{a,b},
#' the \dfn{S-statistic} (Gagolewski, Grzegorzewski, 2010) w.r.t. \eqn{\kappa} for \eqn{x} is defined as
#' \deqn{V_n(x)=\max_{i=1,\dots,n}\{\min\{x_i, \kappa(i/n) \}\}}{V_n(x)=max{ min{\kappa(i/n), x_i} } for i=1,\dots,n}
#'
#' If \code{disable.check} is set to \code{FALSE}, then
#' eventual \code{NA} values are removed from the input vector.
#'
#' If a non-increasingly sorted vector is given as input (set \code{sorted.dec} to \code{TRUE})
#' the result is computed in linear time.
#'
#' @references Gagolewski M., Grzegorzewski P., S-Statistics and Their Basic Properties,
#' In: Borgelt C. et al (Eds.), Combining Soft Computing and Statistical Methods in Data Analysis, (SMPS 2010), Springer-Verlag, 2010, 281-288.
#'
#' @title S-statistic
#' @param x a vector of real numbers.
#' @param kappa an increasing function, \eqn{\kappa} (see Details), a so-called control function.
#' @param sorted.dec logical; \code{TRUE} if the vector has already been sorted non-increasingly; defaults \code{FALSE}.
#' @param disable.check logical; \code{TRUE} to disable some validity checks on the input vector; defaults \code{FALSE}.
#' @return The function returns a single number or NA if improper input has been given.
#' @examples
#' x <- rpareto2(25, 1.05, 1);
#' kappa <- function(x) { pmax(0,pmin(1,x))*25; }
#' Sstat(x, kappa, FALSE, TRUE);
#' @seealso \code{\link{index.h}}, \code{\link{index.g}}, \code{\link{index.rp}}, \code{\link{index.lp}}, \code{\link{Sstat2}}, \code{\link{psstat}}, \code{\link{dsstat}}
#' @export
Sstat <- function(x, kappa, sorted.dec=FALSE, disable.check=FALSE)
{
	if (!disable.check)
	{
		if (length(x) == 0) return(0);
		if (mode(x) != "numeric") return(NA);
		if (any(x < 0)) return(NA);
		x <- x[!is.na(x)];
	}

	n <- length(x);

	if (!sorted.dec)
		x <- sort(x, decreasing=TRUE);

	return (max(pmin(x,kappa((1:n)/n))));
}


