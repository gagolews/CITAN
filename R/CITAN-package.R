## This file is part of the CITAN library.
##
## Copyright 2011-2012 Marek Gagolewski
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


#' \pkg{CITAN} is a library of functions useful in --- but not limited to ---
#' quantitative research in the field of scientometrics.
#' It contains various tools for preprocessing
#' bibliographic data retrieved from e.g. Elsevier's \emph{SciVerse Scopus} and
#' calculating bibliometric impact of individuals.
#' Also, some functions dealing with Pareto-Type II (GPD)
#' and Discretized Pareto-Type II statistical models
#' are included (e.g. Zhang-Stephens and MLE estimators,
#' goodness-of-fit and two-sample tests, confidence intervals
#' for the theoretical Hirsch index etc.).
#' They may be used to describe and analyze many phenomena encountered
#' in the social sciences.
#'
#'
#' Fair and objective assessment methods of individual scientists
#' had become the focus of scientometricians' attention since the
#' very beginning of their discipline. A quantitative expression
#' of some publication-citation process'
#' characteristics is assumed to be a predictor of broadly conceived
#' scientific competence. It may be used e.g. in building decision support
#' systems for scientific quality control.
#'
#' The \eqn{h}-index, proposed by J.E. Hirsch (2005)
#' is among the most popular scientific impact indicators.
#' An author who has published \eqn{n}
#' papers has the Hirsch index equal to \eqn{H}, if each of his \eqn{H}
#' publications were cited at least \eqn{H} times, and each of the
#' remaining \eqn{n-H} items were cited no more than \eqn{H} times. This simple
#' bibliometric tool quickly received much attention in the academic community
#' and started to be a subject of intensive
#' research. It was noted that, contrary to earlier approaches,
#' i.e. publication count, citation count, etc.,
#' this measure  concerns both productivity and impact of an
#' individual.
#' \cr
#'
#' In a broader perspective, this issue is a special case
#' of the so-called \dfn{Producer Assessment Problem}
#' (PAP; see Gagolewski, Grzegorzewski, 2010b).
#'
#' Consider a \emph{producer} (e.g. a writer, scientist, artist,
#' craftsman) and a nonempty set of his \emph{products} (e.g. books,
#' papers, works, goods). Suppose that each product is given a
#' \emph{rating} (of quality, popularity, etc.) which is a single
#' number in \eqn{I=[a,b]}, where \eqn{a} denotes the lowest admissible
#' valuation. We typically choose \eqn{I=[0,\infty]} (an interval
#' in the extended real line).
#' Some instances of the PAP are listed below.
#'
#' \tabular{cllll}{
#'   \tab \strong{Producer}    \tab \strong{Products}   \tab \strong{Rating method} \tab \strong{Discipline}\cr
#' A \tab Scientist            \tab Scientific articles \tab Number of citations    \tab Scientometrics\cr
#' B \tab Scientific institute \tab Scientists          \tab The h-index            \tab Scientometrics\cr
#' C \tab Web server           \tab Web pages           \tab Number of in-links     \tab Webometrics\cr
#' D \tab Artist               \tab Paintings           \tab Auction price          \tab Auctions\cr
#' E \tab Billboard company    \tab Advertisements      \tab Sale results           \tab Marketing\cr
#' }
#'
#' Each possible state of producer's activity can therefore be represented by a point
#' \eqn{x\in I^n} for some \eqn{n}. Our aim is thus to construct
#' and analyze --- both theoretically and empirically ---
#' aggregation operators (cf. Grabisch et al, 2009) which can be used for rating
#' producers. A family of such functions should take  the two
#' following aspects of producer's quality into account:
#' \itemize{
#'     \item the ability to make highly-rated products,
#'     \item overall productivity, \eqn{n}.
#' }
#' For some more formal considerations please refer to (Gagolewski, Grzegorzewski, 2011).
#' \cr\cr
#'
#'
#'
#'
#' The \pkg{CITAN} package consists of four types of tools.
#'
#' \strong{(1)}
#' Given a numeric vector, the first class of functions
#' \bold{computes the values of certain impact functions}.
#' Among them we have:
#' \enumerate{
#' \item Hirsch's \eqn{h}-index (Hirsch, 2005; see \code{\link{index.h}}),
#' \item Egghe's \eqn{g}-index (Egghe, 2006; see \code{\link{index.g}}),
#' \item the \eqn{r_p} and \eqn{l_p} indices
#'       (Gagolewski, Grzegorzewski, 2009; Gagolewski, Debski, Nowakiewicz, 2009;
#'       see \code{\link{index.rp}} and \code{\link{index.lp}}), which
#'       generalize the \eqn{h}-index, the \eqn{w}-index (Woeginger, 2008), and
#'       the MAXPROX-index (Kosmulski, 2007),
#' \item S-statistics (Gagolewski, Grzegorzewski, 2010a, 2011;
#'       see \code{\link{Sstat}} and \code{\link{Sstat2}}),
#'       which generalize the OWMax operators (Dubois et al, 1988)
#'       and the \eqn{h}- and \eqn{r_\infty}-indices.
#' }
#'
#' \strong{(2)}
#' To \bold{preprocess and analyze bibliometric data} (cf. Gagolewski, 2011) retrieved
#' from e.g.  Elsevier's \emph{SciVerse Scopus}
#' we need the \pkg{RSQLite} package. It is an interface to the free
#' SQLite DataBase Management System (see \url{http://www.sqlite.org/}).
#' All data is stored in a so-called Local Bibliometric Storage (\acronym{LBS}),
#' created with the \code{\link{lbsCreate}} function.
#'
#' The data frames \code{\link{Scopus_ASJC}} and \code{\link{Scopus_SourceList}}
#' contain various information on current source coverage of SciVerse Scopus.
#' They may be needed during the creation of the LBS and \code{\link{lbsCreate}}
#' for more details.
#' \emph{License information: this data are publicly available
#'       and hence no special permission is needed to redistribute them
#'       (information from Elsevier).}
#'
#' \pkg{CITAN} is able to import publication data from Scopus CSV files
#' (saved with settings "Output: complete format" or "Output: Citations only",
#' see \code{\link{Scopus_ReadCSV}}). Note that the output limit in Scopus
#' is 2000 entries per file. Therefore, to perform
#' bibliometric research we often need to divide the query results
#' into many parts. \pkg{CITAN} is able to merge them back even if
#' records are repeated.
#'
#' The data may be accessed via functions from the \pkg{DBI} interface.
#' However, some typical tasks may be automated using
#' e.g. \code{\link{lbsDescriptiveStats}} (basic description of the whole sample
#' or its subsets, called \sQuote{Surveys}),
#' \code{\link{lbsGetCitations}} (gather citation sequences selected
#' authors), and \code{\link{lbsAssess}} (mass-compute impact functions'
#' values for given citation sequences).
#'
#' There are also some helpful functions (in **EXPERIMENTAL** stage) which use
#' the \pkg{RGtk2} library (see Lawrence, Lang, 2010)
#' to display some suggestions on which documents or authors should be
#' merged, see \code{\link{lbsFindDuplicateTitles}} and
#' \code{\link{lbsFindDuplicateAuthors}}.
#'
#' \strong{(3)}
#' Additionally, a set of \bold{functions dealing with stochastic aspects
#' of S-statistics (generalized OWMax operators), the \eqn{h}-index} and the Pareto type-II statistical models
#' is included (Gagolewski, Grzegorzewski, 2010a). We have the following.
#' \itemize{
#' \item Functions that work for any continuous distribution
#' (cf. Gagolewski, Grzegorzewski, 2010a):
#' \enumerate{
#'    \item \code{\link{psstat}}, \code{\link{dsstat}} for computing
#'         the distribution of S-statistics generated by some control function,
#'    \item \code{\link{phirsch}}, \code{\link{dhirsch}} for computing
#'         the distribution of the Hirsch index,
#'    \item \code{\link{rho.get}} for computing the so-called \eqn{\rho}-index
#'         (\eqn{\rho_\kappa}), which is a particular location characteristic
#'         of a given probability distribution depending on
#'         a control function \eqn{\kappa}.
#' }
#'
#' \item Tools for the Pareto-type II family:
#' \enumerate{
#'    \item \code{\link{ppareto2}}, \code{\link{dpareto2}},
#'        \code{\link{qpareto2}}, \code{\link{rpareto2}} for general functions
#'        dealing with the Pareto distribution of the second kind,
#'        including the c.d.f., p.d.f, quantiles and random deviates,
#'    \item \code{\link{pareto2.phirsch}}, \code{\link{pareto2.dhirsch}} for
#'        computing the distribution of the Hirsch index (much faster than
#'        the above general versions),
#'    \item \code{\link{pareto2.htest}} --- two-sample \eqn{h}-test for
#'        equality of shape parameters based on the difference
#'        of \eqn{h}-indices,
#'    \item \code{\link{pareto2.htest.approx}} --- two-sample asymptotic
#'        (approximate) \eqn{h}-test,
#'    \item \code{\link{pareto2.ftest}} --- two-sample exact F-test for
#'        equality of shape parameters,
#'    \item \code{\link{pareto2.zsestimate}} --- estimation of parameters
#'        using the Bayesian method (MMSE) developed by
#'        Zhang and Stevens (2009),
#'    \item \code{\link{pareto2.mlekestimate}},
#'        \code{\link{pareto2.mleksestimate}} --- estimation of parameters
#'        using the MLE,
#'    \item \code{\link{discrpareto2.mlekestimate}},
#'        \code{\link{discrpareto2.mleksestimate}} --- estimation of parameters
#'        of the Discretized Pareto-type II distribution using the MLE,
#'    \item \code{\link{pareto2.goftest}}, \code{\link{discrpareto2.goftest}}
#'        --- goodness-of-fit tests,
#'    \item \code{\link{pareto2.confint.rho}},
#'        \code{\link{pareto2.confint.rho.approx}} --- exact and
#'        approximate (asymptotic) confidence intervals for
#'        the \eqn{\rho}-index basing on S-statistics,
#'    \item \code{\link{pareto2.confint.h}} --- exact confidence intervals
#'        for the theoretical \eqn{h}-index.
#' }
#' }
#'
#'
#' \strong{(4)}
#' Moreover, we have implemented some \bold{simple graphical methods}
#' than may be used to illustrate various aspects of data being analyzed,
#' see \code{\link{plot.citfun}}, \code{\link{curve.add.rp}},
#' and \code{\link{curve.add.lp}}.
#' \cr\cr
#'
#'
#'
#'
#' Please feel free to send any comments and suggestions (e.g.
#' to include some new bibliometric impact indices) to the author
#' (see also \url{http://www.ibspan.waw.pl/~gagolews}).
#'
#' For a complete list of functions, call \code{library(help="CITAN")}.
#' \cr\cr
#'
#' \bold{Keywords}: Hirsch's h-index, Egghe's g-index, L-statistics,
#' S-statistics, bibliometrics, scientometrics, informetrics,
#' webometrics, aggregation operators, arity-monotonicity,
#' impact functions, impact assessment.
#'
#' @name CITAN-package
#' @aliases CITAN
#' @docType package
#' @title CITation ANalysis toolpack
#' @author Marek Gagolewski \email{gagolews@@ibspan.waw.pl}
#' @references
#' CITAN homepage, \url{http://www.ibspan.waw.pl/~gagolews/CITAN/}\cr
#' GTK+ Project, \url{http://www.gtk.org/download.html}\cr
#' SQLite DBMS, \url{http://www.sqlite.org/}\cr
#' Dubois D., Prade H., Testemale C. (1988). Weighted fuzzy pattern matching,
#'    Fuzzy Sets and Systems 28, s. 313-331.\cr
#' Egghe L. (2006). Theory and practise of the g-index, Scientometrics 69(1),
#'    131-152.\cr
#' Gagolewski M., Grzegorzewski P. (2009). A geometric approach to the construction of
#'    scientific impact indices, Scientometrics 81(3), 617-634.\cr
#' Gagolewski M., Debski M., Nowakiewicz M. (2009). Efficient algorithms for computing
#'    ''geometric'' scientific impact indices, Research Report of Systems
#'    Research Institute, Polish Academy of Sciences RB/1/2009.\cr
#' Gagolewski M., Grzegorzewski P. (2010a). S-statistics and their basic properties,
#'    In: Borgelt C. et al (Eds.), Combining Soft Computing and Statistical
#'    Methods in Data Analysis, Springer-Verlag, 281-288.\cr
#' Gagolewski M., Grzegorzewski P. (2010b). Arity-monotonic extended aggregation
#'    operators, In: Hullermeier E., Kruse R., Hoffmann F. (Eds.),
#'    Information Processing and Management of Uncertainty in Knowledge-Based
#'    Systems, CCIS 80, Springer-Verlag, 693-702.\cr
#' Gagolewski M. (2011). Bibliometric Impact Assessment with R and the CITAN Package,
#' Journal of Informetrics 5(4), 678-692.\cr
#' Gagolewski M., Grzegorzewski P. (2011a). Axiomatic Characterizations of (quasi-)
#' L-statistics and S-statistics and the Producer Assessment Problem,
#; In: Galichet S., Montero J., Mauris G. (Eds.), Proc. 7th conf. European Society
#' for Fuzzy Logic and Technology (EUSFLAT/LFA 2011), Atlantic Press, 53-58.
#' Grabisch M., Pap E., Marichal J.-L., Mesiar R. (2009). Aggregation functions,
#'    Cambridge.\cr
#' Gagolewski M., Grzegorzewski P. (2011b). Possibilistic analysis of arity-monotonic
#'    aggregation operators and its relation to bibliometric impact assessment
#'    of individuals, International Journal of Approximate Reasoning 52(9), 1312-1324.\cr
#' Hirsch J.E. (2005). An index to quantify individual's scientific research output,
#'    Proceedings of the National Academy of Sciences 102(46),
#'    16569-16572.\cr
#' Kosmulski M. (2007). MAXPROD - A new index for assessment of the scientific output
#'    of an individual, and a comparison with the h-index, Cybermetrics 11(1).\cr
#' Lawrence M., Lang D.T. (2010). RGtk2: A graphical user interface toolkit for R,
#'    Journal of Statistical Software 37(8), 1-52.\cr
#' Woeginger G.J. (2008). An axiomatic characterization of the Hirsch-index,
#'    Mathematical Social Sciences 56(2), 224-232.\cr
#' Zhang J., Stevens M.A. (2009). A New and Efficient Estimation Method for the
#'    Generalized Pareto Distribution, Technometrics 51(3), 316-325.\cr
NA



.onLoad <- function(lib, pkg)
{
   library.dynam("CITAN", pkg, lib);
   packageStartupMessage("CITAN loaded. For more information please visit http://www.ibspan.waw.pl/~gagolews/CITAN/.");
}
