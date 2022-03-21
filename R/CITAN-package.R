## This file is part of the CITAN package for R
##
## Copyright 2011-2022 Marek Gagolewski
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
#'
#' The package is deprecated, see \pkg{agop} instead.
#'
#'
#' For the complete list of functions, call \code{library(help="CITAN")}.
#' \cr\cr
#'
#' @name CITAN-package
#' @aliases CITAN
#' @docType package
#' @title CITation ANalysis toolpack
#' @author Marek Gagolewski
#' @references
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
#' for Fuzzy Logic and Technology (EUSFLAT/LFA 2011), Atlantis Press, 53-58.
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
#' Woeginger G.J. (2008). An axiomatic characterization of the Hirsch-index,
#'    Mathematical Social Sciences 56(2), 224-232.\cr
#' Zhang J., Stevens M.A. (2009). A New and Efficient Estimation Method for the
#'    Generalized Pareto Distribution, Technometrics 51(3), 316-325.\cr
#'
#' @importFrom stringi stri_trim_both
#' @importFrom stringi stri_replace_all_fixed
#' @importFrom RSQLite dbGetInfo
#' @importFrom RSQLite dbGetQuery
#' @importFrom RSQLite dbCommit
#' @importFrom DBI dbDriver
#' @importFrom RSQLite dbConnect
#' @importFrom DBI dbDisconnect
#' @importFrom RSQLite dbListTables
#' @importFrom grDevices as.graphicsAnnot
#' @importFrom grDevices dev.interactive
#' @importFrom grDevices devAskNewPage
#' @importFrom graphics barplot
#' @importFrom graphics boxplot
#' @importFrom graphics mtext
#' @importFrom graphics par
#' @importFrom graphics pie
#' @importFrom stats na.omit
#' @importFrom utils read.csv
invisible(NULL)
