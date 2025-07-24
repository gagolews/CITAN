## This file is part of the CITAN package for R
##
## Copyright 2011-2025 Marek Gagolewski
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


#' \pkg{CITAN} is a library of functions useful in
#' quantitative research in the field of scientometrics.
#'
#' The package is deprecated, see \pkg{agop} instead.
#'
#'
#' For the complete list of functions, call \code{library(help="CITAN")}.
#'
#'
#' @references
#' Gagolewski M. (2011). Bibliometric Impact Assessment with R and the CITAN Package,
#'    Journal of Informetrics 5(4), 678-692.
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
#' @importFrom methods is
#' @keywords internal
"_PACKAGE"
