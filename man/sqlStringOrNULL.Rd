\name{sqlStringOrNULL}
\alias{sqlStringOrNULL}
\title{Get nullable character string values for use in an SQL query}
\usage{
  sqlStringOrNULL(value)
}
\arguments{
  \item{value}{character vector to be processed.}
}
\value{
  The function returns a character vector containing a
  processed input values and/or \code{"NULL"}s.
}
\description{
  Returns trimmed and single-quote-escaped (see
  \code{\link{sqlEscapeTrim}}) character strings or
  \code{"NULL"}s on empty input.
}
\seealso{
  \code{\link{sqlNumericOrNULL}},
  \code{\link{sqlSwitchOrNULL}},
  \code{\link{sqlEscapeTrim}}
}

