\name{sqlSwitchOrNULL}
\alias{sqlSwitchOrNULL}
\title{Switch-like nullable construct for use in an SQL query}
\usage{
  sqlSwitchOrNULL(value, search, replace)
}
\arguments{
  \item{value}{character vector to be processed.}

  \item{search}{character vector of length \code{n}.}

  \item{replace}{character vector of length \code{n}.}
}
\value{
  The function returns a character vector containing values
  from \code{replace} that correspond to elements of
  \code{search} and/or \code{"NULL"}s.
}
\description{
  Given a character vector, the function tries to match
  each of its values up in the array \code{search}. On
  success, corresponding value of \code{replace} will is
  returned. Otherwise, it outputs \code{"NULL"}.
}
\seealso{
  \code{\link{sqlNumericOrNULL}},
  \code{\link{sqlStringOrNULL}}, \code{\link{sqlTrim}}
}

