\name{sqlEscape}
\alias{sqlEscape}
\title{Escape character strings for use in an SQL query}
\usage{
  sqlEscape(str, useBytes = FALSE)
}
\arguments{
  \item{str}{a character vector where matches are sought,
  or an object which can be coerced by \code{as.character}
  to a character vector.}

  \item{useBytes}{logical.  If \code{TRUE} the matching is
  done byte-by-byte rather than character-by-character.
  See man page for \code{\link{gsub}}.}
}
\value{
  See 'Value' for \code{\link{gsub}}.
}
\description{
  Escapes character strings for use in an SQL query.
}
\details{
  The SQL standard specifies that single-quotes in strings
  should be escaped by putting two single quotes in a row.
  This function repeats the quotes using
  \code{\link{gsub}}.
}
\seealso{
  \code{\link{gsub}}, \code{\link{sqlTrim}},
  \code{\link{sqlEscapeTrim}}
}

