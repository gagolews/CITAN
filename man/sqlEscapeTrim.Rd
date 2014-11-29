\name{sqlEscapeTrim}
\alias{sqlEscapeTrim}
\title{Escape character strings for use in an SQL query and trim white-spaces on both sides}
\usage{
  sqlEscapeTrim(str, useBytes = FALSE)
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
  Escapes given character strings for use in an SQL query
  and trims white-spaces on both sides.
}
\details{
  The SQL standard specifies that single-quotes in strings
  should be escaped by putting two single quotes in a row.
  This function repeats the quotes using
  \code{\link{gsub}}.
}
\seealso{
  \code{\link{gsub}}, \code{\link{sqlEscape}},
  \code{\link{sqlTrim}}
}

