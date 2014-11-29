\docType{data}
\name{Scopus_ASJC}
\alias{Scopus_ASJC}
\title{Scopus ASJC (All Science. Journals Classification) classification codes}
\description{
  List of Elsevier's \emph{SciVerse Scopus} ASJC (All
  Science. Journals Classification) source classification
  codes.
}
\details{
  Last update: October 2011. The data file is based on the
  official and publicly available (no permission needed as
  stated by Elsevier) Scopus list of covered titles, see
  \url{http://www.info.sciverse.com/documents/files/scopus-training/resourcelibrary/xls/title_list.xls}.

  It consists of 334 ASJC 4-digit integer codes (column
  \code{ASJC}) together with their group identifiers
  (column \code{ASJC_Parent}) and descriptions (column
  \code{Description}).

  ASJC codes are used to classify Scopus sources (see
  \code{\link{Scopus_SourceList}}).
}
\author{
  Marek Gagolewski \email{gagolews@ibspan.waw.pl}
}
\references{
  \url{http://www.info.sciverse.com/scopus/scopus-in-detail/facts/}
}
\seealso{
  \code{\link{Scopus_SourceList}},
  \code{\link{Scopus_ReadCSV}},
  \code{\link{Scopus_ImportSources}}
}
\keyword{ASJC,}
\keyword{datasets}
\keyword{journal}
\keyword{Scopus,}

