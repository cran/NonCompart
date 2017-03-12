\name{loadEXPC}
\alias{loadEXPC}
\title{Load EX and PC domain files in folders}
\description{This loads and returns EX and PC domain files in the specified folders}
\usage{loadEXPC(folders)}
\arguments{
  \item{folders}{folders where to find EX and PC domain files}
}
\details{
  This reads EX and PC domain files in the specified folder.
  This calls \code{readEX} and \code{readPC} functions.
}
\value{
  \item{EX}{combined EX domain data}
  \item{PC}{combined PC doamin data}
}
\author{Kyun-Seop Bae <k@acr.kr>}

\seealso{
\code{\link{help}}, \code{\link{readEX}}, \code{\link{readPC}}
}
\keyword{ XPT }
