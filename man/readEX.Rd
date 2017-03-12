\name{readEX}
\alias{readEX}

\title{Read EX domain files}
\description{This reads EX domain files from the specified folders.}
\usage{
readEX(folders)
}
\arguments{
  \item{folders}{folders where to find EX doamin files}
}
\details{
This calls \code{combXPT} function.
This is called by \code{loadEXPC} function.
}
\value{This returns combined table of EX doamin.}
\author{Kyun-Seop Bae <k@acr.kr>}
\seealso{
\code{\link{help}}, \code{\link{combXPT}}, \code{\link{loadEXPC}}
}
\keyword{ XPT }
