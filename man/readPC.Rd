\name{readPC}
\alias{readPC}
\title{Read PC domain files}
\description{This reads PC domain files from the specified folders.}
\usage{
readPC(folders)
}
\arguments{
  \item{folders}{folders where to find PC doamin files}
}
\details{
This calls \code{combXPT} function.
This is called by \code{loadEXPC} function.
}
\value{This returns combined table of PC doamin.}
\author{Kyun-Seop Bae <k@acr.kr>}
\seealso{
\code{\link{help}}, \code{\link{combXPT}}, \code{\link{loadEXPC}}
}
\keyword{ XPT }