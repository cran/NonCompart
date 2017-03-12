\name{NCA0}
\alias{NCA0}
\title{NCA of SDTM data for single subject}
\description{
This performs Noncompartmental Analysis(NCA) for only one subject from the CDISC EX and PC domain.
}
\usage{NCA0(EX0, PC0)}
\arguments{
  \item{EX0}{Data of one subject from EX domain}
  \item{PC0}{Data of one subject from PC domain}
}
\details{
This calls IndiNCA function.
This is called by rNCA function.
}
\value{
This returns NCA results vector.
}
\author{Kyun-Seop Bae <k@acr.kr>}

\seealso{
\code{\link{help}}, \code{\link{rNCA}}, \code{\link{IndiNCA}}
}
\keyword{ NCA }
