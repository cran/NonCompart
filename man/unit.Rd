\name{unit}
\alias{unit}
\title{
Returns unit of CDISC PP domain PPTESTCD
}
\description{
This function returns the unit for the PPTESTCD like CMAX, CMAXD, AUCLST, MRTEVLST.
}
\usage{
unit(code, uTime = "h", uConc = "ng/mL", uDose = "mg")
}
\arguments{
  \item{code}{PPTESTCD}
  \item{uTime}{unit of time}
  \item{uConc}{unit of concentration}
  \item{uDose}{unit of dose}
}
\details{
It is called by \code{IndiNCA}.
}
\value{
  \item{[1]}{unit character}
  \item{[2]}{conversion factor}
}
\author{
Kyun-Seop Bae <k@acr.kr>
}

\seealso{
\code{\link{IndiNCA}}
}
\examples{
unit("AUCLST")
unit("CMAXD")
unit(code="MRTEVLST", uTime="h", uConc="ug/L", uDose="mg")
}
\keyword{unit}