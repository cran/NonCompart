\name{rNCA}
\alias{rNCA}
\title{Do NCA for review}
\description{This performs NCA from the CDISC EX and PC datasets.}
\usage{
rNCA(ex, pc, study = "", trt = "", id = "", analyte = "", 
     codeBQL = c("< 0", "NQ", "BQL", "BQoL", "<LOQ"), MinPoints = 5)
}
\arguments{
  \item{ex}{EX domain data, usually from the \code{loadEXPC}}
  \item{pc}{PC domain data, usually form the \code{loadEXPC}}
  \item{study}{vector of study names in EX and PC domain to do NCA}
  \item{trt}{vector of treatment names in EXTRT to do NCA}
  \item{id}{vector of subject IDs in USUBJID to do NCA}
  \item{analyte}{vector of molecular species in PCTESTCD to do NCA}
  \item{codeBQL}{symbols of below the quantitation limit}
  \item{MinPoints}{minimum number of sampling points for NCA}
}
\details{
This calls NCA0.
Results of this can be further processed by foreNCA to plot and compare between studies and dose groups.
}
\value{
This returns a table of NCA results
}
\author{Kyun-Seop Bae <k@acr.kr>}
\seealso{
\code{\link{help}}, \code{\link{NCA0}}, \code{\link{loadEXPC}}, \code{\link{foreNCA}} 
}
\keyword{NCA}
