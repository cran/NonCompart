\name{foreNCA}
\alias{foreNCA}
\title{Forest plot to compare NCA results}
\description{This function compares NCA results usually from rNCA function}
\usage{
foreNCA(NCAres = "", PPTESTCD = "", PCTESTCD = "", title = "", ...)
}
\arguments{
  \item{NCAres}{NCA results from rNCA function}
  \item{PPTESTCD}{CDISC SDTM PP domain Test Code to coompare}
  \item{PCTESTCD}{Molecular species to compare specified in PCTESTCD of CDISC SDTM PC domain}
  \item{title}{Title of the plot}
  \item{\dots}{further aguments to pass to the forestplot function}
}
\details{This functio calls \code{forestplot} in forest package.}
\value{Currently, this just plots.}
\author{Kyun-Seop Bae <k@acr.kr>}
\seealso{
\code{\link{help}}, \code{\link{rNCA}}
}
\keyword{Forest Plot}% use one of  RShowDoc("KEYWORDS")
