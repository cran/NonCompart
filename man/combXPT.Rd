\name{combXPT}
\alias{combXPT}
\title{Combine XPT files}
\description{This function combines specified CDISC domain XPT files across the folders.}
\usage{combXPT(folders, domain)}
\arguments{
  \item{folders}{where to find specified CDISC domain XPT files}
  \item{domain}{domain XPT files to be comined across the folders}
}
\details{
You need to designate only one CDISC domain name.
You may specify one or more folders to find the domain XPT files.
}
\value{
  \item{XPT}{combined table}
}
\author{Kyun-Seop Bae <k@acr.kr>}
\seealso{\code{\link{help}},\code{\link{readEX}},\code{\link{readPC}}}
\keyword{XPT}
