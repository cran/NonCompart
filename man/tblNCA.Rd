\name{tblNCA}
\alias{tblNCA}
\title{Table output NCA}
\description{Do multiple NCA and returns a result table. See \code{sNCA} for more detail i.e. \code{iAUC}}
\usage{
tblNCA(concData, key = "Subject", colTime = "Time", colConc = "conc", dose = 0, 
       adm = "Extravascular", dur = 0, doseUnit = "mg", timeUnit = "h", 
       concUnit = "ug/L", down = "Linear", R2ADJ = 0, MW = 0, SS = FALSE, 
       iAUC = "", excludeDelta = 1)
}
\arguments{
  \item{concData}{concentration data table}
  \item{key}{column names of concData to be shown in the output table}
  \item{colTime}{column name for time}
  \item{colConc}{column name for concentration}
  \item{dose}{administered dose}
  \item{adm}{one of \code{"Bolus"} or \code{"Infusion"} or \code{"Extravascular"} to indicate drug administration mode}
  \item{dur}{duration of infusion}
  \item{doseUnit}{unit of dose}
  \item{timeUnit}{unit of time}
  \item{concUnit}{unit of concentration}
  \item{down}{method to calculate AUC, \code{"Linear"} or \code{"Log"}}
  \item{R2ADJ}{Lowest threshold of adjusted R-square value to do manual slope determination}
  \item{MW}{molecular weight of drug}
  \item{SS}{if steady-state, this should be TRUE. AUCLST (AUClast) is used instead of AUCIFO (AUCinf) for the calculation of Vz (VZFO, VZO), CL (CLFO, CLO), and Vdss (VSSO).}
  \item{iAUC}{data.frame for interval AUC}
  \item{excludeDelta}{Improvement of R2ADJ larger than this value could exclude the last point. Default value 1 is for the compatibility with other software.}
}
\details{
  Tabular output of NCA with many subjects. Author recommends to use \code{excludeDelta} option with about 0.3.
}
\value{Basically same with \code{\link{sNCA}}}
\author{Kyun-Seop Bae <k@acr.kr>}
\seealso{
\code{\link{help}}, \code{\link{sNCA}}
}
\examples{
tblNCA(Theoph, key="Subject", dose=320, concUnit="mg/L")
tblNCA(Indometh, key="Subject", colTime="time", colConc="conc", dose=25, 
       adm="Infusion", dur=0.5, concUnit="mg/L")
}
\keyword{Output Form}
