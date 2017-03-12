\name{plotPK}
\alias{plotPK}
\title{
Plot concentration vs. time curve for individuals and collectively.
}
\description{
Generates individual and superposed concentration vs. time curve and save it in pdf files.
}
\usage{
plotPK(concData, id, Time, conc, unitTime = "hr", unitConc = "ng/mL", trt = "", 
       fit = "Linear", dose = 0, adm = "Extravascular", dur = 0, outdir = "Output")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{concData}{name of data table containing time-concentration data of multiple subjects}
  \item{id}{column name for subject ID}
  \item{Time}{column name for the time}
  \item{conc}{column name for the concentration}
  \item{unitTime}{unit for the time}
  \item{unitConc}{unit for the concentration}
  \item{trt}{column name for the treatment code. This is useful for crossover study like bioequivalence trial.}
  \item{fit}{one of \code{"Linear"} or \code{"Log"} to indicate the way to calculate AUC}
  \item{dose}{administered dose. One should be careful for the unit. This can be a vector containing dose for each subject in order.}
  \item{adm}{one of \code{"Bolus"} or \code{"Infusion"} or \code{"Extravascular"} to indicate drug administration mode}
  \item{dur}{infusion duration for constant infusion, otherwise 0.  This can be a vector containing values for each subject in order.}
  \item{outdir}{name of the folder to be used for the output files}
}
\details{
This function generates plots for individual and summary concentration vs. time curve. This function calles \code{NCA()}.
}
\value{
This function saves pdf files and tiff files in the \code{outdir} folder.
}
\author{Jee Eun Lee}

\seealso{\code{\link{NCA}}}
\examples{
plotPK(Theoph, "Subject", "Time", "conc", unitTime="hr", unitConc="mg/L", dose=320)
plotPK(Indometh, "Subject", "time", "conc", unitTime="hr", unitConc="mg/L", adm="Bolus", dose=25)
}
\keyword{Plot}
