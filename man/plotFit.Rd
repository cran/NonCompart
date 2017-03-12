\name{plotFit}
\alias{plotFit}
\title{
Plot best fit slope
}
\description{
Automatically select best fit slope for the given x(usually time) and log(y)(usually concentration) values.
}
\usage{
plotFit(concData, id, Time, conc, mol = "", adm = "Extravascular", ID = "", Mol = "")
}
\arguments{
  \item{concData}{name of data table containing time-concentration data of multiple subjects}
  \item{id}{column name for subject ID}
  \item{Time}{column name for the time}
  \item{conc}{column name for the concentration}
  \item{mol}{column name for molecular species}
  \item{adm}{one of \code{"Bolus"} or \code{"Infusion"} or \code{"Extravascular"} to indicate drug administration mode}
  \item{ID}{Subject ID for this plot}
  \item{Mol}{the name of molecular species to see}
}
\details{
Find the best fit slope then plot it.
Currently this function uses ordinary least square method(OLS) only.
This function calles \code{BestSlope} function.
}
\value{
  \item{R2}{R-squared}
  \item{R2ADJ}{adjusted R-squared}
  \item{LAMZNPT}{number of points used for slope}
  \item{LAMZ}{negative of slope, lambda_z}
  \item{b0}{intercept of regression line}
  \item{CORRXY}{correlation of log(y) and x}
  \item{LAMZLL}{earliest x for lambda_z}
  \item{LAMZUL}{last x for lambda_z}
  \item{CLSTP}{predicted y value at last point, predicted concentration for the last time point}
}

\author{
Jee Eun Lee
}


\seealso{
\code{\link{BestSlope}}
}
\examples{
plotFit(Theoph, "Subject", "Time", "conc", ID="1")
plotFit(Indometh, "Subject", "time", "conc", adm="Bolus", ID="1")
}
\keyword{Plot}
