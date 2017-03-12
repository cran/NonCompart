\name{NCA}
\alias{NCA}
\title{Noncompartmental analysis for a dataset with multiple subjects}
\description{
conduct noncompartmental analysis for many subjects in a data table
}
\usage{NCA(concData, id, Time, conc, trt="", fit = "Linear", dose = 0, 
           adm = "Extravascular", dur = 0, report = "Table", iAUC = "", 
           uTime = "h", uConc = "ug/L", uDose = "mg")}
\arguments{
  \item{concData}{name of data table containing time-concentration data of multiple subjects}
  \item{id}{column name for subject ID}
  \item{Time}{column name for the time}
  \item{conc}{column name for the concentration}
  \item{trt}{column name for the treatment code. This is useful for crossover study like bioequivalence trial.}
  \item{fit}{one of \code{"Linear"} or \code{"Log"} to indicate the way to calculate AUC}
  \item{dose}{administered dose. One should be careful for the unit. This can be a vector containing dose for each subject in order.}
  \item{adm}{one of \code{"Bolus"} or \code{"Infusion"} or \code{"Extravascular"} to indicate drug administration mode}
  \item{dur}{infusion duration for constant infusion, otherwise 0.  This can be a vector containing values for each subject in order.}
  \item{report}{either of \code{"Table"} or \code{"Text"} to specify the type of return value}
  \item{iAUC}{data.frame with three columns, "Name", "Start", "End" to specify partial interval AUC}
  \item{uTime}{unit of time}
  \item{uConc}{unit of concentration}
  \item{uDose}{unit of dose}
}
\details{
  This function calls \code{IndiNCA} repeatedly to do NCA for each subject.
  If you specify \code{Report="Text"}, this function returns in free text format to be used in a report file.
}
\value{
  \item{CMAX}{maximum concentration, Cmax}
  \item{CMAXD}{dose normalized Cmax, CMAX / Dose, Cmax / Dose}
  \item{TMAX}{time of maximum concentration, Tmax}
  \item{TLAG}{time to observe the first non-zero concentration, for extravascular administration only}
  \item{CLST}{last positive concentration observed, Clast}
  \item{CLSTP}{last positive concentration predicted, Clast_pred}
  \item{TLST}{time of last positive concentration, Tlast}
  \item{LAMZHL}{half-life by lambda z, ln(2)/LAMZ}
  \item{LAMZ}{lambda_z negative of best fit terminal slope}
  \item{LAMZLL}{earliest time for LAMZ}
  \item{LAMZUL}{last time for LAMZ}
  \item{LAMZNPT}{number of points for LAMZ}
  \item{CORRXY}{correlation of log(concentration) and time}
  \item{R2}{R-squared}
  \item{R2ADJ}{R-squared adjusted}
  \item{C0}{back extrapolated concentration at time 0, for bolus intravascular administration only}
  \item{AUCLST}{AUC from 0 to TLST}
  \item{AUCALL}{AUC using all the given points, including trailing zero concentrations}
  \item{AUCIFO}{AUC infinity observed}
  \item{AUCIFOD}{AUCIFO / Dose}
  \item{AUCIFP}{AUC infinity predicted using CLSTP instead of CLST}
  \item{AUCIFPD}{AUCIFP / Dose}
  \item{AUCPEO}{AUC \% extrapolation observed}
  \item{AUCPEP}{AUC \% extrapolated for AUCIFP}
  \item{AUCPBEO}{AUC \% back extrapolation observed, for bolus IV administration only}
  \item{AUCPBEP}{AUC \% back extrapolation predicted with AUCIFP, for bolus IV administration only}
  \item{AUMCLST}{AUMC to the TLST}
  \item{AUMCIFO}{AUMC infinity observed using CLST}
  \item{AUMCIFP}{AUMC infinity determined by CLSTP}
  \item{AUMCPEO}{AUMC \% extrapolated observed}
  \item{AUMCPEP}{AUMC \% extrapolated predicted}
  \item{MRTIVLST}{mean residence time (MRT) to TLST, for intravascular administration}
  \item{MRTIVIFO}{mean residence time (MRT) infinity using CLST, for intravascular administration}
  \item{MRTIVIFP}{mean residence time (MRT) infinity using CLSTP, for intravascular administration}
  \item{MRTEVLST}{mean residence time (MRT) to TLST, for extravascular administration}
  \item{MRTEVIFO}{mean residence time (MRT) infinity using CLST, for extravascular administration}
  \item{MRTEVIFP}{mean residence time (MRT) infinity using CLSTP, for extravascular administration}
  \item{VZO}{volume of distribution determined by LAMZ and AUCIFO, for intravascular administration}
  \item{VZP}{volume of distribution determined by LAMZ and AUCIFP, for intravascular administration}
  \item{VZFO}{VZO for extravascular administration, VZO/F, F is bioavailability}
  \item{VZFP}{VZP for extravascular administration, VZP/F, F is bioavailability}
  \item{CLO}{clearance using AUCIFO, for intravascular administration}
  \item{CLP}{clearance using AUCIFP, for intravascular administration}
  \item{CLFO}{CLO for extravascular administration, CLO/F, F is bioavailability}
  \item{CLFP}{CLP for extravascular administration, CLP/F, F is bioavailability}
  \item{VSSO}{volume of distribution at steady state using CLST, for intravascular administration only}
  \item{VSSP}{volume of distribution at stead state using CLSTP, for intravascular administration only}
}
\references{
\enumerate{
\item Gabrielsson J, Weiner D. Pharmacokinetic and Pharmacodynamic Data Analysis - Concepts and Applications. 5th ed. 2016.
\item Shargel L, Yu A. Applied Biopharmaceutics and Pharmacokinetics. 7th ed. 2015. 
\item Rowland M, Tozer TN. Clinical Pharmacokinetics and Pharmacodynamics - Concepts and Applications. 4th ed. 2011.
\item Gibaldi M, Perrier D. Pharmacokinetics. 2nd ed. revised and expanded. 1982.
}
}
\author{Kyun-Seop Bae <k@acr.kr>}
\seealso{\code{\link{IndiNCA}}}
\examples{
# Theoph and Indometh data: dose in mg, conc in mg/L, time in h
NCA(Theoph, "Subject", "Time", "conc", dose=320, uConc="mg/L")
NCA(Indometh, "Subject", "time", "conc", dose=25, adm="Bolus", uConc="mg/L")

iAUC = data.frame(Name=c("AUC[0-12h]","AUC[0-24h]"), Start=c(0,0), End=c(12,24)) ; iAUC
NCA(Theoph, "Subject", "Time", "conc", dose=320, iAUC=iAUC, uConc="mg/L")
NCA(Indometh, "Subject", "time", "conc", dose=25, adm="Bolus", iAUC=iAUC, uConc="mg/L")

writeLines(NCA(Theoph, "Subject", "Time", "conc", dose=320, report="Text", uConc="mg/L"),
           "Theoph_Linear_CoreOutput.txt")
writeLines(NCA(Theoph, "Subject", "Time", "conc", dose=320, fit="Log", report="Text", 
           uConc="mg/L"), "Theoph_Log_CoreOutput.txt")
writeLines(NCA(Indometh, "Subject", "time", "conc", dose=25, adm="Bolus", report="Text", 
           uConc="mg/L"), "Indometh_Bolus_Linear_CoreOutput.txt")
writeLines(NCA(Indometh, "Subject", "time", "conc", dose=25, adm="Bolus", fit="Log",
           report="Text", uConc="mg/L"), "Indometh_Bolus_Log_CoreOutput.txt")
writeLines(NCA(Indometh, "Subject", "time", "conc", dose=25, adm="Infusion", dur=0.25,
           report="Text", uConc="mg/L"), "Indometh_Infusion_Linear_CoreOutput.txt")
writeLines(NCA(Indometh, "Subject", "time", "conc", dose=25, adm="Infusion", dur=0.25,
           fit="Log", report="Text", uConc="mg/L"), "Indometh_Infusion_Log_CoreOutput.txt")
}
\keyword{NCA}
