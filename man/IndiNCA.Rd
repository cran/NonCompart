\name{IndiNCA}
\alias{IndiNCA}
\title{Noncompartmental Analysis for an Individual}
\description{It performs a noncompartmental analysis with one subject data}
\usage{IndiNCA(x, y, dose = 0, fit = "Linear", adm = "Extravascular", dur = 0,
               report = "Table", iAUC = "", uTime = "h", uConc = "ug/L", uDose = "mg")}
\arguments{
  \item{x}{vector values of independent variable, usually time}
  \item{y}{vector values of dependent variable, usually concentration}
  \item{dose}{administered dose for the subject}
  \item{fit}{either of \code{"Linear"} or \code{"Log"} to indicate the way to calculate AUC and AUMC}
  \item{adm}{one of \code{"Bolus"} or \code{"Infusion"} or \code{"Extravascular"} to indicate drug administration mode}
  \item{dur}{infusion duration for constant infusion, otherwise 0}
  \item{report}{either of \code{"Table"} or \code{"Text"} to specify the type of return value}
  \item{iAUC}{data.frame with three columns, "Name", "Start", "End" to specify the invervals for partial (interval) AUC}
  \item{uTime}{unit of time}
  \item{uConc}{unit of concentration}
  \item{uDose}{unit of dose}
}
\details{
  This performs a noncompartmental analysis for a subject. It returns practically the same result with the most popular commercial software.
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
\seealso{\code{\link{AUC}}, \code{\link{BestSlope}}}
\examples{
IndiNCA(Theoph[Theoph$Subject==1,"Time"], Theoph[Theoph$Subject==1, "conc"], dose=320, uConc="mg/L")
IndiNCA(Indometh[Indometh$Subject==1,"time"], Indometh[Indometh$Subject==1, "conc"], dose=25,
        adm="Bolus", uConc="mg/L")
IndiNCA(Indometh[Indometh$Subject==1,"time"], Indometh[Indometh$Subject==1, "conc"], dose=25,
        adm="Infusion", dur=0.25, uConc="mg/L")

IndiNCA(Theoph[Theoph$Subject==1,"Time"], Theoph[Theoph$Subject==1, "conc"], dose=320,
        report="Text", uConc="mg/L")
IndiNCA(Indometh[Indometh$Subject==1,"time"], Indometh[Indometh$Subject==1, "conc"], dose=25,
        adm="Bolus", report="Text", uConc="mg/L")
IndiNCA(Indometh[Indometh$Subject==1,"time"], Indometh[Indometh$Subject==1, "conc"], dose=25,
        adm="Infusion", dur=0.25, report="Text", uConc="mg/L")

iAUC = data.frame(Name=c("AUC[0-12h]","AUC[0-24h]"), Start=c(0,0), End=c(12,24)) ; iAUC
IndiNCA(Theoph[Theoph$Subject==1,"Time"], Theoph[Theoph$Subject==1, "conc"], dose=320,
        iAUC=iAUC, uConc="mg/L")
IndiNCA(Indometh[Indometh$Subject==1,"time"], Indometh[Indometh$Subject==1, "conc"], dose=25,
        adm="Bolus", iAUC=iAUC, uConc="mg/L")
IndiNCA(Indometh[Indometh$Subject==1,"time"], Indometh[Indometh$Subject==1, "conc"], dose=25,
        adm="Infusion", dur=0.25, iAUC=iAUC, uConc="mg/L")
}
\keyword{NCA}
\keyword{noncompartmenal analysis}
