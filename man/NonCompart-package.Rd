\name{NonCompart-package}
\alias{NonCompart-package}
\alias{NonCompart}
\docType{package}
\title{Noncompartmental Analysis for Pharmacokinetic Data}
\description{
It conducts a noncompartmental analysis(NCA) as closely as possible to the most widely used commercial pharmacokinetic analysis software.
}

\details{
The main functions are
\preformatted{
NCA	     to perform NCA for many subjects.

IndiNCA  to perform NCA for one subject.
}
}
\author{
Kyun-Seop Bae <k@acr.kr>, Jee Eun Lee <JeeEun.Lee@fda.hhs.gov>
}

\references{
\enumerate{
\item Gabrielsson J, Weiner D. Pharmacokinetic and Pharmacodynamic Data Analysis - Concepts and Applications. 5th ed. 2016.
\item Shargel L, Yu A. Applied Biopharmaceutics and Pharmacokinetics. 7th ed. 2015.
\item Rowland M, Tozer TN. Clinical Pharmacokinetics and Pharmacodynamics - Concepts and Applications. 4th ed. 2011.
\item Gibaldi M, Perrier D. Pharmacokinetics. 2nd ed. revised and expanded. 1982.
}
}
\keyword{ package }
\keyword{ NCA }
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
