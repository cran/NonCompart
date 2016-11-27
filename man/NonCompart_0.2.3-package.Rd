\name{NonCompart-package}
\alias{NonCompart-package}
\alias{NonCompart}
\docType{package}
\title{Noncompartmental Analysis for Pharmacokinetics}
\description{
Conduct noncompartmental analysis(NCA) as close as possible to the most widely used commercial pharmacokinetic analysis software.
}
\details{
The main functions are
\preformatted{
NCA	     to perform NCA for many subjects.

IndiNCA  to perform NCA for one subject.
}
}
\author{{Kyun-Seop Bae <k@acr.kr>}
Maintainer:{Kyun-Seop Bae <k@acr.kr>}
}
\references{
Rowland M, Tozer TN. Clinical Pharmacokinetics and Pharmacodynamics - Concepts and Applications. 4e. 2011.
}
\keyword{ package }
\keyword{ NCA }
\examples{
# Theoph and Indometh data: dose in mg, conc in mg/L, time in h
NCA(Theoph, "Subject", "Time", "conc", Dose=320)
NCA(Indometh, "Subject", "time", "conc", Dose=25, AdmMode="Bolus")

iAUC = data.frame(Name=c("AUC[0-12h]","AUC[0-24h]"), Start=c(0,0), End=c(12,24)) ; iAUC
NCA(Theoph, "Subject", "Time", "conc", Dose=320, iAUC=iAUC)
NCA(Indometh, "Subject", "time", "conc", Dose=25, AdmMode="Bolus", iAUC=iAUC)

writeLines(NCA(Theoph, "Subject", "Time", "conc", Dose=320, Report="Text"),
           "Theoph_Linear_CoreOutput.txt")
writeLines(NCA(Theoph, "Subject", "Time", "conc", Dose=320, Method="Log", Report="Text"),
           "Theoph_Log_CoreOutput.txt")
writeLines(NCA(Indometh, "Subject", "time", "conc", Dose=25, AdmMode="Bolus", Report="Text"),
           "Indometh_Bolus_Linear_CoreOutput.txt")
writeLines(NCA(Indometh, "Subject", "time", "conc", Dose=25, AdmMode="Bolus", Method="Log",
           Report="Text"), "Indometh_Bolus_Log_CoreOutput.txt")
writeLines(NCA(Indometh, "Subject", "time", "conc", Dose=25, AdmMode="Infusion", TimeInfusion=0.25,
           Report="Text"), "Indometh_Infusion_Linear_CoreOutput.txt")
writeLines(NCA(Indometh, "Subject", "time", "conc", Dose=25, AdmMode="Infusion", TimeInfusion=0.25,
           Method="Log", Report="Text"), "Indometh_Infusion_Log_CoreOutput.txt")

IndiNCA(Theoph[Theoph$Subject==1,"Time"], Theoph[Theoph$Subject==1, "conc"], Dose=320)
IndiNCA(Indometh[Indometh$Subject==1,"time"], Indometh[Indometh$Subject==1, "conc"], Dose=25,
        AdmMode="Bolus")
IndiNCA(Indometh[Indometh$Subject==1,"time"], Indometh[Indometh$Subject==1, "conc"], Dose=25,
        AdmMode="Infusion", TimeInfusion=0.25)

IndiNCA(Theoph[Theoph$Subject==1,"Time"], Theoph[Theoph$Subject==1, "conc"], Dose=320,
        Report="Text")
IndiNCA(Indometh[Indometh$Subject==1,"time"], Indometh[Indometh$Subject==1, "conc"], Dose=25,
        AdmMode="Bolus", Report="Text")
IndiNCA(Indometh[Indometh$Subject==1,"time"], Indometh[Indometh$Subject==1, "conc"], Dose=25,
        AdmMode="Infusion", TimeInfusion=0.25, Report="Text")

iAUC = data.frame(Name=c("AUC[0-12h]","AUC[0-24h]"), Start=c(0,0), End=c(12,24)) ; iAUC
IndiNCA(Theoph[Theoph$Subject==1,"Time"], Theoph[Theoph$Subject==1, "conc"], Dose=320,
        iAUC=iAUC)
IndiNCA(Indometh[Indometh$Subject==1,"time"], Indometh[Indometh$Subject==1, "conc"], Dose=25,
        AdmMode="Bolus", iAUC=iAUC)
IndiNCA(Indometh[Indometh$Subject==1,"time"], Indometh[Indometh$Subject==1, "conc"], Dose=25,
        AdmMode="Infusion", TimeInfusion=0.25, iAUC=iAUC)
}
