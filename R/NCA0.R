NCA0 = function(EX0, PC0)
{
  cDose = as.numeric(EX0["EXDOSE"])
  cUnitDose = EX0["EXDOSU"]
  cRoute = EX0["EXROUTE"]
  cStart = EX0["EXSTDTC"]
  cEnd = EX0["EXENDTC"]
  if (!is.na(cEnd) & cEnd != "NA" & cEnd != "") {
    if (cEnd > cStart) {
      admMethod = "Infusion"
      cDur = as.numeric(difftime(strptime(cEnd,"%Y-%m-%dT%H:%M:%S"), strptime(cStart, "%Y-%m-%dT%H:%M:%S"), units="hours"))
    } else if (toupper(cRoute) == "INTRAVENOUS") {
      admMethod = "Bolus"
      cDur = 0
    } else {
      admMethod = "Extravascular"
      cDur = 0
    }
  } else {
    admMethod = "Extravascular"
    cDur = 0
  }
  cUnitConc = unique(PC0[,"PCSTRESU"])
  TAD = as.numeric(difftime(strptime(PC0[,"PCDTC"],"%Y-%m-%dT%H:%M:%S"), strptime(cStart, "%Y-%m-%dT%H:%M:%S"), units="hours"))
  TAD[TAD < 0] = 0
  y = as.numeric(PC0[,"PCSTRESN"])
  Res0 = IndiNCA(TAD, y, dose=cDose, adm=admMethod, dur=cDur, uConc=cUnitConc, uDose=cUnitDose)
  tNames = names(Res0)
  Res = c(admMethod, cDose, cUnitDose, cDur, "h", cUnitConc, Res0)
  names(Res) = c("AdmMethod", "Dose", "UnitDose", "Dur", "UnitTime", "UnitConc", tNames)
  return(Res)
}
