foreNCA = function(NCAres="", PPTESTCD="", PCTESTCD="", title="", ...)
{ 
  cPara = PPTESTCD
  MOLEs = sort(as.character(unique(NCAres[,"PCTESTCD"]))) ; nMOLE = length(MOLEs)
  STUDYs = sort(as.character(unique(NCAres[,"STUDYID"]))) ; nSTUDY = length(STUDYs)
  DOSEs = sort(unique(NCAres[,"Dose"])) ; nDOSE = length(DOSEs)

  ForRes = c("Molecule", "Study", "Dose", "Mean", "SD")
  for (i in 1:nMOLE) {
    cMOLE = MOLEs[i]
    for (j in 1:nSTUDY) {
      cSTUDY = STUDYs[j]
      for (k in 1:nDOSE) {
        cDOSE = DOSEs[k]
        cData = NCAres[NCAres[,"PCTESTCD"]==cMOLE & NCAres[,"STUDYID"]==cSTUDY & NCAres[,"Dose"] == cDOSE, cPara]
        if (length(cData) > 0) {
          cMean = mean(cData, na.rm=TRUE)
          cSD = sd(cData, na.rm=TRUE)
          ForRes = rbind(ForRes, c(cMOLE, cSTUDY, cDOSE, cMean, cSD))
        }
      }
    }
  }

  unitDose = unique(as.character(NCAres[,"UnitDose"]))
  unitConc = unique(as.character(NCAres[,"UnitConc"]))
  unitPara = as.character(Unit(code=PPTESTCD, timeUnit="h", concUnit=unitConc, doseUnit=unitDose)[1])

  For1 = ForRes[ForRes[,1] == PCTESTCD,]
  cMean = c(NA, NA, as.numeric(For1[,4]))
  cLL = c(NA, NA, as.numeric(For1[,4]) - 2*as.numeric(For1[,5]))
  cUL = c(NA, NA, as.numeric(For1[,4]) + 2*as.numeric(For1[,5]))
  For1[,4] = Round(as.numeric(For1[,4]),2)
  For1[,5] = Round(as.numeric(For1[,5]),2)
  For1 = rbind(ForRes[1,], c("", "", unitDose, unitPara, ""), For1)
  
  if (title == "") {
    Synonym = RptCfg[RptCfg[,"PPTESTCD"]==toupper(Trim(PPTESTCD)), "NCI"]
  } else {
    Synonym = title 
  }
  xlbl = paste0(Synonym, " (", unitPara,")")

#  require(forestplot)
  forestplot::forestplot(For1[,2:5], cMean, cLL, cUL, txt_gp=forestplot::fpTxtGp(ticks=grid::gpar(cex=1), xlab=grid::gpar(cex=1), title=grid::gpar(cex=2)), title=Synonym, xlab=xlbl, ...)
}
