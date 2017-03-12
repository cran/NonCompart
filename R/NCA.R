NCA = function(concData, id, Time, conc, trt="", fit="Linear", dose=0, adm="Extravascular", dur=0, report="Table", iAUC="", uTime="h", uConc="ug/L", uDose="mg")
{
  if (!is.numeric(dose) | !is.numeric(dur) | !is.character(adm) | !is.character(fit)) stop("Bad Input!")

  colOrd = paste0(adm, "Default")
  ColName00 = RptCfg[RptCfg[,colOrd] > 0, c("PPTESTCD", colOrd)]
  ColName0 = ColName00[order(ColName00[, colOrd]), "PPTESTCD"] ;

  if (!(max(dose) > 0)) ColName0= setdiff(ColName0, c("CMAXD", "AUCIFOD", "AUCIFPD"))

  if (!missing(iAUC)) {
    ColName0 = union(ColName0, as.character(iAUC[,"Name"]))
  }

  SUBJIDs = unique(as.character(concData[,id]))
  nSUBJID = length(SUBJIDs)
  
  if (trt == "") {
    TRTs = ""
    nTRT = 1
  } else {
    TRTs = sort(unique(as.character(concData[,trt])))
    nTRT = length(TRTs)
  }

  if (length(dose) > 1 & length(dose) != nSUBJID*nTRT) stop("dose should be fixed or given for each subject!")
  if (length(dur) > 1 & length(dur) != nSUBJID*nTRT) stop("dur should be fixed or given for each subject!")

  if (trt == "") {
    Res0 = data.frame(SUBJID=character(), stringsAsFactors=FALSE)
    if (report == "Table") {
      Result = data.frame()
    } else {
      Result = vector()
    }
    for (i in 1:nSUBJID) {
      cSUBJID = SUBJIDs[i]
      Dat = concData[concData[,id]==cSUBJID,]
      if (nrow(Dat) > 0) {
        x = as.numeric(Dat[,Time])
        y = as.numeric(Dat[,conc])
        if (length(dose) > 1) {
          cDose = dose[i]
        } else {
          cDose = dose 
        }
        if (length(dur) > 1) {
          cTimeInfusion = dur[i]
        } else {
          cTimeInfusion = dur
        }
        if (adm == "Infusion" & !(cTimeInfusion > 0)) stop("Infusion mode should have dur larger than 0!")

        Res0 = rbind(Res0, data.frame(cSUBJID, stringsAsFactors=FALSE))
        cResult = IndiNCA(x, y, fit=fit, dose=cDose, adm=adm, dur=cTimeInfusion, report=report, iAUC=iAUC, uTime=uTime, uConc=uConc, uDose=uDose)
        if (report == "Table") {
          Result = rbind(Result, cResult)
        } else {
          Result = c(Result, "NCA REPORT", paste0("Subject=", cSUBJID), "", cResult, "", "")
        }
      }
    }
    if (report == "Table") {
      Result = cbind(Res0, Result)
      colnames(Result) = c(id, ColName0)
    }
  } else {
    TRTs = sort(unique(as.character(concData[,trt])))
    nTRT = length(TRTs)
    Res0 = data.frame(SUBJID=character(), TRT=character(), stringsAsFactors=FALSE)
    if (report == "Table") {
      Result = data.frame()
    } else {
      Result = vector()
    }
    for (i in 1:nSUBJID) {
      for (j in 1:nTRT) {
        cSUBJID = SUBJIDs[i]
        cTRT = TRTs[j]
        Dat = concData[concData[,id]==cSUBJID & concData[,trt]==cTRT,]
        if (nrow(Dat) > 0) {
          x = as.numeric(Dat[,Time])
          y = as.numeric(Dat[,conc])
          Res0 = rbind(Res0, data.frame(cSUBJID, cTRT, stringsAsFactors=FALSE))
          cResult = IndiNCA(x, y, fit=fit, dose=dose, adm=adm, dur=dur, report=report, iAUC=iAUC, uTime=uTime, uConc=uConc, uDose=uDose)
          if (report == "Table") {
            Result = rbind(Result, cResult)
          } else {
            Result = c(Result, "NCA REPORT", paste0("Subject=", cSUBJID), paste0("Treatment=", cTRT), cResult, "", "")
          }
        }
      }
    }
    if (report == "Table") {
      Result = cbind(Res0, Result)
      colnames(Result) = c(id, trt, ColName0)
    }
  }

  return(Result)
}
