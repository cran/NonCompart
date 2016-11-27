NCA <-
function(Data, colSubj, colTime, colConc, colTrt, Method="Linear", Dose=0, AdmMode="Extravascular", TimeInfusion=0, Report="Table", iAUC)
{
  if (!is.numeric(Dose) | !is.numeric(TimeInfusion) | !is.character(AdmMode) | !is.character(Method)) stop("Bad Input!")
  if (AdmMode == "Infusion" & !(TimeInfusion > 0)) stop("Infusion mode should have TimeInfusion larger than 0!")

  colOrd = paste0(AdmMode, "Default")
  ColName00 = RptCfg[RptCfg[,colOrd] > 0, c("PPTESTCD", colOrd)]
  ColName0 = ColName00[order(ColName00[, colOrd]), "PPTESTCD"] ;

  if (!(Dose > 0)) ColName0= setdiff(ColName0, c("CMAXD", "AUCIFOD", "AUCIFPD"))

  if (!missing(iAUC)) {
    ColName0 = union(ColName0, as.character(iAUC[,"Name"]))
  }

  SUBJIDs = unique(as.character(Data[,colSubj]))
  nSUBJID = length(SUBJIDs)
  if (missing(colTrt)) {
    Res0 = data.frame(SUBJID=character(), stringsAsFactors=FALSE)
    if (Report == "Table") {
      Result = data.frame()
    } else {
      Result = vector()
    }
    for (i in 1:nSUBJID) {
      cSUBJID = SUBJIDs[i]
      Dat = Data[Data[,colSubj]==cSUBJID,]
      if (nrow(Dat) > 0) {
        x = Dat[,colTime]
        y = Dat[,colConc]
        Res0 = rbind(Res0, data.frame(cSUBJID, stringsAsFactors=FALSE))
        if (!missing(iAUC)) {
          cResult = IndiNCA(x, y, Method=Method, Dose=Dose, AdmMode=AdmMode, TimeInfusion=TimeInfusion, RetNames=ColName0, Report=Report, iAUC=iAUC)
        } else {
          cResult = IndiNCA(x, y, Method=Method, Dose=Dose, AdmMode=AdmMode, TimeInfusion=TimeInfusion, RetNames=ColName0, Report=Report)
        }
        if (Report == "Table") {
          Result = rbind(Result, cResult)
        } else {
          Result = c(Result, "NCA REPORT", paste0("Subject=", cSUBJID), "", cResult, "", "")
        }
      }
    }
    if (Report == "Table") {
      Result = cbind(Res0, Result)
      colnames(Result) = c(colSubj, ColName0)
    }
  } else {
    TRTs = sort(unique(as.character(Data[,colTrt])))
    nTRT = length(TRTs)
    Res0 = data.frame(SUBJID=character(), TRT=character(), stringsAsFactors=FALSE)
    if (Report == "Table") {
      Result = data.frame()
    } else {
      Result = vector()
    }
    for (i in 1:nSUBJID) {
      for (j in 1:nTRT) {
        cSUBJID = SUBJIDs[i]
        cTRT = TRTs[j]
        Dat = Data[Data[,colSubj]==cSUBJID & Data[,colTrt]==cTRT,]
        if (nrow(Dat) > 0) {
          x = Dat[,colTime]
          y = Dat[,colConc]
          Res0 = rbind(Res0, data.frame(cSUBJID, cTRT, stringsAsFactors=FALSE))
          if (!missing(iAUC)) {
            cResult = IndiNCA(x, y, Method=Method, Dose=Dose, AdmMode=AdmMode, TimeInfusion=TimeInfusion, RetNames=ColName0, Report=Report, iAUC=iAUC)
          } else {
            cResult = IndiNCA(x, y, Method=Method, Dose=Dose, AdmMode=AdmMode, TimeInfusion=TimeInfusion, RetNames=ColName0, Report=Report)
          }
          if (Report == "Table") {
            Result = rbind(Result, cResult)
          } else {
            Result = c(Result, "NCA REPORT", paste0("Subject=", cSUBJID), paste0("Treatment=", cTRT), cResult, "", "")
          }
        }
      }
    }
    if (Report == "Table") {
      Result = cbind(Res0, Result)
      colnames(Result) = c(colSubj, colTrt, ColName0)
    }
  }

  return(Result)
}
