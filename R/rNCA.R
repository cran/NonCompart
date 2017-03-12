rNCA = function(ex, pc, study="", trt="", id="", analyte="", codeBQL=c("< 0", "NQ", "BQL", "BQoL", "<LOQ"), MinPoints=5)
{
  for (i in 1:ncol(ex)) {
    ex[,i] = Trim(as.character(ex[,i]))
    ex[is.na(ex[,i]),i] = ""
  }
  for (i in 1:ncol(pc)) {
    pc[,i] = Trim(as.character(pc[,i]))
    pc[is.na(pc[,i]),i] = ""
  }
  Routes = unique(ex[,"EXROUTE"])
  if (length(Routes) > 1 & "INTRAVENOUS" %in% toupper(Routes)) stop("Intravenous and extravascular routes cannot be mixed!")

  ex[toupper(ex[,"EXTRT"])=="PLACEBO","EXDOSE"] = "0"
  ex = ex[!is.na(ex[,"EXSTDTC"]) & Trim(ex[,"EXSTDTC"]) != "",]

  pc[toupper(pc[,"PCSTRESC"]) %in% toupper(codeBQL), "PCSTRESN"] = "0"
  pc = pc[!is.na(pc[,"PCSTRESN"]) & pc[,"PCSTRESN"] != "",]
  pc[as.numeric(pc[,"PCSTRESN"]) < as.numeric(pc[,"PCLLOQ"]),"PCSTRESN"] = "0"

  if (study[1] == "") study = sort(unique(pc[,"STUDYID"]))
  if (trt[1] == "") trt = sort(unique(ex[,"EXTRT"]))
  trt = setdiff(toupper(trt), "PLACEBO")
  if (id[1] == "") id = sort(unique(pc[,"USUBJID"]))
  if (analyte[1] == "") analyte = sort(unique(pc[,"PCTESTCD"]))
  nAnal = length(analyte)

  EX1 = ex[(ex[,"STUDYID"] %in% study) & (toupper(ex[,"EXTRT"]) %in% trt) & (ex[,"USUBJID"] %in% id) & nchar(ex[,"EXSTDTC"])==19 & nchar(ex[,"EXENDTC"])==19,]
  id = sort(unique(EX1[,"USUBJID"]))
  PC1 = pc[pc[,"STUDYID"] %in% study & pc[,"USUBJID"] %in% id & toupper(pc[,"PCTESTCD"]) %in% toupper(analyte) & toupper(substr(pc[,"PCSPEC"],1,6)) == "PLASMA" & nchar(pc[,"PCDTC"])==19,]

  if ("PCTPTNUM" %in% colnames(PC1)) {
    colPC = c("PCTESTCD", "PCDTC", "PCSTRESN", "PCSTRESU", "PCTPTNUM")
  } else {
    colPC = c("PCTESTCD", "PCDTC", "PCSTRESN", "PCSTRESU")
  }
          
  IDs = sort(unique(PC1[,"USUBJID"]))
  nID = length(IDs)
  if (nID == 0) stop("No subject with the given conditon.")

  Res = vector()
  for (i in 1:nID) {
    cID = IDs[i]
    EXi = unique(EX1[EX1[,"USUBJID"]==cID, c("STUDYID", "USUBJID", "EXTRT", "EXDOSE", "EXDOSU", "EXROUTE", "EXSTDTC", "EXENDTC")])
    nEXi = nrow(EXi)
    PCi = PC1[PC1[,"USUBJID"]==cID,]
    if (nEXi == 1) {
      cStudy = EXi[1,"STUDYID"]
      cSubj = EXi[1, "USUBJID"]
      cTrt = EXi[1, "EXTRT"]
      cRefTime = EXi[1, "EXSTDTC"]
      PCy = PCi[,c("PCTESTCD", "PCDTC", "PCSTRESN", "PCSTRESU")]
      if (nAnal == 1) {
        if (sum(as.numeric(PCy[,"PCSTRESN"]) > 0) > MinPoints) {
#print(paste(i))
          Res0 = NCA0(EXi[1,], PCy)
          tName = names(Res0)
          Res1 = c(cStudy, cSubj, cTrt, analyte[1], cRefTime, Res0)
          names(Res1) = c("STUDYID", "USUBJID", "EXTRT", "PCTESTCD", "PCRFTDTC", tName)
          Res = rbind(Res, Res1)
        }
      } else {
        for (j in 1:nAnal) {
          cAnalyte = analyte[j]
          PCy2 = PCy[PCy[,"PCTESTCD"] == cAnalyte,]
          if (sum(as.numeric(PCy2[,"PCSTRESN"]) > 0) > MinPoints) {
#print(paste(i,j))
            Res0 = NCA0(EXi[1,], PCy2)
            tName = names(Res0)
            Res1 = c(cStudy, cSubj, cTrt, cAnalyte, cRefTime, Res0)
            names(Res1) = c("STUDYID", "USUBJID", "EXTRT", "PCTESTCD", "PCRFTDTC", tName)
            Res = rbind(Res, Res1)
          }
        }
      }
    } else {
      for (j in 1:nEXi) {
        cStudy = EXi[j,"STUDYID"]
        cSubj = EXi[j, "USUBJID"]
        cTrt = EXi[j, "EXTRT"]
        cRefTime = EXi[j, "EXSTDTC"]
        if (j < nEXi) {
          if (nAnal == 1) {
            cEXSTDTC = EXi[j, "EXSTDTC"]
            PCy = PCi[PCi[,"PCDTC"] > EXi[j, "EXSTDTC"] & PCi["PCDTC"] < EXi[j+1, "EXSTDTC"], colPC]
            PCy0 = PCi[substr(PCi[,"PCDTC"], 1, 10) == substr(cEXSTDTC, 1, 10), colPC]
            loc = which(PCy0[,"PCDTC"] == PCy[1,"PCDTC"])
            if (length(loc) > 0) {
              if (loc[1] > 1) PCy = rbind(PCy0[loc[1] - 1,], PCy)
            }
            if ("PCTPTNUM" %in% colnames(PCi) & nrow(PCy) > 0) {
              if (PCy[nrow(PCy),"PCTPTNUM"] == PCy[1,"PCTPTNUM"] & nrow(PCy) > 1) PCy = PCy[1:(nrow(PCy)-1),]
            }
            if (sum(as.numeric(PCy[,"PCSTRESN"]) > 0) > MinPoints) {
#print(paste(i,j))
              Res0 = NCA0(EXi[j,], PCy)
              tName = names(Res0)
              Res1 = c(cStudy, cSubj, cTrt, analyte[1], cRefTime, Res0)
              names(Res1) = c("STUDYID", "USUBJID", "EXTRT", "PCTESTCD", "PCRFTDTC", tName)
              Res = rbind(Res, Res1)
            }
          } else {
            for (k in 1:nAnal) {
              cAnalyte = analyte[k]
              cEXSTDTC = EXi[j, "EXSTDTC"]
              PCy = PCi[PCi[,"PCDTC"] > EXi[j, "EXSTDTC"] & PCi["PCDTC"] < EXi[j+1, "EXSTDTC"] & PCi[,"PCTESTCD"] == cAnalyte, colPC]
              PCy0 = PCi[substr(PCi[,"PCDTC"], 1, 10) == substr(cEXSTDTC, 1, 10) & PCi["PCTESTCD"] == cAnalyte, colPC]
              loc = which(PCy0[,"PCDTC"] == PCy[1,"PCDTC"])
              if (length(loc) > 0) {
                if (loc[1] > 1) PCy = rbind(PCy0[loc[1] - 1,], PCy)
              }
              if ("PCTPTNUM" %in% colnames(PCi) & nrow(PCy) > 0) {
                if (PCy[nrow(PCy),"PCTPTNUM"] == PCy[1,"PCTPTNUM"] & nrow(PCy) > 1) PCy = PCy[1:(nrow(PCy)-1),]
              }
              if (sum(as.numeric(PCy[,"PCSTRESN"]) > 0) > MinPoints) {
#print(paste(i,j,k))
                Res0 = NCA0(EXi[j,], PCy)
                tName = names(Res0)
                Res1 = c(cStudy, cSubj, cTrt, cAnalyte, cRefTime, Res0)
                names(Res1) = c("STUDYID", "USUBJID", "EXTRT", "PCTESTCD", "PCRFTDTC", tName)
                Res = rbind(Res, Res1)
              }
            }
          }
        } else {
          if (nAnal == 1) {
            cEXSTDTC = EXi[j, "EXSTDTC"]
            PCy = PCi[PCi[,"PCDTC"] > cEXSTDTC, c("PCTESTCD", "PCDTC", "PCSTRESN", "PCSTRESU")]
            PCy0 = PCi[substr(PCi[,"PCDTC"], 1, 10) == substr(cEXSTDTC, 1, 10), c("PCTESTCD", "PCDTC", "PCSTRESN", "PCSTRESU")]
            loc = which(PCy0[,"PCDTC"] == PCy[1,"PCDTC"])
            if (length(loc) > 0) {
              if (loc[1] > 1) PCy = rbind(PCy0[loc[1] - 1,], PCy)
            }
            if (sum(as.numeric(PCy[,"PCSTRESN"]) > 0) > MinPoints) {
#print(paste(i,j))
              Res0 = NCA0(EXi[j,], PCy)
              tName = names(Res0)
              Res1 = c(cStudy, cSubj, cTrt, analyte[1], cRefTime, Res0)
              names(Res1) = c("STUDYID", "USUBJID", "EXTRT", "PCTESTCD", "PCRFTDTC", tName)
              Res = rbind(Res, Res1)
            }
          } else {
            for (k in 1:nAnal) {
              cAnalyte = analyte[k]
              cEXSTDTC = EXi[j, "EXSTDTC"]
              PCy = PCi[PCi[,"PCDTC"] > cEXSTDTC & PCi[,"PCTESTCD"] == cAnalyte, c("PCTESTCD", "PCDTC", "PCSTRESN", "PCSTRESU")]
              PCy0 = PCi[substr(PCi[,"PCDTC"], 1, 10) == substr(cEXSTDTC, 1, 10) & PCi["PCTESTCD"] == cAnalyte, c("PCTESTCD", "PCDTC", "PCSTRESN", "PCSTRESU")]
              loc = which(PCy0[,"PCDTC"] == PCy[1,"PCDTC"])
              if (length(loc) > 0) {
                if (loc[1] > 1) PCy = rbind(PCy0[loc[1] - 1,], PCy)
              }
              if (sum(as.numeric(PCy[,"PCSTRESN"]) > 0) > MinPoints) {
#print(paste(i,j,k))
                Res0 = NCA0(EXi[j,], PCy)
                tName = names(Res0)
                Res1 = c(cStudy, cSubj, cTrt, cAnalyte, cRefTime, Res0)
                names(Res1) = c("STUDYID", "USUBJID", "EXTRT", "PCTESTCD", "PCRFTDTC", tName)
                Res = rbind(Res, Res1)
              }
            }
          }
        }
      }
    }
  }
  Res = Res[order(as.character(Res[,"STUDYID"]), as.character(Res[,"USUBJID"]), as.character(Res[,"PCRFTDTC"]), as.character(Res[,"PCTESTCD"])),]
  Res = conv.pp(Res)
  attr(Res, "NCA") = "ncaRes"
  class(Res) = union(class(Res), "pp")
  return(Res)
}
