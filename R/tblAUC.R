tblAUC = function(Data, key = "Subject", colX = "Time", colY = "Y", iAUC = "", Ymax = "Emax", XofYmax = "TEmax", AUCname = "AUEClast", Outer = "NEAREST")
{
  class(Data) = "data.frame"
  nKey = length(key)
  for (i in 1:nKey) {
    if (sum(is.na(Data[, key[i]])) > 0) stop(paste(key[i], "has NA value, which is not allowed!"))
  }

  IDs = unique(as.data.frame(Data[, key], ncol=nKey))
  nID = nrow(IDs)

  Res = data.frame()
  for (i in 1:nID) {
    mask = Data[, key[1]] == IDs[i, 1]
    if (nKey > 1) {
      for (j in 2:nKey) {
        mask = mask & Data[, key[j]] == IDs[i, j]
      }
    }
    tData = Data[mask, , drop=FALSE]
    if (nrow(tData) > 0) {
      tRes = gAUC(tData[, colX], tData[, colY], Ymax=Ymax, XofYmax=XofYmax, AUCname=AUCname, iAUC=iAUC)
      Res = rbind(Res, as.data.frame(tRes))
    }
  }
  Res = cbind(IDs, Res)
  colnames(Res)[1:nKey] = key
  rownames(Res) = NULL
  return(Res)
}
