conv.pp = function(nca)
{
  numericCols = c("Dose", "Dur", "CMAX", "CMAXD", "TMAX", "TLAG", "CLST", "CLSTP", "TLST", "LAMZHL", "LAMZ",
                  "LAMZLL", "LAMZUL" , "LAMZNPT","CORRXY", "R2", "R2ADJ", "AUCLST", "AUCALL", "AUCIFO",
                  "AUCIFOD", "AUCPEO", "AUCIFP", "AUCIFPD", "AUCPEP", "AUMCLST", "AUMCIFO",
                  "AUMCPEO", "AUMCIFP", "AUMCPEP",
                  "C0", "MRTIVLST", "MRTIVIFO", "MRTIVIFP", "MRTEVLST", "MRTEVIFO", "MRTEVIFP",
                  "VZFO", "VZFP", "CLFO", "CLFP", "VZO", "VZP", "CLO", "CLP", "VSSO", "VSSP")

  colNames = colnames(nca)
  DIM = dim(nca)
  nca =  matrix(unlist(nca), nrow=DIM[1], ncol=DIM[2])
  colnames(nca) = colNames

  Res = data.frame(nca[,1])
  for (i in 2:ncol(nca)) {
    cCol = colNames[i]
    if (cCol %in% numericCols) {
      Res = cbind(Res, as.numeric(nca[,cCol]))
    } else {
      Res = cbind(Res, nca[,cCol])
    }
  }
  colnames(Res) = colnames(nca)
  attr(Res, "NCA") = "pp"
  return(Res)
}
