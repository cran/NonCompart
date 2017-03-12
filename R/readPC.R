readPC = function(folders)
{
  PC = combXPT(folders, "PC")

# If DTC does not sec, attch :29 at the end
  for (i in 1:nrow(PC)) {
    if (nchar(PC[i, "PCDTC"]) == 16) PC[i, "PCDTC"] = paste0(PC[i, "PCDTC"], ":29")
    else PC[i, "PCDTC"] = ""
  }

# Numeric type column will be set
  colNum = intersect(c("PCSTRESN", "VISITNUM", "PCTPTNUM"), colnames(PC))
  nCol = length(colNum)
  for (i in 1:nCol) {
    PC[,colNum[i]] = as.double(PC[,colNum[i]])
  }

  if (!("PCSPEC" %in% colnames(PC))) PC[,"PCSPEC"] == "PLASMA"
  if (!("PCLLOQ" %in% colnames(PC))) PC[,"PCLLOQ"] == "0"

  return(PC)
}
