combXPT = function(folders, domain)
{
  nFolder = length(folders)
  if (length(nFolder) == 0) stop("You did not specify any folder!")
  if (length(domain) != 1) stop("One domain name needs.")
  folders = Trim(folders)
  domain = toupper(Trim(domain))

  for (i in 1:nFolder) {
    cFolder = folders[i]
    if (!(substr(cFolder, nchar(cFolder), nchar(cFolder)) %in% c("/", "\\"))) folders[i] = paste0(cFolder, "/")
  }

#  require(foreign)
  for (i in 1:nFolder) {
    cFileName = paste0(folders[i], domain, ".XPT")
    if (file.exists(cFileName)) {
      cXPT = foreign::read.xport(cFileName)
      for (j in 1:ncol(cXPT)) {
        cXPT[,j] = as.character(cXPT[,j])
        cXPT[is.na(cXPT[,j]),j] = ""
      }
      if (i == 1) {
        colNames = colnames(cXPT)
        XPT = cXPT
      } else {
        colNames = intersect(colNames, colnames(cXPT))
        XPT = rbind(XPT[colNames], cXPT[colNames])
      }
    }
  }
  return(XPT)
}
