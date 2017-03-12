unit = function(code, uTime="h", uConc="ng/mL", uDose="mg")
{
  fAmt = c(1, 1e3, 1e6, 1e9, 1e12)
  names(fAmt) = c("g", "mg", "ug", "ng", "pg")

  code = toupper(code)
  uTime = tolower(uTime)
  uDose = tolower(uDose)
  if (toupper(uConc) == toupper("mg/mL")) uConc = "g/L"
  if (toupper(uConc) == toupper("ug/mL")) uConc = "mg/L"
  if (toupper(uConc) == toupper("ng/mL")) uConc = "ug/L"
  if (toupper(uConc) == toupper("pg/mL")) uConc = "ng/L"

  if (code %in% c("CMAX", "CLST", "CLSTP", "C0")) return(c(uConc, 1))
  if (code == "CMAXD") return(c(paste0(uConc,"/",uDose), 1))

  if (code %in% c("TMAX", "TLAG", "TLST", "LAMZHL", "LAMZLL", "LAMZUL", "MRTIVLST",
                  "MRTIVIFO", "MRTIVIFP", "MRTEVLST", "MRTEVIFO", "MRTEVIFP")) {
    return(c(uTime, 1))
  }

  if (code == "LAMZ") return(c(paste0("/",uTime), 1))
  if (code %in% c("LAMZNPT", "CORRXY", "R2", "R2ADJ")) return(c("", 1))
  if (code %in% c("AUCLST", "AUCALL", "AUCIFO", "AUCIFP")) return(c(paste0(uTime,"*",uConc), 1))
  if (code %in% c("AUCIFOD", "AUCIFPD")) return(c(paste0(uTime,"*",uConc,"/",uDose), 1))
  if (code %in% c("AUCPEO", "AUCPEP", "AUCPBEO", "AUCPBEP", "AUMCPEO", "AUMCPEP")) return(c("%", 1))
  if (code %in% c("AUMCLST", "AUMCIFO", "AUMCIFP")) return(c(paste0(uTime,"2*",uConc), 1))

  tConc = strsplit(uConc, "/")[[1]]
  uAmt = tConc[1]
  uVol = tConc[2]

  if (code %in% c("VZO", "VZP", "VZFO", "VZFP", "VSSO", "VSSP")) {
    return(c(uVol, fAmt[uAmt]/fAmt[uDose]))
  }
  if (code %in% c("CLO", "CLP", "CLFO", "CLFP")) {
    return(c(paste0(uVol,"/",uTime), fAmt[uAmt]/fAmt[uDose]))
  }

  stop("Error: unknown PPTESTCD code!")
}
