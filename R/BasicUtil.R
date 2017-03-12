#BasicUtil.R

Trim  = function(x) gsub("^\\s+|\\s+$", "", x)

toproper = function(string)
{
  return(paste0(toupper(substr(string, 1, 1)), tolower(substr(string, 2, nchar(string)))))
}

Note = function() file.show(system.file("NOTE.txt", package="NonCompart"))

PrepPDF = function(FileName, Paper="letter", FontFamily="Courier")
{
  pdf(FileName, paper=Paper, width=8.5, height=11, family=FontFamily, title="NCA Report")
}

ClosePDF = function()
{
  dev.off()
}

AddPage = function(Cex=0.8, Header1="", Header2="", Header3="", Footer1="", Footer2="", Footer3="", PrintRowNum=FALSE, StartRowNum=1)
{
  if (Cex < 0.1 | Cex > 1.2) stop("Too small or too large cex!")
  nCol = trunc(242 - 190*Cex)
  nRow = trunc(155 - 125*Cex)
  options(width=max(10, trunc(300 - 250*Cex)))

  par(oma=c(0,0,0,0), mfrow=c(1,1), mar=c(0, 0, 0, 0), adj=0, cex=Cex)
  plot(0, 0, type="n", ylim=c(2*nRow-1, 0), xlim=c(0, nCol-1), xaxt="n", yaxt="n", ylab="", xlab="", bty="n")
  if (PrintRowNum == TRUE) for(j in 1:nRow) text(-2 - 1 + 3 - floor(log10(j + StartRowNum - 1)), 2*j - 1, paste0(j + StartRowNum - 1, ":"), offset=0)
  text(0, -3, Header1, pos=4)
  text(nCol/2, -3, Header2, pos=3)
  text(nCol, -3, Header3, pos=2)

  text(0, 2*nRow + 1, Footer1, pos=4)
  text(nCol/2, 2*nRow + 1, Footer2, pos=1)
  text(nCol, 2*nRow + 1, Footer3, pos=2)
}

Text1 = function(Row, Col, Text, Cex=0.8)
{
  if (Cex < 0.1 | Cex > 1.2) stop("Too small or too large cex!")
  nCol = trunc(242 - 190*Cex)
  nRow = trunc(155 - 125*Cex)
  if (Col > nCol | Row > nRow) stop("Text seems out of paper!")

  text(Col - 1, 2*Row - 1, Text, cex=Cex, offset=0)
}

TextM = function(MTxt, Cex=0.8,  Header1="", Header2="", Header3="", Footer2="", Footer1="", Footer3="", PrintRowNum=FALSE)
{
  if (Cex < 0.1 | Cex > 1.2) stop("Too small or too large cex!")
  nRow = trunc(155 - 125*Cex)

  for (i in 1:length(MTxt)) {
    if (i%%nRow == 1) AddPage(Cex=Cex, Header1=Header1, Header2=Header2, Header3=Header3, Footer2=Footer2, Footer1=Footer1, Footer3=Footer3, PrintRowNum=PrintRowNum)
    Text1((i - 1)%%nRow + 1, 1, MTxt[i])
  }
}

Unit = function(code="", timeUnit="h", concUnit="ng/mL", doseUnit="mg", MW=0)
{
  if (length(strsplit(doseUnit, "/")[[1]]) != 1) stop("Dose unit should not be based on body weight or BSA!")
  if (!is.numeric(MW)) stop("Molecular weight should be positive number!")
  if (MW < 0) stop("Molecular weight should be positive number!")
  
  rGram = c(1, 1e3, 1e6, 1e9, 1e12)
  names(rGram) = c("g", "mg", "ug", "ng", "pg")

  rMol = c(1, 1e3, 1e6, 1e9, 1e12)
  names(rMol) = c("mol", "mmol", "umol", "nmol", "pmol")

  doseUnit = tolower(doseUnit)
  timeUnit = tolower(timeUnit)

  if (toupper(concUnit) == toupper("mg/mL")) concUnit = "g/L"
  if (toupper(concUnit) == toupper("ug/mL")) concUnit = "mg/L"
  if (toupper(concUnit) == toupper("ng/mL")) concUnit = "ug/L"
  if (toupper(concUnit) == toupper("pg/mL")) concUnit = "ng/L"

  if (toupper(concUnit) == toupper("mmol/mL")) concUnit = "mol/L"
  if (toupper(concUnit) == toupper("umul/mL")) concUnit = "mmol/L"
  if (toupper(concUnit) == toupper("nmol/mL")) concUnit = "umol/L"
  if (toupper(concUnit) == toupper("pmol/mL")) concUnit = "nmol/L"

  tConc = strsplit(concUnit, "/")[[1]]
  uAmt = tConc[1]
  uVol = tConc[2]
  
  if ((uAmt %in% names(rMol) & doseUnit %in% names(rGram)) | (uAmt %in% names(rGram) & doseUnit %in% names(rMol))) {
    if (is.na(MW)) warning("Molecular weight should be given for more informative results!")
    if (MW <= 0) warning("Molecular weight should be given for more informative results!")
  } 

  TestCD = c("b0", "CMAX", "CMAXD", "TMAX", "TLAG", "CLST", "CLSTP", "TLST", "LAMZHL", "LAMZ",
             "LAMZLL", "LAMZUL", "LAMZNPT", "CORRXY", "R2", "R2ADJ", "C0", "AUCLST", "AUCALL",
             "AUCIFO", "AUCIFOD", "AUCIFP", "AUCIFPD", "AUCPEO", "AUCPEP", "AUCPBEO", "AUCPBEP",
             "AUMCLST", "AUMCIFO", "AUMCIFP", "AUMCPEO", "AUMCPEP",
             "MRTIVLST", "MRTIVIFO", "MRTIVIFP", "MRTEVLST", "MRTEVIFO", "MRTEVIFP",
             "VZO", "VZP", "VZFO", "VZFP", "CLO", "CLP", "CLFO", "CLFP", "VSSO", "VSSP")
  nTestCD = length(TestCD)
  Res = data.frame(Unit=rep("",nTestCD), Factor=rep(1,nTestCD), stringsAsFactors = FALSE)
  rownames(Res) = TestCD

  for (i in 1:nTestCD) {
    Code = TestCD[i]

    if (Code %in% c("CMAX", "CLST", "CLSTP", "C0")) Res[i, 1] = concUnit
    if (Code == "CMAXD") Res[i, 1] = paste0(concUnit,"/",doseUnit)

    if (Code %in% c("TMAX", "TLAG", "TLST", "LAMZHL", "LAMZLL", "LAMZUL", "MRTIVLST",
                    "MRTIVIFO", "MRTIVIFP", "MRTEVLST", "MRTEVIFO", "MRTEVIFP")) Res[i, 1] = timeUnit

    if (Code == "LAMZ") Res[i, 1] = paste0("/",timeUnit)
    if (Code %in% c("b0", "LAMZNPT", "CORRXY", "R2", "R2ADJ")) Res[i, 1] = ""
    if (Code %in% c("AUCLST", "AUCALL", "AUCIFO", "AUCIFP")) Res[i, 1] = paste0(timeUnit,"*",concUnit)
    if (Code %in% c("AUCIFOD", "AUCIFPD")) Res[i, 1] = paste0(timeUnit,"*",concUnit,"/",doseUnit)
    if (Code %in% c("AUCPEO", "AUCPEP", "AUCPBEO", "AUCPBEP", "AUMCPEO", "AUMCPEP")) Res[i, 1] = "%"
    if (Code %in% c("AUMCLST", "AUMCIFO", "AUMCIFP")) Res[i, 1] = paste0(timeUnit,"2*",concUnit)

    if (Code %in% c("VZO", "VZP", "VZFO", "VZFP", "VSSO", "VSSP")) {
      if (uAmt %in% names(rMol) & doseUnit %in% names(rGram)) Res[i, ] = c(uVol, rMol[uAmt]/rGram[doseUnit] / MW)
      else if (uAmt %in% names(rGram) & doseUnit %in% names(rMol)) Res[i, ] = c(uVol, rGram[uAmt]/rMol[doseUnit] * MW)
      else if (uAmt %in% names(rGram) & doseUnit %in% names(rGram)) Res[i, ] = Res[i, ] = c(uVol, rGram[uAmt]/rGram[doseUnit])
      else Res[i, ] = c(uVol, rMol[uAmt]/rMol[doseUnit])
    }
    if (Code %in% c("CLO", "CLP", "CLFO", "CLFP")) {
      if (uAmt %in% names(rMol) & doseUnit %in% names(rGram)) Res[i, ] = c(paste0(uVol,"/",timeUnit), rMol[uAmt]/rGram[doseUnit] / MW)
      else if (uAmt %in% names(rGram) & doseUnit %in% names(rMol)) Res[i, ] = c(paste0(uVol,"/",timeUnit), rGram[uAmt]/rMol[doseUnit] * MW)
      else if (uAmt %in% names(rGram) & doseUnit %in% names(rGram)) Res[i, ] = Res[i, ] = c(paste0(uVol,"/",timeUnit), rGram[uAmt]/rGram[doseUnit])
      else Res[i, ] = c(paste0(uVol,"/",timeUnit), rMol[uAmt]/rMol[doseUnit])
    }
  }

  Res[,2] = as.numeric(Res[,2])
  Res[Res[,2] == 0 | Res[,2] == Inf, 2] = NA
  
  if (code == "") return(Res)
  else return(Res[code,])

}
