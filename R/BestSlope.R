BestSlope = function(x, y, adm="Extravascular")
{
# Last Modified: 2017.7.19
# Author: Kyun-Seop Bae k@acr.kr
# Calls : Slope
# Called by : sNCA

  RetNames = c("R2", "R2ADJ", "LAMZNPT", "LAMZ", "b0", "CORRXY", "LAMZLL", "LAMZUL", "CLSTP")
  Result = rep(NA_real_, length(RetNames))
  names(Result) = RetNames

  n = length(x)
  if (n != length(y) | !is.numeric(x) | !is.numeric(y) | length(y[y < 0]) > 0) {
    #    warning("Check Input!")
    Result["LAMZNPT"] = 0
    return(Result)
  }

  if (length(unique(y)) == 1) { # Case of all the same values
    Result["LAMZNPT"] = 0
    Result["b0"] = unique(y)
    return(Result)
  }

  if (toupper(adm) == "BOLUS") {
    locStart = which.max(y)      # From Tmax (for Bolus)
  } else {
    locStart = which.max(y) + 1  # From next to Tmax (for the others)
  }
  locLast = max(which(y > 0))    # Till non-zero concentration

  if (locLast - locStart < 2) {  # Too few to fit, R2ADJ becomes NaN.
    Result["LAMZNPT"] = 0
    return(Result)
  }

  tmpMat = matrix(nrow=(locLast - locStart - 1), ncol=length(RetNames)) # Slope function returns 9 columns
  colnames(tmpMat) = RetNames
  for (i in locStart:(locLast - 2)) {
    tmpMat[i - locStart + 1,] = Slope(x[i:locLast], log(y[i:locLast]))
  }
  tmpMat = tmpMat[tmpMat[,"LAMZNPT"] > 2,]
  if (nrow(tmpMat) > 0) {
    maxAdjRsq = max(tmpMat[,"R2ADJ"]) # The second column is "Rsq_adjusted" which is the criterion
    OKs = ifelse(abs(maxAdjRsq - tmpMat[,"R2ADJ"]) < 1e-4, TRUE, FALSE) # Tolerance is 1e-4, Phoneix WinNonlin 6.4 User's Guide p33
    nMax = max(tmpMat[OKs,"LAMZNPT"])   # Third column is "No_points_lambda_z" or "n"
    Result = tmpMat[OKs & tmpMat[,"LAMZNPT"]==nMax,]
  } else {
    Result["LAMZNPT"] = 0
  }

  return(Result)
}
