BestSlope = function(x, y, adm="Extravascular", TOL=1e-4, excludeDelta=1)
{
  Result = c(R2      = NA_real_, # R square
             R2ADJ   = NA_real_, # R square adjusted
             LAMZNPT = 0,        # Number of points for Lambda z
             LAMZ    = NA_real_, # Lambda z, terminal slope as a positive number
             b0      = NA_real_, # intercept from OLS, i.e. simple linear regression
             CORRXY  = NA_real_, # Correlation of x, y
             LAMZLL  = NA_real_, # Lower time for lambda z
             LAMZUL  = NA_real_, # Upper time for lambda z
             CLSTP   = NA_real_) # Clast predicted
# Input Check
  if (excludeDelta < 0) stop("Option excludeDelta should be non-negative!")
  n = length(x)
  if (n == 0 | n != length(y) | !is.numeric(x) | !is.numeric(y) | length(y[y < 0]) > 0) {
    Result["LAMZNPT"] = 0
    return(Result)
  }

  if (length(unique(y)) == 1) { # Case of all the same values
    Result["LAMZNPT"] = 0
    Result["b0"] = unique(y)
    return(Result)
  }

  admUpper = toupper(trimws(adm))
  r0 = .findBestSlope(x, y, admUpper, TOL, Result)
  r0["CLSTP"] = if (r0["LAMZNPT"] > 0) exp(r0["b0"] - r0["LAMZ"]*max(x[is.finite(y)])) else NA

  if (excludeDelta < 1) {
    r1 = .findBestSlope(x[-n], y[-n], admUpper, TOL, Result)
    r1["CLSTP"] = if (r1["LAMZNPT"] > 0) exp(r1["b0"] - r1["LAMZ"]*max(x[is.finite(y)])) else NA

    if (is.na(r1["R2ADJ"])) {
      Result = r0
    } else if (is.na(r0["R2ADJ"])) {
      Result = r1
    } else if (r1["R2ADJ"] - r0["R2ADJ"] > excludeDelta) {
      Result = r1
    } else {
      Result = r0
    }
  } else {
    Result = r0
  }

  if (Result["LAMZNPT"] > 0) {
    attr(Result, "UsedPoints") = which(x==Result["LAMZLL"]):which(x==Result["LAMZUL"])
  } else {
    attr(Result, "UsedPoints") = NULL
  }

  return(Result)
}

# Helper: find best slope for a given set of points
.findBestSlope = function(x, y, admUpper, TOL, defaultResult)
{
  r = defaultResult

  if (admUpper == "BOLUS") {
    locStart = which.max(y)      # From Tmax (for Bolus)
  } else {
    locStart = which.max(y) + 1  # From next to Tmax (for the others)
  }
  locLast = max(which(y > 0))    # Till non-zero concentration

  if (is.na(locStart) | is.na(locLast)) {
    r["LAMZNPT"] = 0
    return(r)
  }

  if (locLast - locStart < 2) {  # Too few to fit, if this is 2, R2ADJ becomes NaN.
    r["LAMZNPT"] = 0
    return(r)
  }

  nSlopes = locLast - locStart - 1
  tmpMat = matrix(nrow=nSlopes, ncol=length(r))
  colnames(tmpMat) = names(r)
  for (i in locStart:(locLast - 2)) {
    tmpMat[i - locStart + 1, 1:8] = Slope(x[i:locLast], log(y[i:locLast]))
  }
  tmpMat = tmpMat[is.finite(tmpMat[, "R2ADJ"]) & tmpMat[, "LAMZNPT"] > 2, , drop=FALSE]

  if (nrow(tmpMat) > 0) {
    maxAdjRsq = max(tmpMat[, "R2ADJ"])
    OKs = abs(maxAdjRsq - tmpMat[, "R2ADJ"]) < TOL
    nMax = max(tmpMat[OKs, "LAMZNPT"])
    r = tmpMat[OKs & tmpMat[, "LAMZNPT"] == nMax,]
  } else {
    r["LAMZNPT"] = 0
  }

  return(r)
}
