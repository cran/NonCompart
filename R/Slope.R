Slope = function(x, y)
{
# Author: Kyun-Seop Bae k@acr.kr
# Called by: BestSlope
# Calls: none except base
# Last modification: 2017.7.20
# INPUT
#    x: time
#    y: natural log of concentration
# RETURN
  Result = c(R2=NA,     # R square
             R2ADJ=NA,  # R square adjusted
             LAMZNPT=0, # Number of points for Lambda z
             LAMZ=NA,   # Lambda z, terminal slope as a positive number
             b0=NA,     # intercept from OLS, i.e. simple linear regression
             CORRXY=NA, # Correlation of x, y
             LAMZLL=NA, # Lower time for lambda z
             LAMZUL=NA, # Upper time for lambda z
             CLSTP=NA)  # Concentration last predicted in original scale
#
  n = length(x)
  if (n == 1 | n != length(y) | !is.numeric(x) | !is.numeric(y)) {
    return(Result)  # return default
  }

  mx  = mean(x)
  my  = mean(y)
  Sxx = sum((x - mx)*(x - mx))
  Sxy = sum((x - mx)*(y - my))
  Syy = sum((y - my)*(y - my))
  b1  = Sxy/Sxx

  if (is.nan(b1) | b1 > 0) {
    return(Result)   # return default
  }

# Expectedly
  Result["LAMZNPT"] = n
  Result["LAMZ"]    = -b1
  Result["b0"]      = my - b1*mx
  Result["R2"]      = b1 * Sxy/Syy
  Result["R2ADJ"]   = 1 - (1 - Result["R2"])*(n - 1)/(n - 2)
  Result["CORRXY"]  = sign(b1)*sqrt(Result["R2"])
  Result["LAMZLL"]  = x[1]
  Result["LAMZUL"]  = x[n]
  Result["CLSTP"]   = exp(Result["b0"] + b1 * x[n])
  return(Result)
}
