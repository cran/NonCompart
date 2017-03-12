IntAUC = function(x, y, t1, t2, Res, fit="Linear")
{
  n = length(x)
  if (n != length(y) | !is.numeric(x) | !is.numeric(y)) stop("Bad Input!")
  if (t1 > Res["TLST"]) stop("Start time of interval AUC is after Tlast.")

  tL = Res["TLST"]
  if (t2 > tL & is.na(Res["LAMZ"])) return(NA)

  newSeries = Interpol(x, y, t1, Res["LAMZ"], Res["b0"], fit=fit)
  newSeries = Interpol(newSeries[[1]], newSeries[[2]], t2, Res["LAMZ"], Res["b0"], fit=fit)
  x = newSeries[[1]]
  y = newSeries[[2]]

  if (fit=="Linear") {
    if (t2 <= tL) {
      ResIntAUC = LinAUC(x[x>=t1 & x<=t2], y[x>=t1 & x<=t2])[[1]]
    } else {
      ResIntAUC = LinAUC(x[x>=t1 & x<=tL], y[x>=t1 & x<=tL])[[1]] + LogAUC(x[x>=tL & x<=t2], y[x>=tL & x<=t2])[[1]]
    }
  } else if (fit=="Log") {
    ResIntAUC = LogAUC(x[x>=t1 & x<=t2], y[x>=t1 & x<=t2])[[1]]
  } else stop("Unknown fit!")

  return(ResIntAUC)
}
