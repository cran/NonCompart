plotFit = function(concData, id, Time, conc, mol="", adm="Extravascular", ID="", Mol="")
{
  if (mol == "") {
    x = concData[concData[,id] == ID, Time]
    y = concData[concData[,id] == ID, conc]
  } else {
    x = concData[concData[,id] == ID & concData[,mol] == Mol, Time]
    y = concData[concData[,id] == ID & concData[,mol] == Mol, conc]
  }

  finalMat = BestSlope(x, y, adm=adm)

  y1 = ifelse(y==0, 0.1, y)
  plot(x, log(y), yaxt = "n", ylim=c(min(log(y1)), finalMat["b0"]), xlab="Time", ylab="Concentration", main=paste("Best Fit ID:", ID))
  yticks = seq(round(min(log(y1))), ceiling(max(log(y1))))
  ylabels = sapply(yticks, function(i) as.expression(bquote(e^ .(i))))
  axis(2, at=exp(yticks), labels=ylabels)
  abline(a=finalMat["b0"],b=-finalMat["LAMZ"], untf=TRUE, col="blue")

  return(finalMat)
}

