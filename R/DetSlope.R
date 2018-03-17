# interactive selcting points for terminal slope
DetSlope = function(x, y)
{
# Check input
  x = x[y != 0]
  y = y[y != 0]
  y = log(y)           # for concentration-time plot only
  n1 = length(x)
  if (n1 != length(y)) stop("Length of A and B should be same.")
  if (sum(c(is.na(x), is.na(y))) > 0) stop("NAs are not allowed.")
  if ((is.numeric(x) & is.numeric(y)) != TRUE) stop("Only numeric vectors are allowed")

# Preparation of environment and variables
  OldOpt = options("locatorBell")
  options(locatorBell = FALSE) # locator bell is annoying

  dev.new()            # not to overwrite previous plot
  DefPar = par(bg="white")
  plot(x, y, xlab="Time", ylab="log(Concentration)", main="Choose points for terminal slope") # for time-concentration only

  sel = rep(FALSE, n1) # selection indictator
  sel.1 = 0            # minimum from which(sel)
  sel.2 = 0            # maximum from which(sel)
  a = 0                # intercept of regression line
  b = 0                # slope of regression line

  while (TRUE) {
# Receive a point
    ans = identify(x, y, plot=FALSE, n=1)
    if (length(ans) == 0) break

    if (sum(sel) == 1 & sel.1 != ans) { # ans is the second point
      sel[sel.1:ans] = TRUE
    } else {
      sel[ans] = !sel[ans]
    }

# Clear and Redraw
    if (a != 0 | b != 0) abline(a=a, b=b, col="white")
    points(x, y, pch=21, col="white", bg="white", cex=2)

    if (sum(sel) > 1) {          # 2 points or more selected
      sel.1 = min(which(sel == TRUE))
      sel.2 = max(which(sel == TRUE))
      z = lm(y[sel] ~ x[sel])
      a = z$coefficients[1]
      b = z$coefficients[2]
      abline(a=a, b=b)
      legend("topright", paste0("Adj. R-square=", format(summary(z)$adj.r.squared, digits=3)), box.col="white", inset=0.01)
      points(x[sel], y[sel], pch=16)
      points(x[!sel], y[!sel], pch=4)
    } else if (sum(sel) == 1) {  # 1 point selected
      sel.1 = which(sel == TRUE)
      sel.2 = 0
      points(x[sel], y[sel], pch=16)
      points(x[!sel], y[!sel])
    } else {                     # no point selected -> reset
      sel.1 = 0
      sel.2 = 0
      points(x, y)
    }
  }

# Restore and Return
  par(DefPar)
  options(locatorBell = OldOpt)
  Res = Slope(x[sel], y[sel])
  attr(Res, "Points") = which(sel)
  return(Res)
}
