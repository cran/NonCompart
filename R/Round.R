#BasicUtil.R

Round = function(x, n=0) sign(x)*trunc(abs(x)*10^n + 0.5)/10^n

