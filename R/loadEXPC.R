loadEXPC = function(folders)
{
# How to get nominal time without PCTPT or PCTPTNUM?
# if possible, utilize VISIT, VISITNUM, PCTPT, PCTPTNUM -> use intersect of columns

  EX = readEX(folders)
  PC = readPC(folders)

  Res = list(EX, PC)
  names(Res) = c("EX", "PC")
  attr(Res, "NCA") = "EXPC"
  return(Res)
}
