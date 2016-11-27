\name{BestSlope}
\alias{BestSlope}
\title{Choose best fit slope for the log(y) and x regression by the criteria of adjusted R-square}
\description{
It sequentially fits (log(y) ~ x) from the last point of x to the previous points with at least 3 points.
It chooses the highest adjusted R-square.
If the difference is less then 1e-4, it chooses longer slope.
}
\usage{BestSlope(x, y, AdmMode = "Extravascular")}
\arguments{
  \item{x}{vector values of x-axis, usually time}
  \item{y}{vector values of y-axis, usually concentration}
  \item{AdmMode}{one of \code{"Bolus"} or \code{"Infusion"} or \code{"Extravascular"} to indicate drug administration mode}
}
\details{
Choosing the best terminal slope (y in log scale) in pharmacokinetic analysis is somewhat challenging, and it could varies by analyist.
This function uses the same method which the most popular software uses.
Currently this function uses ordinary least square method(OLS) only.
}
\value{
  \item{Rsq}{R-squared}
  \item{adjRsq}{adjusted R-squared}
  \item{n}{number of points used for slope}
  \item{Lambda_z}{negative of slope, lambda_z}
  \item{b0}{intercept of regression line}
  \item{Corr_XY}{correlation of log(y) and x}
  \item{Lambda_z_lower}{earlist x for lambda_z}
  \item{Lambda_z_upper}{last x for lambda_z}
  \item{Clast_pred}{predicted y value at last point, concentration of last_predicted}
}
\author{Kyun-Seop Bae <k@acr.kr>}
\seealso{\code{\link{Slope}}}
\examples{
BestSlope(Theoph[Theoph$Subject==1, "Time"],Theoph[Theoph$Subject==1, "conc"])
BestSlope(Indometh[Indometh$Subject==1, "time"],Indometh[Indometh$Subject==1, "conc"],
          AdmMode="Bolus")
}
\keyword{Slope}
\keyword{best fit slope}
