\name{IntAUC}
\alias{IntAUC}
\title{Calculate interval AUC}
\description{calculate interval AUC}
\usage{IntAUC(x, y, t1, t2, Res, Method = "Linear")}
\arguments{
  \item{x}{vector values of x-axis, usually time}
  \item{y}{vector values of y-axis, usually concentration}
  \item{t1}{start time for AUC}
  \item{t2}{end time for AUC}
  \item{Res}{result from \code{IndiNCA} function}
  \item{Method}{one of \code{"Linear"} or \code{"Log"} to indicate the way to calculate AUC}
}
\details{
This calculate interval (partial) AUC (from t1 to t2) with the given series of x and y.
If t1 and/or t2 cannot be found within x vector, it interpolates according to the Method.
}
\value{return interval AUC value (scalar)}
\author{Kyun-Seop Bae <k@acr.kr>}
\seealso{\code{\link{AUC}}, \code{\link{Interpol}}}
\examples{
Res = IndiNCA(Theoph[Theoph$Subject==1,"Time"], Theoph[Theoph$Subject==1, "conc"], Dose=320)
IntAUC(Theoph[Theoph$Subject==1, "Time"], Theoph[Theoph$Subject==1, "conc"], t1=0.5, t2=11, Res)
}
\keyword{AUC}
\keyword{interval AUC}
\keyword{partial AUC}
