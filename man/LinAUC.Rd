\name{LinAUC}
\alias{LinAUC}
\title{Area Under the Curve(AUC) and Area Under the first Moment Curve(AUMC) by linear trapezoidal method}
\description{calculate AUC and AUMC using linear trapezoidal method}
\usage{LinAUC(x, y)}
\arguments{
  \item{x}{vector values of x-axis, usually time}
  \item{y}{vector values of y-axis, usually concentration}
}
\details{
This function returns AUC and AUMC by linear trapezoidal method.
}
\value{
  \item{AUC}{area under the curve}
  \item{AUMC}{area under the first moment curve}
}
\references{Gibaldi M, Perrier D. Pharmacokinetics 2e revised and expanded. pp 409-416. 1982}
\author{Kyun-Seop Bae <k@acr.kr>}
\seealso{\code{\link{LogAUC}},\code{\link{AUC}}}
\examples{
LinAUC(Theoph[Theoph$Subject==1, "Time"],Theoph[Theoph$Subject==1, "conc"])
AUC(Theoph[Theoph$Subject==1, "Time"],Theoph[Theoph$Subject==1, "conc"]) # compare the last line
}
\keyword{AUC}
