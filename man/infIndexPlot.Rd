\name{infIndexPlot}
\alias{infIndexPlot}
\alias{influenceIndexPlot}
\alias{infIndexPlot.lm}
\title{Influence Index Plot}
\description{
Provides index plots of Cook's distances,
leverages, Studentized residuals, and outlier significance levels for a 
regression object.
}
\usage{
infIndexPlot(model, ...)

influenceIndexPlot(model, ...)

\method{infIndexPlot}{lm}(model, vars=c("Cook", "Studentized", "Bonf", "hat"), 
    id=TRUE, grid=TRUE, main="Diagnostic Plots", ...)

}
\arguments{
  \item{model}{A regression object of class \code{lm} or \code{glm}.}
  \item{vars}{
All the quantities listed in this argument are plotted.  Use \code{"Cook"}
for Cook's distances, \code{"Studentized"} for Studentized
residuals, \code{"Bonf"} for Bonferroni p-values for an outlier test, and 
and \code{"hat"} for hat-values (or leverages).   Capitalization is optional.
All may be abbreviated by the first one or more letters.
}
  \item{main}{main title for graph}
  \item{id}{a list of named values controlling point labelling. The default, \code{TRUE}, is
    equivalent to \code{id=list(method="y", n=2, cex=1, col=palette()[1], location="lr")};
    \code{FALSE} suppresses point labelling. See \code{\link{showLabels}} for details.}
  \item{grid}{If TRUE, the default, a light-gray background grid is put on the graph.}
  \item{\dots}{Arguments passed to \code{plot}}
}

\value{
Used for its side effect of producing a graph.  Produces four index plots
of Cook's distance, Studentized Residuals, the corresponding Bonferroni 
p-values for outlier tests, and leverages.
}
\references{ 
  Cook, R. D. and Weisberg, S. (1999)
  \emph{Applied Regression, Including Computing and Graphics.} Wiley.

  Fox, J. (2008)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Second Edition. Sage.  
  
  Fox, J. and Weisberg, S. (2011) 
  \emph{An R Companion to Applied Regression}, Second Edition, Sage.
  
  Weisberg, S. (2014) 
  \emph{Applied Linear Regression}, Fourth Edition, Wiley.
}

\author{Sanford Weisberg, \email{sandy@umn.edu}}

\seealso{ \code{\link{cooks.distance}}, \code{\link{rstudent}},
  \code{\link{outlierTest}}, \code{\link{hatvalues}} }
\examples{
influenceIndexPlot(lm(prestige ~ income + education + type, Duncan))
}

\keyword{ hplot }% at least one, from doc/KEYWORDS
\keyword{ regression }% __ONLY ONE__ keyword per line
