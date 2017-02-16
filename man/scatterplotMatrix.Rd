\name{scatterplotMatrix}
\alias{scatterplotMatrix}
\alias{scatterplotMatrix.formula}
\alias{scatterplotMatrix.default}
\alias{spm}
\title{Scatterplot Matrices}
\description{
  Enhanced scatterplot matrices with univariate displays down the diagonal;
  \code{spm} is an abbreviation for \code{scatterplotMatrix}.
  This function just sets up a call to \code{pairs} with custom panel functions.
}
\usage{

scatterplotMatrix(x, ...)

\method{scatterplotMatrix}{formula}(formula, data=NULL, subset, id=FALSE, ...)

\method{scatterplotMatrix}{default}(x, smooth=TRUE, id=FALSE, legend=TRUE, 
    var.labels=colnames(x), 
    diagonal=c("density", "boxplot", "histogram", "oned", "qqplot", "none"), 
    adjust=1, nclass, plot.points=TRUE, spread = !by.groups, reg.line=lm,
    transform=FALSE, family=c("bcPower", "yjPower"), ellipse=FALSE, 
    groups=NULL, by.groups=FALSE, use=c("complete.obs", "pairwise.complete.obs"),
    col=if (n.groups == 1) palette()[3:1] else rep(palette(), length=n.groups),
    pch=1:n.groups, lwd=1, lty=1, 
    cex=par("cex"), cex.axis=par("cex.axis"), cex.labels=NULL, 
    cex.main=par("cex.main"), row1attop=TRUE, ...)

spm(x, ...)
}

\arguments{
  \item{x}{a data matrix, numeric data frame.}
  \item{formula}{a one-sided \dQuote{model} formula, of the form 
    \code{ ~ x1 + x2 + ... + xk} or \code{ ~ x1 + x2 + ... + xk | z} where \code{z}
    evaluates to a factor or other variable to divide the data into groups.}
  \item{data}{for \code{scatterplotMatrix.formula},
    a data frame within which to evaluate the formula.}
  \item{subset}{expression defining a subset of observations.}
  \item{smooth}{specifies the smoother to be used along with its arguments; if \code{FALSE}, no smoother is shown;
    can be a list giving the smoother function and its named arguments; \code{TRUE}, the default, is equivalent to
    \code{list(smoother=loessLine)}. See \code{\link{ScatterplotSmoothers}} for the smoothers supplied by the
    \pkg{car} package and their arguments.}
  \item{id}{controls point identification; if \code{FALSE} (the default), no points are identified;
    can be a list of named arguments to the \code{\link{showLabels}} function;
    \code{TRUE} is equivalent to \code{list(method="mahal", n=2, cex=1, col=palette()[1], location="lr")},
    which identifies the 2 points (in each group, if \code{by.groups=TRUE}) with the largest Mahalanobis distances from the center
    of the data; \code{list(method="identify")} for interactive point identification is not allowed.}
  \item{legend}{controls placement of a legend if the plot is drawn by groups; if \code{FALSE}, the legend
    is suppressed. Can be a list with the named element\code{coords} specifying the position of the legend
    in any form acceptable to the \code{\link{legend}} function; 
    \code{TRUE} (the default) is equivalent to \code{list(coords=NULL)}, for which placement will vary
    by the the value of the \code{diagonal} argument---e.g., \code{"topright"} for \code{diagonal="density"}.}
  \item{var.labels}{variable labels (for the diagonal of the plot).}
  \item{diagonal}{contents of the diagonal panels of the plot. If plotting by groups, a different
    univariate display (with the exception of \code{"histogram"}) will be drawn for each group.}
  \item{adjust}{relative bandwidth for density estimate, passed to 
    \code{density} function.}
  \item{nclass}{number of bins for histogram, passed to \code{hist}
    function.}
  \item{plot.points}{if \code{TRUE} the points are plotted in each 
    off-diagonal panel.}
  \item{spread}{if TRUE, estimate the (square root) of the variance function.  For \code{loessLine} and
    for \code{gamLine}, this is done by separately smoothing the squares of the postive and negative
    residuals from the mean fit, and then adding the square root of the fitted values to the mean fit.  For
    \code{quantregLine}, fit the .25 and .75 quantiles with a quantile regression additive model.
    The default is TRUE if \code{by.groups=FALSE} and FALSE is \code{by.groups=TRUE}.}
  \item{reg.line}{if not \code{FALSE} a line is plotted using the
    function given by this argument; e.g., using \code{rlm} in
    package \code{MASS} plots a robust-regression line.}
  \item{transform}{if \code{TRUE}, multivariate normalizing power transformations
    are computed with \code{\link{powerTransform}}, rounding the estimated powers to `nice' values for plotting; 
    if a vector of powers, one for each variable, these are applied prior to plotting. If  there are \code{groups}
    and \code{by.groups} is \code{TRUE}, then the transformations are estimated \emph{conditional} on the
    \code{groups} factor.}
  \item{family}{family of transformations to estimate: \code{"bcPower"} for the Box-Cox family or 
    \code{"yjPower"} for the Yeo-Johnson family (see \code{\link{powerTransform}}).}
  \item{ellipse}{controls plotting data-concentration ellipses. If \code{FALSE} (the default), no ellipses are plotted.
    Can be a list of named values giving \code{levels}, a vector of one or more bivariate-normal probability-contour levels at which to
    plot the ellipses; and \code{robust}, a logical value determing whether to use the \code{\link{cov.trob}} function in the \pkg{MASS} package
    to calculate the center and covariance matrix for the data ellipses. \code{TRUE} is equivalent to \code{list(levels=c(.5, .95), robust=TRUE)}.}
  \item{groups}{a factor or other variable dividing the data into groups; groups are
    plotted with different colors and plotting characters.}
  \item{by.groups}{if \code{TRUE}, regression lines are fit by groups.}
  \item{use}{if \code{"complete.obs"} (the default), cases with missing data are omitted; if \code{"pairwise.complete.obs"), all valid cases are used
    in each panel of the plot.}}
  \item{pch}{plotting characters for points; default is the plotting characters in
    order (see \code{\link[graphics]{par}}).}
  \item{col}{colors for lines and points; the default is taken from the color palette, 
  	with \code{palette()[3]} for linear regression lines, \code{palette()[2]} for nonparametric regression lines, and \code{palette()[1]} 
  	for points if there are
  	no groups, and successive colors for the groups if there are groups.}
  \item{lwd}{width of linear-regression lines (default \code{1}).}
  \item{lty}{type of linear-regression lines (default \code{1}, solid line).}
  \item{cex, cex.axis, cex.labels, cex.main}{set sizes of various graphical elements
    (see \code{\link[graphics]{par}}).}
  \item{row1attop}{If \code{TRUE} (the default) the first row is at the top, as in a matrix, as
  	opposed to at the bottom, as in graph (argument suggested by Richard Heiberger).}
  \item{...}{arguments to pass down.}
}

\value{
  \code{NULL}. This function is used for its side effect: producing
  a plot.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\references{
  Fox, J. and Weisberg, S. (2011) 
  \emph{An R Companion to Applied Regression}, Second Edition, Sage.
}

\seealso{\code{\link[graphics]{pairs}}, \code{\link{scatterplot}},
  \code{\link{dataEllipse}}, \code{\link{powerTransform}}, 
  \code{\link{bcPower}}, \code{\link{yjPower}}, \code{\link[MASS]{cov.trob}},
  \code{\link{showLabels}}, \code{\link{ScatterplotSmoothers}}.}

\examples{
scatterplotMatrix(~ income + education + prestige | type, data=Duncan)
scatterplotMatrix(~ income + education + prestige, 
    transform=TRUE, data=Duncan, id=TRUE, smooth=list(smoother=gamLine))
scatterplotMatrix(~ income + education + prestige | type, smooth=FALSE, 
	by.group=TRUE, transform=TRUE, ellipse=TRUE, data=Duncan)
}

\keyword{hplot}
