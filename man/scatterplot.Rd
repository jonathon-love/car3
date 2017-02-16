\name{scatterplot}
\alias{scatterplot}
\alias{scatterplot.formula}
\alias{scatterplot.default}
\alias{sp}

\title{Scatterplots with Boxplots}
\description{
  Makes enhanced scatterplots, with boxplots in the margins, a nonparametric regression smooth, 
  smoothed conditional spread, outlier identification, and a regression line; 
  \code{sp} is an abbreviation for \code{scatterplot}.
}
\usage{
scatterplot(x, ...)

\method{scatterplot}{formula}(formula, data, subset, xlab, ylab, id=FALSE, 
    legend=TRUE, ...)

\method{scatterplot}{default}(x, y, smooth=TRUE, id=FALSE, legend=TRUE,
    spread=!by.groups, reg.line=lm,
    boxplots=if (by.groups) "" else "xy",
    xlab=deparse(substitute(x)), ylab=deparse(substitute(y)), las=par("las"),
    lwd=1, lty=1, 
    log="", jitter=list(), xlim=NULL, ylim=NULL,
    cex=par("cex"), cex.axis=par("cex.axis"), cex.lab=par("cex.lab"), 
    cex.main=par("cex.main"), cex.sub=par("cex.sub"), 
    groups, by.groups=!missing(groups), ellipse=FALSE,
    col=if (n.groups == 1) palette()[3:1] else rep(palette(), length=n.groups),
    pch=1:n.groups, 
    reset.par=TRUE, grid=TRUE, ...)

sp(x, ...)
}

\arguments{
  \item{x}{vector of horizontal coordinates (or first argument of generic function).}
  \item{y}{vector of vertical coordinates.}
  \item{formula}{a ``model'' formula, of the form \code{y ~ x} or 
    (to plot by groups) \code{y ~ x | z}, where \code{z} evaluates to a factor 
    or other variable dividing the data into groups. If \code{x} is a factor, then parallel boxplots
    are produced using the \code{\link{Boxplot}} function.}
  \item{data}{data frame within which to evaluate the formula.}
  \item{subset}{expression defining a subset of observations.}
  \item{smooth}{specifies the smoother to be used along with its arguments; if \code{FALSE}, no smoother is shown;
    can be a list giving the smoother function and its named arguments; \code{TRUE}, the default, is equivalent to
    \code{list(smoother=loessLine)}. See \code{\link{ScatterplotSmoothers}} for the smoothers supplied by the
    \pkg{car} package and their arguments.}
  \item{id}{controls point identification; if \code{FALSE} (the default), no points are identified;
    can be a list of named arguments to the \code{\link{showLabels}} function;
    \code{TRUE} is equivalent to \code{list(method="mahal", n=2, cex=1, col=palette()[1], location="lr")},
    which identifies the 2 points (in each group) with the largest Mahalanobis distances from the center
    of the data.}
  \item{legend}{controls placement and properties of a legend if the plot is drawn by groups; if \code{FALSE}, the legend
    is suppressed. Can be a list of named arguments, as follows: \code{title} for the legend; \code{inset}, giving space
    (as a proportion of the axes) to offset the legend from the axes; \code{coords} specifying the position of the legend
    in any form acceptable to the \code{\link{legend}} function or, if not given, the legend is placed \emph{above}
    the plot in the upper margin; \code{columns} for the legend, determined automatically to prefer a horizontal layout if not given explicitly.
    \code{TRUE} (the default) is equivalent to \code{list(title=deparse(substitute(groups)), inset=0.02)}.}
  \item{spread}{if TRUE, estimate the (square root) of the variance function.  For the \code{\link{loessLine}} and
    for \code{\link{gamLine}} smoothers, this is done by separately smoothing the squares of the postive and negative
    residuals from the mean fit, and then adding the square root of the fitted values to the mean fit.  For
    \code{\link{quantregLine}}, the functions shows the .25 and .75 quantiles from a quantile regression additive model.
    The default is \code{TRUE} if \code{by.groups=FALSE} and \code{FALSE} is \code{by.groups=TRUE}.}
  \item{reg.line}{function to draw a regression line on the plot
    or \code{FALSE} not to plot a regression line.}
  \item{boxplots}{if \code{"x"} a boxplot for \code{x} is drawn below the plot;
    if \code{"y"} a boxplot for \code{y} is drawn to the left of the plot;
    if \code{"xy"} both boxplots are drawn; set to \code{""} or \code{FALSE} to
    suppress both boxplots.}
  \item{xlab}{label for horizontal axis.}
  \item{ylab}{label for vertical axis.}
  \item{las}{if \code{0}, ticks labels are drawn parallel to the
    axis; set to \code{1} for horizontal labels (see \code{\link[graphics]{par}}).}
  \item{lwd}{width of linear-regression lines (default \code{1}).}
  \item{lty}{type of linear-regression lines (default \code{1}, solid line).}
  \item{log}{same as the \code{log} argument to \code{\link{plot}}, to produce log axes.}
  \item{jitter}{a list with elements \code{x} or \code{y} or both, specifying jitter factors
    for the horizontal and vertical coordinates of the points in the scatterplot. The 
    \code{\link[base]{jitter}} function is used to randomly perturb the points; specifying a
    factor of \code{1} produces the default jitter. Fitted lines are unaffected by the jitter.}
  \item{xlim}{the x limits (min, max) of the plot; if \code{NULL}, determined from the data.}
  \item{ylim}{the y limits (min, max) of the plot; if \code{NULL}, determined from the data.}
  \item{groups}{a factor or other variable dividing the data into groups; groups are
    plotted with different colors and plotting characters.}
  \item{by.groups}{if \code{TRUE}, regression lines are fit by groups.}
  \item{ellipse}{controls plotting data-concentration ellipses. If \code{FALSE} (the default), no ellipses are plotted.
    Can be a list of named values giving \code{levels}, a vector of one or more bivariate-normal probability-contour levels at which to
    plot the ellipses; and \code{robust}, a logical value determing whether to use the \code{\link{cov.trob}} function in the \pkg{MASS} package
    to calculate the center and covariance matrix for the data ellipses. \code{TRUE} is equivalent to \code{list(levels=c(.5, .95), robust=TRUE)}.}
  \item{col}{colors for lines and points; the default is taken from the color palette, 
  	with \code{palette()[3]} for linear regression lines, \code{palette()[2]} for nonparametric regression lines, and \code{palette()[1]} 
  	for points if there are
  	no groups, and successive colors for the groups if there are groups.}
  \item{pch}{plotting characters for points; default is the plotting characters in
    order (see \code{\link[graphics]{par}}).}
  \item{cex, cex.axis, cex.lab, cex.main, cex.sub}{set sizes of various graphical elements; 
    (see \code{\link[graphics]{par}}).}
  \item{reset.par}{if \code{TRUE} (the default) then plotting parameters are reset to their previous values
    when \code{scatterplot} exits; if \code{FALSE} then the \code{mar} and \code{mfcol} parameters are
    altered for the current plotting device. Set to \code{FALSE} if you want to add graphical elements
    (such as lines) to the plot.}
  \item{grid}{If TRUE, the default, a light-gray background grid is put on the graph}
  \item{\dots}{other arguments passed down and to \code{plot}.}
}

\value{
  If points are identified, their labels are returned; otherwise \code{NULL} is returned invisibly.
}

\author{John Fox \email{jfox@mcmaster.ca}}


\seealso{\code{\link[graphics]{boxplot}},
  \code{\link[base]{jitter}}, \code{\link{legend}},
  \code{\link{scatterplotMatrix}}, \code{\link{dataEllipse}}, \code{\link{Boxplot}},
  \code{\link[MASS]{cov.trob}},
  \code{\link{showLabels}}, \code{\link{ScatterplotSmoothers}}.}

\examples{
scatterplot(prestige ~ income, data=Prestige, ellipse=TRUE)

scatterplot(prestige ~ income, data=Prestige, smooth=list(smoother=quantregLine))

scatterplot(prestige ~ income | type, data=Prestige, 
            smooth=list(smoother=loessLine, span=1))

scatterplot(prestige ~ income | type, data=Prestige, legend=list(coords="topleft"))

scatterplot(vocabulary ~ education, jitter=list(x=1, y=1), 
            data=Vocab, smooth=FALSE, lwd=3)

scatterplot(infant.mortality ~ gdp, log="xy", data=UN, id=list(n=5))

scatterplot(income ~ type, data=Prestige)

if(interactive()){ # remember to exit from point-identification mode
    scatterplot(infant.mortality ~ gdp, id=list(method="identify"), data=UN)
}
}


\keyword{hplot}
