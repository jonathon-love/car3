\name{Summarize}
\alias{Summarize}
\alias{Vcov}
\alias{Confint}
\alias{Summarize.lm}
\alias{Summarize.default}
\alias{print.Summarize.lm}
\alias{Vcov.lm}
\alias{Vcov.default}
\alias{Confint.lm}
\alias{Confint.default}

\title{
Modified Functions for Summarizing Linear Models
}

\description{
\pkg{car} package replacements for the \code{\link{summary}}, \code{\link{vcov}} and \code{\link{confint}} functions for \code{lm} objects, with additional arguments and different defaults.
Default methods that call the standard R functions are provided for the \code{Summarize}, \code{Vcov}, and \code{Confint} generics, so the \pkg{car} functions should be safe to use in general.
}

\usage{

Summarize(object, ...)

\S3method{Summarize}{lm}(object, correlation = FALSE, symbolic.cor = FALSE, 
  vcov. = Vcov, conf.intervals = FALSE, ...)

\S3method{print}{Summarize.lm}(x, digits = max(3, getOption("digits") - 3), 
  symbolic.cor = x$symbolic.cor, 
  signif.stars = getOption("show.signif.stars"), 
  header=TRUE, resid.summary=FALSE, adj.r.squared=FALSE, brief=FALSE, ...)
 
\method{Summarize}{default}(object, ...)

Vcov(object, ...)

\method{Vcov}{lm}(object, ...)

\method{Vcov}{default}(object, ...)

Confint(object, ...)

\method{Confint}{lm}(object, parm, level = 0.95, vcov.=Vcov, ...)

\method{Confint}{default}(object, ...)
}

\arguments{
  \item{object}{
an object of class \code{"lm"}, usually, a result of a call to \code{\link{lm}}.
}
  \item{correlation, symbolic.cor}{see \code{\link{summary.lm}}}
  \item{x,digits,signif.stars}{see \code{\link{summary.lm}}}
  \item{parm,level}{see \code{\link{confint}}}
  \item{vcov.}{
either a function for estimating the covariance matrix of the regression coefficients, e.g., \code{\link{hccm}}, or an estimated covariance matrix for the coefficients. See discussion below.
}
  \item{conf.intervals}{
if \code{TRUE}, then the output of \code{\link{confint}} is appended to the summary output.  His option is added for use when \code{vcov. = Boot}; see discussion below.
}
  \item{header}{
if \code{TRUE}, print the header for the summary output
  }
  \item{resid.summary}{
if \code{TRUE}, print the five-number summary of the residuals in the summary
  }
  \item{adj.r.squared}{
if \code{TRUE}, print the adjusted r-squared in the summary
  }
  \item{brief}{
if \code{TRUE}, set \code{header}, \code{resid.summary} and \code{adj.r.squared} to \code{FALSE}
  }
  \item{\dots}{
additional arguments, for consistency with other summary methods
}
}
\details{
All these functions mimic functions in the \pkg{stats} package for summarizing aspects of a linear model.

The \code{Summarize} and \code{Confint} functions add support for the \code{vcov.} argument, which allows specifying a covariance matrix for the regression coefficients other than the usual covariance matrix returned by the function \code{\link{vcov}}.  This argument may be either the name of a function, so that the call to \code{vcov.(object)} returns a covariance matrix, or else \code{vcov.} is set equal to a covariance matrix.  Setting \code{vcov.=hccm} uses 'proposal 3' described by Long and Ervin (2000) for a sandwich coefficient-variance estimator that may be robust against nonconstant variance.  Setting \code{vcov. = hccm(object, type = "hc2")} would use the matrix returned by the \code{hccm} function using proposal 2.  If \code{vcov. = Boot}, with a required capital '\code{B}', then the covariance matrix is computed using a case-resampling bootstrap by the call \code{cov(Boot(object)$t)}; see \code{\link{Boot}} and \code{\link{hccm}}. 

When using \code{vcov. = Boot}, we recommend setting \code{conf.intervals = TRUE}, so confidence intervals based on the Bca method will be based on the same bootstrap as used in the summary. Alternatively, users may prefer to save the bootstrap, and this can be accomplished by first setting \code{b1 <- Boot(object)}, and then calling \code{Summarize(object, vcov. = cov(b1$t))}. 

The overall F-test in the \code{Summarize} output uses the supplied covariance matrix in a call to the \code{\link{linearHypothesis}} function.

The print method has additional arguments to customize the standard \code{summary} output.  Standard output is obtained by setting \code{residual.summary=TRUE, adjusted.r.squared=TRUE} with \code{vcov. = vcov}.
}
\value{
The \code{Summarize.lm} function returns a list with all the elements shown at \code{\link{summary.lm}} plus two additional items:
\item{vcov.}{The name of the method used to estimate the covariance matrix}
\item{conf.intervals}{if these were requested, the confidence intervals for coefficient estimates}

The \code{Vcov.lm} function returns the usual covariance matrix returned by \code{\link{vcov}}.

The \code{Confint.lm} function returns either (1) the output from \code{\link{confint.lm}} if \code{vcov. = vcov}; (2) the output from \code{confint(Boot(object))} if \code{vcov. = Boot}; or (3) Wald-type confidence intervals using the supplied covariance matrix for any other choice of \code{vcov.}.
}


\references{
  Long, J. S. and Ervin, L. H. (2000)
  Using heteroscedasity consistent standard errors in the linear regression model. 
  \emph{The American Statistician} \bold{54}, 217--224.
  
  White, H. (1980)
  A heteroskedastic consistent covariance matrix estimator and a direct test of heteroskedasticity.
  \emph{Econometrica} \bold{48}, 817--838.
}
\author{
Sanford Weisberg \email{sandy@umn.edu}
}


\seealso{\code{\link{summary.lm}}, \code{\link{confint}}, \code{\link{vcov.lm}}, \code{\link{hccm}}, 
\code{\link{Boot}}, \code{\link{linearHypothesis}}}

\examples{
mod.prestige <- lm(prestige ~ education + income + type, Prestige)
Summarize(mod.prestige, vcov.=hccm)
print(Summarize(mod.prestige, conf.intervals=TRUE), brief=TRUE)
\dontrun{
Summarize(mod.prestige, vcov.=Boot)
}
Confint(mod.prestige, vcov.=hccm)
\dontrun{
Confint(mod.prestige, vcov.=Boot)
}
}


\keyword{models}
\keyword{regression}
