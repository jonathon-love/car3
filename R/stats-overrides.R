# Alternatives to `stats` functions for various reasons 12/26/2017
# Vcov.lm:  The original calls summary.lm.  The overrride is free-standing
# Summarize.lm:  Adds new argument vcov.=vcov to specify a covariance matrix.  The default reproduces the
#              current output.  The linearHypothesis function is used to compute the overall F-test.
# print.Summarize.lm: new arguments:
#              header=TRUE prints or suppresses the header
#              resid.summary=TRUE prints or suppresses the residual summary 
#              adj.r.squared=TRUE prints or suppresses printing of the adjusted r.squared
#              brief=FALSE if TRUE sets header=resid.summary=adj.r.squared=FALSE
#     In addition output is modified to include the vcov. argument if it is not set to vcov
# Confint.lm:  new argument vcov.=vcov where vcov. is either a matrix of the right size or
#     a fuction so that vcov.(object) returns an estmated covariance matrix.
#     If cov.=Boot, the returned intervals are Confint(Boot(object)), which uses the `car` function

# 2016-12-27 For now, the override functions start with a Capital letter
              
#  The original vcov.lm simply calls summary(lm-mod) to get sigmahat and (X'X)^(-1)
#  Using the Capitalized version of these functions this override is not really needed
#  2017-02-10: Renamed using uc letters; introduced default methods. J. Fox

Vcov <- function(object, ...){
  UseMethod("Vcov")
}

Vcov.default <- function (object, ...) vcov(object, ...)

Vcov.lm <- function(object, ...){
  z <- object
  p <- z$rank
  if (p == 0) return(numeric(0))
  rdf <- z$df.residual
  if (is.null(z$terms)) 
    stop("invalid 'lm' object:  no 'terms' component")
  if (!inherits(object, "lm")) 
    warning("calling summary.lm(<fake-lm-object>) ...")
  Qr <- qr(object)
  n <- NROW(Qr$qr)
  if (is.na(z$df.residual) || n - p != z$df.residual) 
    warning("residual degrees of freedom in object suggest this is not an \"lm\" fit")
  r <- z$residuals
  w <- z$weights
  rss <- if (is.null(w)) sum(r^2) else sum(w * r^2)
  resvar <- rss/rdf
  p1 <- 1L:p
  R <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
  dimnames(R) <- list(names(z$coefficients)[Qr$pivot[p1]],
                      names(z$coefficients)[Qr$pivot[p1]])
  resvar * R
}

# Summarize adds vcov. and conf.intervals arguments

Summarize <- function(object, ...){
  UseMethod("Summarize")
}

Summarize.default <- function(object, ...) summary(object, ...)

Summarize.lm <- function (object, correlation = FALSE, symbolic.cor = FALSE, 
                        vcov. = Vcov, conf.intervals=FALSE,...) {
    z <- object
    p <- z$rank
    rdf <- z$df.residual
    if (p == 0) {
      r <- z$residuals
      n <- length(r)
      w <- z$weights
      if (is.null(w)) {
        rss <- sum(r^2)
      }
      else {
        rss <- sum(w * r^2)
        r <- sqrt(w) * r
      }
      resvar <- rss/rdf
      ans <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
      class(ans) <- "Summarize.lm"
      ans$aliased <- is.na(coef(object))
      ans$residuals <- r
      ans$df <- c(0L, n, length(ans$aliased))
      ans$coefficients <- matrix(NA, 0L, 4L)
      dimnames(ans$coefficients) <- list(NULL, c("Estimate", 
                                                 "Std. Error", "t value", "Pr(>|t|)"))
      ans$sigma <- sqrt(resvar)
      ans$r.squared <- ans$adj.r.squared <- 0
      return(ans)
    }
    if (is.null(z$terms)) 
      stop("invalid 'lm' object:  no 'terms' component")
    if (!inherits(object, "lm")) 
      warning("calling summary.lm(<fake-lm-object>) ...")
#    Qr <- stats:::qr.lm(object)
    Qr <- qr(object)
    n <- NROW(Qr$qr)
    if (is.na(z$df.residual) || n - p != z$df.residual) 
      warning("residual degrees of freedom in object suggest this is not an \"lm\" fit")
    r <- z$residuals
    f <- z$fitted.values
    w <- z$weights
    if (is.null(w)) {
      mss <- if (attr(z$terms, "intercept")) 
        sum((f - mean(f))^2)
      else sum(f^2)
      rss <- sum(r^2)
    }
    else {
      mss <- if (attr(z$terms, "intercept")) {
        m <- sum(w * f/sum(w))
        sum(w * (f - m)^2)
      }
      else sum(w * f^2)
      rss <- sum(w * r^2)
      r <- sqrt(w) * r
    }
    resvar <- rss/rdf
    if (is.finite(resvar) && resvar < (mean(f)^2 + var(f)) * 
        1e-30) 
      warning("essentially perfect fit: summary may be unreliable")
    p1 <- 1L:p
    R <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
#    se <- sqrt(diag(R) * resvar)
    V <- if(is.matrix(vcov.)) vcov. else
         if(deparse(substitute(vcov.) == "Boot")) cov((b1 <- Boot(object))$t) else
         vcov.(object)
    se <- sqrt(diag(V))
    est <- z$coefficients[Qr$pivot[p1]]
    tval <- est/se
    ans <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
    ans$residuals <- r
    ans$coefficients <- cbind(est, se, tval, 2 * pt(abs(tval), 
                                                    rdf, lower.tail = FALSE))
    dimnames(ans$coefficients) <- list(names(z$coefficients)[Qr$pivot[p1]], 
                                       c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    ans$aliased <- is.na(coef(object))
    ans$sigma <- sqrt(resvar)
    ans$df <- c(p, rdf, NCOL(Qr$qr))
    if (p != attr(z$terms, "intercept")) {
      df.int <- if (attr(z$terms, "intercept")) 
        1L
      else 0L
      ans$r.squared <- mss/(mss + rss)
      ans$adj.r.squared <- 1 - (1 - ans$r.squared) * ((n - df.int)/rdf)
#      ans$fstatistic <- c(value = (mss/(p - df.int))/resvar, 
#                          numdf = p - df.int, dendf = rdf)
# linearHypothesis computes overall F test allowing for alternative covariance matrices
    mat <- diag(p - df.int)
    if(df.int==1) mat <- cbind(0, mat)
    lh <- linearHypothesis(z, mat, vcov.=V)
    ans$fstatistic <- c(value = lh$F[2], numdf = lh$Df[2], dendf = lh$Res.Df[2])
    }
    else ans$r.squared <- ans$adj.r.squared <- 0
    ans$cov.unscaled <- R
    dimnames(ans$cov.unscaled) <- dimnames(ans$coefficients)[c(1, 1)]
    if (correlation) {
      ans$correlation <- (R * resvar)/outer(se, se)
      dimnames(ans$correlation) <- dimnames(ans$cov.unscaled)
      ans$symbolic.cor <- symbolic.cor
    }
    if (!is.null(z$na.action)) 
      ans$na.action <- z$na.action
    ans$vcov. <- deparse(substitute(vcov.))
    if(conf.intervals){ 
      ans$conf.intervals <- if(ans$vcov. == "Boot") confint(b1) else Confint(object, vcov.=V)
    }
    class(ans) <- "Summarize.lm"
    ans
  } 

print.Summarize.lm <- function(x, digits = max(3, getOption("digits") - 3), symbolic.cor = x$symbolic.cor, 
          signif.stars = getOption("show.signif.stars"), 
          header=TRUE, resid.summary=FALSE, adj.r.squared=FALSE, brief=FALSE, ...) {
  if (brief) header <- resid.summary <- adj.r.squared <- FALSE
  if (header) {
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n", sep = "")
    if(x$vcov. != "vcov"){
      cat("Covariances estimated using:\n", x$vcov., "\n", sep="")}
  }
  resid <- x$residuals
  df <- x$df
  rdf <- df[2L]
  if (resid.summary) {
    cat('\n', if (!is.null(x$weights) && diff(range(x$weights))) 
      "Weighted ", "Residuals:\n", sep = "")
    if (rdf > 5L) {
      nam <- c("Min", "1Q", "Median", "3Q", "Max")
      rq <- if (length(dim(resid)) == 2L) 
        structure(apply(t(resid), 1L, quantile), dimnames = list(nam, 
                                                                 dimnames(resid)[[2L]]))
      else {
        zz <- zapsmall(quantile(resid), digits + 1)
        structure(zz, names = nam)
      }
      print(rq, digits = digits, ...)
    }
    else if (rdf > 0L) {
      print(resid, digits = digits, ...)
    }
    else {
      cat("ALL", df[1L], "residuals are 0: no residual degrees of freedom!\n")
    }
  }
  if (length(x$aliased) == 0L) {
    cat("\nNo Coefficients\n")
  }
  else {
    if (nsingular <- df[3L] - df[1L]) 
      cat("\nCoefficients: (", nsingular, " not defined because of singularities)\n", 
          sep = "")
    else cat("\nCoefficients:\n")
    coefs <- x$coefficients
    if (!is.null(aliased <- x$aliased) && any(aliased)) {
      cn <- names(aliased)
      coefs <- matrix(NA, length(aliased), 4, dimnames = list(cn, 
                                                              colnames(coefs)))
      coefs[!aliased, ] <- x$coefficients
    }
    printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
                 na.print = "NA", ...)
  }
  cat("\nResidual standard error:", format(signif(x$sigma, 
                                                  digits)), "on", rdf, "degrees of freedom\n")
  if (nzchar(mess <- naprint(x$na.action))) 
    cat("  (", mess, ")\n", sep = "")
  if (!is.null(x$fstatistic)) {
    cat("Multiple R-squared:", formatC(x$r.squared, digits = digits))
    if (adj.r.squared) {
      cat(",\tAdjusted R-squared:", formatC(x$adj.r.squared, 
                                            digits = digits))
    }
    cat("\nF-statistic:", formatC(x$fstatistic[1L], digits = digits), 
        "on", x$fstatistic[2L], "and", x$fstatistic[3L], 
        "DF,  p-value:", format.pval(pf(x$fstatistic[1L], 
                                        x$fstatistic[2L], x$fstatistic[3L], lower.tail = FALSE), 
                                     digits = digits), "\n")
  }
  correl <- x$correlation
  if (!is.null(correl)) {
    p <- NCOL(correl)
    if (p > 1L) {
      cat("\nCorrelation of Coefficients:\n")
      if (is.logical(symbolic.cor) && symbolic.cor) {
        print(symnum(correl, abbr.colnames = NULL))
      }
      else {
        correl <- format(round(correl, 2), nsmall = 2, 
                         digits = digits)
        correl[!lower.tri(correl)] <- ""
        print(correl[-1, -p, drop = FALSE], quote = FALSE)
      }
    }
  }
  cat("\n")
  if(!is.null(x$conf.intervals)){ 
    cat("Coefficient Confidence Intervals\n")
    print(x$conf.intervals)
    cat("\n")
  }
  invisible(x)
} 

Confint <- function(object, ...){
  UseMethod("Confint")
}

Confint.default <- function(object, ...) confint(object, ...)

Confint.lm <- function(object, parm, level = 0.95, vcov.=Vcov, ...) {
    vcov.type <- deparse(substitute(vcov.))
    if(vcov.type == "Boot"){
      return(confint(Boot(object)))
    }
    # the following function borrowed from stats:::format.perc()
    format.perc <- function (probs, digits){
        paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits), "%")
    }
    cf <- coef(object)
    pnames <- names(cf)
    if (missing(parm)) 
      parm <- pnames
    else if (is.numeric(parm)) 
      parm <- pnames[parm]
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    fac <- qt(a, object$df.residual)
    pct <- format.perc(a, 3)
    ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm,  pct))
#    ses <- sqrt(diag(vcov(object)))[parm]
    ses <- sqrt(diag(if(is.matrix(vcov.)) vcov. else vcov.(object)))[parm]
    ci[] <- cf[parm] + ses %o% fac
    ci
  }  
