#'Standard accessor functions for ADMB model fits
#'
#'Extract standard information such as log-likelihood, AIC, coefficients, etc.
#'from ADMB model fits
#'
#'
#' @usage \method{AIC}{admb}(object,...,k=2)
#'           \method{vcov}{admb}(object,type=c("par","extra","all"),...)
#'           \method{logLik}{admb}(object,...)
#'           \method{summary}{admb}(object,correlation=FALSE,symbolic.cor = FALSE,...)
#'           \method{stdEr}{admb}(object,type=c("par","extra","all"),...)
#'           \method{print}{admb}(x,verbose=FALSE,...)
#'           \method{coef}{admb}(object,type=c("par","extra","all"),...)
#'           \method{confint}{admb}(object, parm, level=0.95, method="default",...)
#'           \method{deviance}{admb}(object,...)
#' @S3method print admb
#' @S3method print summary.admb
#' @S3method summary admb
#' @S3method AIC admb
#' @S3method vcov admb
#' @S3method logLik admb
#' @S3method stdEr admb
#' @S3method coef admb
#' @S3method confint admb
#' @S3method deviance admb
#'@export stdEr
#'@aliases AIC.admb vcov.admb logLik.admb coef.admb confint.admb deviance.admb
#'stdEr stdEr.admb summary.admb print.admb 
#'@param x an ADMB model fit (of class "admb")
#'@param object an ADMB model fit (of class "admb")
#'@param k penalty value for AIC fits
#'@param type which type of parameters report. "par": parameters only; "extra":
#'sdreport variables; "all": both
#'@param parm (currently ignored: FIXME) select parameters
#'@param level alpha level for confidence interval
#'@param method (character): "default" or "quad", quadratic (Wald) intervals
#'based on approximate standard errors; "profile", profile CIs (if profile was
#'computed); "quantile", CIs based on quantiles of the MCMC-generated posterior
#'density (if MCMC was computed); "HPDinterval", CIs based on highest posterior
#'density (ditto)
#' @param correlation currently unused parameter
#' @param symbolic.cor currently unused parameter
#' @param verbose show messages
#'@param \dots other parameters (for S3 generic compatibility)
#'@return Extracts appropriate values: numeric (scalar) for AIC, type logLik
#'for logLik, numeric vector of coefficients, numeric variance-covariance
#'matrix of parameter estimates
#'@author Ben Bolker
#'@keywords misc
#'@examples
#'
#'  admbex <- system.file("doc","Reedfrog_runs.RData",package="R2admb")
#'  load(admbex)
#'  m1
#'  coef(m1)
#'  summary(m1)
#'  coef(summary(m1)) ## returns just z-table
#'  AIC(m1)
#'  vcov(m1)
#'  logLik(m1)
#'  deviance(m1)
#'
AIC.admb <- function(object,...,k=2) {
	if (length(list(...))>0) stop("multi-object AIC not yet implemented")
	deviance(object)+k*length(coef(object))
}
## confint.default works
confint.admb <- function(object, parm, level=0.95, method="default", ...) {
	if (method %in% c("default","quad")) {
		tab <- confint.default(object)
	} else if (method=="profile") {
		vals <- object[["prof"]]
		if (is.null(vals)) stop("model not fitted with profile=TRUE")
		if (!level %in% c(0.9,0.95,975)) stop("arbitrary levels not yet implemented:",
					"level must be in (0.9,0.95,0.975)")
		tab <- t(sapply(vals,function(x) {
							x$ci[x$ci[,"sig"]==level,c("lower","upper")]
						}))
		colnames(tab) <- paste(c((1-level)/2,(1+level)/2)*100,"%")
		tab
	} else if (method %in% c("quantile","HPDinterval")) {
		vals <- object[["mcmc"]]
		if (is.null(vals)) stop("model not fitted with mcmc=TRUE")
		if (method=="quantile") {
			tab <- t(apply(vals,2,quantile,c((1-level)/2,(1+level)/2)))
		} else {
			require(coda)
			tab <- HPDinterval(as.mcmc(vals))
			colnames(tab) <- paste(c((1-level)/2,(1+level)/2)*100,"%")
		}
	}
	tab
}
print.admb <- function(x, verbose=FALSE, ...) {
	cat("Model file:",x$fn,"\n")
	if (is.null(x$loglik)) {
		cat("No fit\n")
		return(invisible(NULL))
	}
	cat("Negative log-likelihood:",-x$loglik,"\n")
	cat("Coefficients:\n")
	print(coef(x))
	## FIXME: indicate extra parameters?
	if (!is.null(x$mcmc)) {
		mcpar <- attr(x$mcmc,"mcpar")
		cat("MCMC parameters: start=",mcpar[1],", end=",mcpar[2],", thin=",mcpar[3],"\n",sep="")
	}
	if (verbose) cat(x$txt,sep="\n")
}      

summary.admb <- function(object, correlation=FALSE, symbolic.cor = FALSE, ...) {
	p1 <- 1:length(object$coefficients)
	coef.p <- unlist(object$coefficients)
	covmat <- object$vcov
	var.cf <- diag(covmat)
	s.err <- sqrt(var.cf)
	tvalue <- coef.p/s.err
	dn <- c("Estimate", "Std. Error")
	pvalue <- 2 * pnorm(-abs(tvalue))
	coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
	dimnames(coef.table) <- list(names(coef.p), c(dn, 
					"z value", "Pr(>|z|)"))
	ans <- list(coefficients=coef.table,loglik=object$loglik,fn=object$fn)
	class(ans) <- "summary.admb"
	ans
}

print.summary.admb <- function(x,
		digits = max(3, getOption("digits") - 3),
		symbolic.cor = x$symbolic.cor, 
		signif.stars = getOption("show.signif.stars"), ...) {
	coefs <- x$coefficients
	cat("Model file:",x$fn,"\n")
	cat("Negative log-likelihood:",-x$loglik,"\n")
	cat("Coefficients:\n")
	printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
			na.print = "NA", ...)
}


coef.admb <- function(object,type=c("par","extra","all"),...) {
	type <- match.arg(type)
	x <- object$coefficients
	n <- object$npar
	if (is.null(n)) n <- length(x)
	switch(type,par=x[1:n],
			extra=x[-(1:n)],
			all=x)
}
logLik.admb <- function(object,...) object$loglik
vcov.admb <- function(object,type=c("par","extra","all"),...) {
	type <- match.arg(type)
	v <- object$vcov
	n <- object$npar
	if (is.null(n)) n <- ncol(v)
	switch(type,par=v[1:n,1:n,drop=FALSE],
			extra=v[-(1:n),-(1:n),drop=FALSE],
			all=v)
}

stdEr <- function(object, ...) {
	UseMethod("stdEr")
}
stdEr.admb <- function(object,type=c("par","extra","all"),...) {
	type <- match.arg(type)  
	s <- sqrt(diag(object$vcov))
	n <- object$npar
	switch(type,par=s[1:n],
			extra=s[-(1:n)],
			all=s)
}

deviance.admb <- function(object,...) -2*object$loglik
