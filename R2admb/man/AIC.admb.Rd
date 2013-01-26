\name{AIC.admb}
\alias{AIC.admb}
\alias{coef.admb}
\alias{confint.admb}
\alias{deviance.admb}
\alias{logLik.admb}
\alias{print.admb}
\alias{stdEr}
\alias{stdEr.admb}
\alias{summary.admb}
\alias{vcov.admb}
\title{Standard accessor functions for ADMB model fits}
\usage{
  \method{AIC}{admb}(object,...,k=2)
    \method{vcov}{admb}(object,type="fixed",...)
    \method{logLik}{admb}(object,...)
    \method{summary}{admb}(object,correlation=FALSE,symbolic.cor
    = FALSE,...)
    \method{stdEr}{admb}(object,type="fixed",...)
    \method{print}{admb}(x,verbose=FALSE,...)
    \method{coef}{admb}(object,type="fixed",...)
    \method{confint}{admb}(object, parm, level=0.95,
    method="default", type="fixed",...)
    \method{deviance}{admb}(object,...)
}
\arguments{
  \item{x}{an ADMB model fit (of class "admb")}

  \item{object}{an ADMB model fit (of class "admb")}

  \item{k}{penalty value for AIC fits}

  \item{type}{which type of parameters to report. Character
  vector, including one or more of "fixed" (standard,
  fixed-effect parameters); "random" (random effect
  parameters); "extra" (report/sdreport variables); "all"
  (all of the above). sdreport variables; "all": both}

  \item{parm}{(currently ignored: FIXME) select parameters}

  \item{level}{alpha level for confidence interval}

  \item{method}{(character): "default" or "quad", quadratic
  (Wald) intervals based on approximate standard errors;
  "profile", profile CIs (if profile was computed);
  "quantile", CIs based on quantiles of the MCMC-generated
  posterior density (if MCMC was computed); "HPDinterval",
  CIs based on highest posterior density (ditto)}

  \item{correlation}{currently unused parameter}

  \item{symbolic.cor}{currently unused parameter}

  \item{verbose}{show messages}

  \item{\dots}{other parameters (for S3 generic
  compatibility)}
}
\value{
  Extracts appropriate values: numeric (scalar) for AIC,
  type logLik for logLik, numeric vector of coefficients,
  numeric variance-covariance matrix of parameter estimates
}
\description{
  Extract standard information such as log-likelihood, AIC,
  coefficients, etc. from ADMB model fits
}
\examples{
admbex <- system.file("doc","Reedfrog_runs.RData",package="R2admb")
 load(admbex)
 m1
 coef(m1)
 summary(m1)
 coef(summary(m1)) ## returns just z-table
 AIC(m1)
 vcov(m1)
 logLik(m1)
 deviance(m1)
}
\author{
  Ben Bolker
}
\keyword{misc}

