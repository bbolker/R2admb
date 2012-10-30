\name{read_pars}
\alias{read_pars}
\alias{read_psv}
\title{Read in parameters from an AD Model Builder run}
\usage{
  read_pars(fn,drop_phase=TRUE)

  read_psv(fn,names=NULL)
}
\arguments{
  \item{fn}{(character) Base name of AD Model Builder}

  \item{names}{(character) Names of variables}

  \item{drop_phase}{(logical) drop negative-phase (fixed)
  parameters from results?}
}
\value{
  List containing the following elements 1) coefficients:
  parameters estimates, 2) coeflist parameter estimates in
  list format, with proper shape (vectors, matrices, etc.)
  3) se estimated standard errors of coefficients 4) loglik
  log-likelihood 5) maxgrad maximum gradient of
  log-likelihood surface 6) Object cor correlation matrix
  7) vcov variance-covariance matrix 8) npar number of
  parameters
}
\description{
  Reads coefficients, standard errors, log-likelihoods,
  maximum gradients, correlation and variance-covariance
  matrices from AD Model Builder output files
}
\details{
  Given the output from an ADMB run on FOO.tpl,
  \code{read_pars} reads the files FOO.par (parameters,
  log-likelihood, max gradient); FOO.std (standard
  deviations); and FOO.cor (correlations).  \code{read_psv}
  reads the output of MCMC runs
}
\section{Warning}{
  The \code{coeflist} component is untested for data
  structures more complicated than scalars, vectors or
  matrices (i.e. higher-dimensional or ragged arrays)
}
\author{
  Ben Bolker
}
\seealso{
  \code{\link{write_pin}}, \code{\link{write_dat}}
}
\keyword{misc}

