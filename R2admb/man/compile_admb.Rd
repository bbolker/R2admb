\name{compile_admb}
\alias{compile_admb}
\alias{read_admb}
\alias{run_admb}
\title{Compile ADMB files, run, read output}
\usage{
  compile_admb(fn,safe=FALSE,re=FALSE, verbose=FALSE,
    admb_errors=c("stop","warn","ignore"))

  run_admb(fn,verbose=FALSE,mcmc=FALSE,
    mcmc.opts=mcmc.control(),profile=FALSE,
    extra.args="",admb_errors=c("stop","warn","ignore"))

  read_admb(fn,verbose=FALSE,profile=FALSE,
    mcmc=FALSE,mcmc.opts=NULL,admbOut=NULL,checkterm=TRUE)
}
\arguments{
  \item{fn}{(character) name of TPL file, without
  extension}

  \item{safe}{(logical) Compile in safe mode?}

  \item{re}{(logical) Compile in random effects (ADMB-RE)
  mode?}

  \item{profile}{(logical) Run likelihood profiles?}

  \item{extra.args}{(character) extra arguments for ADMB
  run}

  \item{mcmc}{(logical) run post-hoc MCMC?}

  \item{mcmc.opts}{options for MCMC run (see
  \code{\link{mcmc.control}})}

  \item{verbose}{(logical) Verbose output?}

  \item{admb_errors}{(character) how to handle
  compilation/linking errors?}

  \item{admbOut}{(character) ADMB run output for inclusion
  in \code{admb} object (for internal use)}

  \item{checkterm}{(logical) compute termination criteria
  (ratio of min/max eigenvalue) and include it in the saved
  object?}
}
\value{
  \itemize{ \item \code{compile_admb} returns nothing (it
  has the side effect of creating an executable) \item
  \code{run_admb} invisibly returns the output produced by
  the ADMB run; it also produces output files on disk as a
  side effect \item \code{read_admb} returns an object of
  class \code{admb}, containing as much information as
  possible gleaned from the output files (parameter
  estimates, standard errors, variance-covariance matrix,
  profiles, MCMC output) }
}
\description{
  With various tests, calls the \code{admb} script to
  compile from a TPL file to an executable, or runs the
  resulting executable
}
\note{
  Compiling also sets executable mode.
}
\author{
  Ben Bolker
}
\keyword{misc}

