\name{mcmc.control}
\alias{mcmc.control}
\title{Control options for MCMC after ADMB fitting}
\usage{
  mcmc.control(mcmc = 1000, mcmc2 = 0, mcsave,
    mcnoscale = FALSE, mcgrope = FALSE, mcmult = 1,
    mcmcpars = NULL)
}
\arguments{
  \item{mcmc}{Total number of MCMC steps}

  \item{mcmc2}{MCMC2 steps (see ADMB-RE manual)}

  \item{mcsave}{Thinning interval for values saved in the
  PSV file.  Default is \code{pmax(1,floor(mcmc/1000))},
  i.e. aim to save 1000 steps}

  \item{mcnoscale}{don't rescale step size for mcmc
  depending on acceptance rate}

  \item{mcgrope}{(double) Use a candidate distribution that
  is a mixture of a multivariate normal and a fatter-tailed
  distribution with a proportion \code{mcmcgrope} of the
  fatter-tailed distribution; the ADMB manual suggests
  values of \code{mcgrope} between 0.05 and 0.1}

  \item{mcmult}{Multiplier for the MCMC candidate
  distribution}

  \item{mcmcpars}{(character) vector of parameters to track
  in MCMC run. \emph{At least one must be specified.} ADMB
  produces two kinds of output for MCMC.  For any
  \code{sdreport} parameters it will produce a \code{hst}
  file that contains a summary histogram; \code{mcmcpars}
  constructs appropriate \code{sdreport} parameters in the
  auto-generated TPL file.  Step-by-step output for all
  parameters (regulated by \code{mcsave}) is saved in the
  \code{PSV} file.}
}
\value{
  Returns a list of options suitable for passing as the
  \code{mcmc.opts} argument to \code{\link{do_admb}}
}
\description{
  Determines the options (number of steps, save interval,
  etc.)  for running MCMC based on the estimated mode
  (maximum likelihood estimate) and parameter
  variance-covariance matrix
}
\details{
  See the AD Model Builder reference manual. The
  \code{mcrb} option (reduce correlation of the Hessian
  when constructing the candidate distribution) and the
  \code{mcseed} options (seed for random number generator)
  are not yet implemented; \code{mcnoscale} above may not
  work properly
}
\note{
  Some options (\code{mcmc2}, etc.) that can be used in AD
  Model Builder and ADMB-RE may not be available
}
\examples{
mcmc.control(mcmc=2000)
}
\author{
  Ben Bolker
}
\keyword{misc}

