\name{read_plt}
\alias{read_plt}
\title{Read in ADMB profile file}
\usage{
  read_plt(varname)
}
\arguments{
  \item{varname}{(character) Name of profiled variable
  (base name of \code{.plot} file)}
}
\value{
  List containing the following elements
  list(prof=prof1,ci=ci1,prof_norm=profnorm,ci_norm=cinorm)
  \item{prof}{likelihood profile: a two-column matrix
  containing the parameter value and the corresponding
  likelihood (\emph{not} the log-likelihood or negative
  log-likelihood) \item{ci}{matrix of upper and lower
  confidence intervals at the 0.9, 0.95, and 0.975 levels
  \item{prof_norm}{likelihood profile based on a normal
  approximation} \item{cinorm}{confidence interval matrix
  based on normal approximation}
}
\description{
  Read in the output from ADMB likelihood profiling stored
  in a \code{.plt} file
}

