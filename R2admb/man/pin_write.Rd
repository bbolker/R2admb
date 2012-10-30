\name{write_pin}
\alias{write_pin}
\alias{write_dat}
\alias{dat_write}
\title{
  Write parameter and data files for ADMB
}
\description{
  Given base filenames and lists, write output files for starting
  parameter values and data in a format suitable for input by
  AD Model Builder
}
\usage{
write_pin(name, L)
write_dat(name, L, append=FALSE)
dat_write(name, L, append=FALSE)
}
\arguments{
  \item{name}{(character) the base name of the file}
  \item{L}{a list of objects to be written to file}
  \item{append}{(logical) append to existing file?}
}
\value{
  Returns nothing; creates files in the current working directory as a
  side effect
}
\note{numeric vectors and matrices are the only objects that can be
  written (at present)}
\author{
  Hans Skaug
}
\seealso{
  \code{\link{read_pars}}
}
\keyword{misc}
