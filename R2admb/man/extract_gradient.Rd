\name{extract_gradient}
\alias{extract_gradient}
\title{Extract gradients}
\usage{
  extract_gradient(admbfile)
}
\arguments{
  \item{admbfile}{base name of admb project}
}
\value{
  a dataframe with 3 columns var=variable name, value=
  final parameter value, gradient= gradient value
}
\description{
  Extract gradient values from last iteration of screen
  output and return dataframe with variable names, values
  and gradient, sorted in order of ascending absolute value
  of the gradient.
}
\author{
  Jeff Laake
}

