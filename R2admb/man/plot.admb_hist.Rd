\name{plot.admb_hist}
\alias{plot.admb_hist}
\title{Plot MCMC histogram}
\usage{
  \method{plot}{admb_hist}(x,type=c("lattice","ggplot"),dtype=c("hist","density"),pars,...)
}
\arguments{
  \item{x}{plotting data}

  \item{type}{only "lattice" at present}

  \item{dtype}{either "hist" or "density"}

  \item{pars}{passed to rhist}

  \item{\dots}{additional parameters for compatibility}
}
\value{
  plot object
}
\description{
  Plot MCMC histogram
}

