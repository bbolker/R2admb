\newcommand{\PR}{\Sexpr[results=rd]{tools:::Rd_expr_PR(#1)}}

\name{NEWS}
\title{R2admb News}
\encoding{UTF-8}

\section{Changes in version 0.7.5.3}{
  \itemize{
    \item ignore !! *after* beginning of line in variable definitions
    \item fix mcmc2 bug
  }
}

\section{Changes in version 0.7.5.2}{
  \itemize{
    \item fix Windows bug in run_admb
  }
}

\section{Changes in version 0.7.5.1}{
  \itemize{
    \item export dat_write, add append option, doc tweaks
    \item structure cleanup
  }
}

\section{Changes in version 0.7.5}{
  \itemize{
    \item remove ggplot2 dependency (removed 'ggplot' option from
    plot.admb.hist) [to reduce dependency creep in glmmADMB]
    \item allow for negative-phase in read_par
  }
}

\section{Changes in version 0.7.4}{
  \itemize{
    \item major changes to modularity: exposed compile_admb, run_admb,
    read_admb
    \item changes to clean_admb
    \item improved logic for retrieving profile and MCMC parameter names
    in read_admb; run_admb writes .mcinfo file to disk
    \item write_dat complains if class is not matrix or factor
    \item exposed confint.admb
    \item read_pars reads full std matrix (attempt at
    back-compatibility) [and fixed bug]
    \item fixed bug: std err/estimate names mismatch
    
}}

\section{Changes in version 0.7.2}{
  \itemize{
    \item exposed and documented read_pars, write_pin, write_dat
}}

\section{Changes in version 0.7.1 (2010-08-16)}{
  \itemize{
    \item improved warning handling
    \item took out testing for spurious xxalloc error messages
    \item added ignore_admb_errors
}}

\section{Changes in version 0.7}{
  \itemize{
    \item uncertain: more accessors, general cleanup
}}

\section{Changes in version 0.6 (2010-08-16)}{
  \itemize{
    \item Support for MCMC
    \item Support for profiling
    \item Support for random effects
    \item Expanded vignette
}}

\section{Changes in version 0.5 (2009-11-25)}{
  \itemize{
    \item Initial release.
    Basic stuff works on my system, lightly tested elsewhere.
}}
