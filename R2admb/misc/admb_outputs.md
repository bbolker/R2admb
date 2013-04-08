ADMB files
===============

Rudimentary documentation on AD Model Builder's I/O

## Input files

* `*.tpl` TPL (primary model definition) file
* `*.dat` input data
* `*.pin` initial parameter values (optional??)

## Output files

* `admodel.cov`  (binary) variance-covariance matrix (fixed parameters only)
* `admodel.hes`  (binary) hessian matrix
* `eigv.rpt`
* `fmin.log`
* `hesscheck`
* `hessian.bin`
* `b[1-9]`, `s[1-9]`
* `variance`
* `*.b[0-9]+`
* `*.bar` (binary) Parameter values
* `*.bgs` 
* `*.cor` Parameter correlations
* `*.eva`
* `*.log`
* `*.luu`
* `*.p[0-9]+`
* `*.par` Parameter values (11: log-likelihood; 16: max gradient)
* `*.plt` Likelihood profile information
* `*.psv` (binary) MCMC output
* `*.r[0-9]+`
* `*.rep` ?? report file
* `*.rhes`
* `*.shess`
* `*.std` Standard errors of parameters
