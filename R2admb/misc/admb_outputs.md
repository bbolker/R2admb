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
* `*.bar` (binary) Parameter values
* `*.cor` Parameter correlations
* `*.par` Parameter values (11: log-likelihood; 16: max gradient)
* `*.plt` Likelihood profile information
* `*.psv` (binary) MCMC output
* `*.rep` ?? report file
* `*.std` Standard errors of parameters
