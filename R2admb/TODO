
Clean-up/documentation:

  * make sure strategy for copying modifying TPL files is robust:
should work reliably even for mixed-case input file names, and should
always survive a crash well.
  e.g. rather than moving the original tpl file out of the way, use
the _gen file, then move _gen results when finished.

  * clean up argument list

  * make sure strategy for guessing/specifying mode of data (esp. numeric/integer)
is robust

  * modularize run: should be able to do any one of
    = write data and input files
    = write modified TPL file based on data/input
    = compile
    = run
    = read results back into R
 
  * implement: phases, 'fixed' argument (set phase < 0); check bounds relative to starting values

  * accessors: confint, coeftab, coefplot

* attach data as well when checking bounds in DATA section

* weaken/eliminate ggplot2 dependency

* comparisons with PBSadmb?

* support for transformations (on the fly)?

* better setup detection, error trapping:
   comments on adcomp/adlink/admb scripts

* test/document examples of random effects in vignette:
  
* ?? formula interface???

* read bar (binary) rather than par file?

* read hessian file?

* evaluate objective function by fixing vals

* finish bounds support: check bounds on initial parameters?

* better strategy for saving originals when changing case/auto-generating
inputs -- temporary filenames, or temp directory??

* TEST check for variables with "." in the name: stop? change to "_" and warn?      

* check for random effects vectors (don't redefine if already there)

* more checks/flags of compiling step -- stop if error

* more info on thinning etc. of MCMC saved/printed in final object


