## ADMB to do list

Clean-up/documentation:
-----------
* make sure strategy for copying modifying `tpl` files is robust: should work reliably even for mixed-case input file names, and should always survive a crash well. e.g. rather than moving the original `tpl` file out of the way, use the `_gen` file, then move `_gen` results when finished.
* clean up argument list
* robustify strategy for guessing/specifying mode of data (esp. numeric/integer)
* strategy for ragged arrays vs matrices in `par` files?
* finish bounds support: check bounds on initial parameters?
* migrate this list to "issues" on github?

Misc.
-------------
* implement: phases, 'fixed' argument (set phase < 0); check bounds relative to starting values
* accessors: confint, coeftab, coefplot
* attach data as well when checking bounds in DATA section
* weaken/eliminate `ggplot2` dependency
* compare with `PBSadmb`?
* `write.model` (tricky, need to handle section headings)
* figure out non-standard evaluation to avoid needing variables defined inside and outside of data argument?

Wishlist
--------------
* support for transformations (on the fly)?
* ?? formula interface???
* better setup detection, error trapping: comments on adcomp/adlink/admb scripts
* test/document examples of random effects in vignette:
* read bar (binary) rather than par file?
* evaluate objective function by fixing vals
* better strategy for saving originals when changing case/auto-generating
inputs -- temporary filenames, or temp directory??
* TEST check for variables with "." in the name: stop? change to "_" and warn? 
* check for random effects vectors (don't redefine if already there)
* more checks/flags of compiling step -- stop if error
* more info on thinning etc. of MCMC saved/printed in final object

Done
--------
* modularize run: should be able to do any one of
 * write data and input files
 * write modified TPL file based on data/input
 * compile
 * run
 * read results back into R
* read hessian file
* smarter `do_admb`, figure out whether PARAMETERS and DATA section already exist
