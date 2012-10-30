DATA_SECTION
  init_int nobs                // # of observations
  init_vector nexposed(1,nobs) // # exposed per trial
  init_vector TBL(1,nobs)      // total body length
  init_vector Kill(1,nobs)     // # killed per trial
PARAMETER_SECTION
  init_bounded_number c(0,1) // baseline mort prob
  init_bounded_number d(0,50) // size scaling factor
  init_bounded_number g(-1,25) // size scaling power
  vector prob(1,nobs)    // per capita mort prob
  objective_function_value f
  sdreport_number rc;    // this & following for MCMC
  sdreport_number rd;
  sdreport_number rg;
PROCEDURE_SECTION
  rc = c; rd = d; rg = g;     // set MCMC reporting
  dvariable fpen=0.0;         // penalty variable
  // power-Ricker
  prob = c*pow(elem_prod(TBL/d,exp(1-TBL/d)),g);
  // penalties: constrain 0.001 <= prob <= 0.999
  prob = posfun(prob,0.001,fpen);
  f += 1000*fpen;
  prob = 1-posfun(1-prob,0.001,fpen);
  f += 1000*fpen;
  // binomial negative log-likelihood
  f -= sum( log_comb(nexposed,Kill)+
            elem_prod(Kill,log(prob))+
            elem_prod(nexposed-Kill,log(1-prob)));
