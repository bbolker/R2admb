PARAMETER_SECTION
  vector prob(1,nobs)    // per capita mort prob

PROCEDURE_SECTION

  // power-Ricker
  prob = c*pow(elem_prod(TBL/d,exp(1-TBL/d)),g);
  // binomial negative log-likelihood
  f -= sum( log_comb(nexposed,Kill)+
            elem_prod(Kill,log(prob))+
            elem_prod(nexposed-Kill,log(1-prob)));
