DATA_SECTION
  init_int nsite
  init_int nyear
  init_int ntot
  init_vector shells(1,ntot)
  init_vector year(1,ntot)
  init_vector prev(1,ntot)
  init_vector site(1,ntot)
PARAMETER_SECTION
  init_bounded_number a(0,1)
  init_bounded_number interc(0,1)
  init_bounded_number sdind(0,50)
  init_bounded_number sdsite(0,50)
  init_bounded_number sdyear(0,50)
  objective_function_value f
  sdreport_number ra;    // this & following for MCMC
  sdreport_number rsdind;
  sdreport_number rsdsite;
  sdreport_number rsdyear;
PROCEDURE_SECTION
  ra = a ;
  // power-Ricker
  predshelldens = exp(interc+a*prev+yreff
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
