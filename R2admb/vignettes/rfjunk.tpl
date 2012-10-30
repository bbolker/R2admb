DATA_SECTION
  init_int nobs                // # of observations
  init_vector nexposed(1,nobs) // # exposed per trial
  init_vector TBL(1,nobs)      // total body length
  init_vector Kill(1,nobs)     // # killed per trial
PARAMETER_SECTION
  // init_bounded_number c(0,1) // baseline mort prob
  // init_bounded_number d(0,50) // size scaling factor
  //  init_bounded_number g(-1,4) // size scaling power
  init_number c
  init_number d
  init_number g
  vector prob(1,nobs)
  objective_function_value f
PROCEDURE_SECTION
  dvariable fpen=0.0;         // penalty variable
  // power-Ricker
  prob = c*pow(elem_prod(TBL/d,exp(1-TBL/d)),g);
  cout << "calc prob\n";
  cout << "c: " << c << " d: " << d << " g: " << g << "\n";
  cout << TBL(1) << " " << pow(TBL(1)/d*exp(1-TBL(1)/d),g) << "\n";
  cout << pow(0.94,g) << "\n";
  cout << pow(0.94,18.0) << "\n";
  cout << prob << "\n";
  // penalties: constrain 0.001 <= prob <= 0.999
  prob = posfun(prob,0.001,fpen);
  f += 1000*fpen;
  prob = 1-posfun(1-prob,0.001,fpen);
  f += 1000*fpen;
  // binomial negative log-likelihood
  f -= sum( log_comb(nexposed,Kill)+
            elem_prod(Kill,log(prob))+
            elem_prod(nexposed-Kill,log(1-prob)));
  cout << "**" << c <<  " " << d << " " << g << " " << fpen << " " << f << "\n";


