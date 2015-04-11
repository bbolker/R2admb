DATA_SECTION

  init_matrix X(1,56,1,4)
  init_matrix Zherd(1,56,1,15)
  init_vector incidence(1,56)
  init_vector size(1,56)
  init_int nobs

PARAMETER_SECTION

  objective_function_value f
  init_vector beta(1,4)
  init_bounded_number sigma_herd(1e-04,20)
  vector herdvec(1,nobs)
  vector eta(1,nobs)
  vector mu(1,nobs)

  random_effects_vector u_herd(1,15)
PROCEDURE_SECTION

  herdvec = sigma_herd*(Zherd*u_herd);
  eta = X*beta;                      // form linear predictor 
  eta += herdvec;                    // augment with random effects
  mu = pow(1.0+exp(-eta),-1.0);      // logistic transform
  // binomial log-likelihood (unnormalized)
  f -= sum(elem_prod(incidence,log(mu))+
        elem_prod(size-incidence,log(1.0-mu)));
  
  f+=0.5*norm2(u_herd);  // log-prior (standard normal)
