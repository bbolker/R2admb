DATA_SECTION
   init_int n;                        // number of distances
   init_vector distances(1,n);        // distances
   init_number width;                 // truncation half-width of transect
   init_int ndetfct;                  // type of detection function 1=hn, 2=haz
   init_int pt;                       // number of parameter types
   init_ivector links(1,pt);          // link number but would rather use strings; not sure how 1=identity, 2=log, 3=logit
   init_ivector one(1,pt);            // necessary dummy vector of 1s
   init_ivector k(1,pt);              // vector of number of parameters for each type; cols in design matrix
   init_3darray dm(1,pt,1,n,one,k);   // design matrices - one for each parameter type
PARAMETER_SECTION
   init_matrix beta(1,pt,one,k);      // beta parameters for each parameter type
   matrix parmat(1,pt,1,n);           // matrix of parameter values; 1 to n and 1 to pt types of parameters
   vector par;                        // holds par for a single integral calculation
   number mu;                         // holds single integral
   objective_function_value f;        // negative log-likelihood
PROCEDURE_SECTION
   int i, j;
   dvariable xx;
   cout << "detfct =" << ndetfct << endl;
   cout << "distances =" << distances << endl;
   cout << "pt =" << pt << endl;
   cout << "k =" << k << endl;
   
// Create matrix of real parameter values which is pt rows and n columns
   for (i=1;i<=pt;i++)
   {
     dvar_matrix xdm=dm(i);                     // convert to dvar_ type
     dvar_vector xbeta=beta(i);                 // convert to dvar_ type 
     cout << "beta " << xbeta << endl;
     cout << "width " << width << endl;
     parmat(i)=reals(xdm,xbeta,links(i));
   }
   f=0;
   for (j=1;j<=n;j++)
   {
      par=column(parmat,j);
      cout << "par " << par << endl;
      mu=adromb(&model_parameters::detfct,0,width,8);
	  dvariable distance=distances(j);
	  cout << distance << endl;
	  xx=detfct(distance);
	  cout << "g(x) " << xx << endl;
      f+= -log(detfct(distance)) + log(mu);
	  cout << "mu " << mu << endl;
	  cout << "f" << f << endl;
   }  
// Computes reals from betas
FUNCTION dvar_vector reals(dvar_matrix& dm, dvar_vector& beta, int link)
    dvar_vector tmp;
	if(link==1)
        tmp=dm*beta;
    if(link==2)
        tmp=exp(dm*beta);
    if(link==3)
        tmp=1/(1+exp(-dm*beta));
    return tmp;
	
// Computes detection function
FUNCTION dvariable detfct(const dvariable& z)
   dvariable tmp;
   if(ndetfct==1)
     tmp=exp(-.5*z*z/(par(1)*par(1)));
   if(ndetfct==2)
   {
     if(z<0.0000001)
	   tmp=1;
	 else
       tmp=1-exp(-(pow(z/par(1),-(1+par(2)))));
   }
   cout << "z = " << z << " tmp = " << tmp;
   return tmp;
   
