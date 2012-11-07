DATA_SECTION
   init_int n;                        // number of distances
   init_vector xs(1,n);               // distances
   init_number width;                 // truncation half-width of transect
   init_int ifct;                     // type of detection function 1=hn, 2=haz
   init_int pt;                       // number of parameter types
   init_ivector links(1,pt);          // link number but would rather use strings; not sure how 1=identity, 2=log, 3=logit
   init_ivector k(1,pt);              // vector of number of parameters for each type; cols in design matrix
   init_3darray dm(1,pt,1,n,1,k);     // design matrices - one for each parameter type
PARAMETER_SECTION
   init_matrix beta(1,pt,1,k);        // beta parameters for each parameter type
   matrix parmat(1,pt,1,n);           // matrix of parameter values; 1 to n and 1 to pt types of parameters
   vector par;                        // holds par for an observation likelihood calculation
   number mu;                         // holds single integral fct(x)
   objective_function_value f;        // negative log-likelihood
PROCEDURE_SECTION
   int i, j;
// Create matrix of real parameter values which is pt rows and n columns
   for (i=1;i<=pt;i++)
   {
    parmat(i)=reals(dm(i),beta(i),links(i));
   }
// loop over each observation computing sum of log-likelihood values
   f=0;
   for (j=1;j<=n;j++)
   {
      par=column(parmat,j);
      mu=adromb(&model_parameters::fct,0,width,8);
      f+= -log(fct(xs(j))) + log(mu);
   }  
   
//////////////////////////////   
// Computes reals from betas
//////////////////////////////   
FUNCTION dvar_vector reals(dmatrix& dm, dvar_vector& beta, int ilink)
// dm is the design matrix
// beta is vector of parameters - length macthes ncol(dm)
// ilink is type of link function
    dvar_vector tmp;
	if(ilink==1)
        tmp=dm*beta;
    if(ilink==2)
        tmp=exp(dm*beta);
    if(ilink==3)
        tmp=1/(1+exp(-dm*beta));
    return tmp;
//////////////////////////////////////////////////////////   	
// Computes normalizing constant int fct(x) - 0 to width
//////////////////////////////////////////////////////////   	
FUNCTION dvariable fct(const dvariable& x)
// x is integration variable
// ifct is index for function read from data
   dvariable tmp;
   if(ifct==1)
     tmp=exp(-.5*x*x/(par(1)*par(1)));
   if(ifct==2)
   {
     if(x<0.0000001)
	   tmp=1;
	 else
       tmp=1-exp(-(pow(x/par(1),-(1+par(2)))));
   }
   return tmp;
   