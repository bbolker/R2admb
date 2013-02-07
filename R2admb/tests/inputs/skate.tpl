DATA_SECTION
	init_int nobs;
	init_int a;		//age-at-maturity
	init_matrix data(1,nobs,0,3);
	matrix Iobs(1,3,1,nobs);
	ivector iyr(1,nobs);
	LOC_CALCS
		dmatrix tmp=trans(data).sub(1,3);
		Iobs=tmp;
		iyr=ivector(column(data,0));
	END_CALCS
	matrix sim_nu(1,3,1,nobs);
	!! sim_nu=0;
	//!! cout<<data<<endl;

PARAMETER_SECTION
	init_vector logq(1,3,1);
	init_bounded_matrix z(1,3,1,3,0.00,5.0,1);
	//Changed estimates of variance to the inverse to
	//improve mixing properties (a al. WinBugs)
	init_bounded_vector invSDpro(1,3,0.04,800,1);
	init_bounded_vector invSDobs(1,3,0.25,4,1);
	init_bounded_vector theta(1,2,0,1,1);
	init_bounded_number recrate(0.01,1000,1);
	init_vector logNad(1,a,1);
	init_vector logN69(1,2);
	//init_bounded_matrix nu(1,3,1,35,-5,5,2); // No Random effects
	random_effects_matrix nu(1,3,1,35,2); // for Random effects

	objective_function_value f;
	
	vector q(1,3);
	vector N69(1,2);
	vector Nad(1,a);
	vector sdpro(1,3);
	vector sdobs(1,3);
	
	matrix st(1,3,1,35);
	matrix nmod(1,3,1,35);
	matrix logn(1,3,1,35);
	matrix N(1,3,1,35);
	matrix epsilon(1,3,1,35);
	sdreport_number sd_var;
        sdreport_matrix Ipred(1,3,1,35);

PROCEDURE_SECTION

	calculate_process_model2();
	calculate_observation_model();
	calculate_objective_func();
	sd_var=q(1);
	for(int s=1; s<=3; ++s){
	  Ipred(s) = q(s)*N(s);
	}

FUNCTION calculate_process_model2
	int s,y,d; // stage / year / decade
	dvariable Adults;
	q     = mfexp(logq);
	N69   = mfexp(logN69);
	Nad   = mfexp(logNad);
	sdpro = sqrt(1./invSDpro);
	sdobs = sqrt(1./invSDobs);
	N(1,1) = (N69(1)*(1.-theta(1))+0.5*recrate*Nad(1))*exp(-z(1,1)+nu(1,1)*sdpro(1));
	N(2,1) = (N69(2)*(1.-theta(2))+N69(1)*theta(1))*exp(-z(2,1)+nu(2,1)*sdpro(2));
	N(3,1) = (N69(2)*theta(2)+Nad(a))*exp(-z(3,1)+nu(3,1)*sdpro(3));
	
	d=0;
	for(y=1;y<nobs;y++) /* Update d index for each decade */
	{   
		if( !((y-1)%10) && d<3) d++;
		if(y<a) /* Get adults for new recruitment */
		{
			Adults = Nad(y+1);
		}
		else
		{
			Adults = N(3,y-a+1);
		}
		/* Update numbers at each stage */
		N(1,y+1) = (N(1,y)*(1.-theta(1))+0.5*(recrate*Adults)) 
					* exp(-z(1,d)+nu(1,y+1) 
					* sdpro(1)+sim_nu(1,y+1));
		
		N(2,y+1) = (N(2,y)*(1.-theta(2))+N(1,y)*theta(1)) 
					* exp(-z(2,d)+nu(2,y+1) 
					* sdpro(2)+sim_nu(2,y+1));
		
		N(3,y+1) = (N(3,y) + N(2,y)*theta(2)) 
					* exp(-z(3,d)+nu(3,y+1) 
					* sdpro(3)+sim_nu(3,y+1));
	}
	
FUNCTION calculate_observation_model
	int s,y;

	for(s=1;s<=3;s++)
	{
		epsilon(s) = log(Iobs(s)) - log(q(s)*N(s));
	}

FUNCTION calculate_objective_func
	dvar_vector fvec(1,10);
	fvec.initialize();
	
	/* Priors */
	int s;
	for(s=1; s<=3; s++)
	{
		fvec(1)+=dlnorm(q(s),0.0,0.1);
	}
	fvec(2)=dbeta(theta(1),20,80);
	fvec(2)+=dbeta(theta(2),25,75);
	fvec(3)=dlnorm(N69(1),0.,0.4);
	fvec(3)+=dlnorm(N69(2),-1.203973,0.4);
	for(s=1; s<= a;s++)
	{
		fvec(4)+=dlnorm(Nad(s),log(3.5), 0.3);
	}
	fvec(5)=dlnorm(recrate,log(3.25),0.25);

	/* Likelihood of the survey data*/
	for(s=1;s<=3;s++)
	{
		fvec(6) += dnorm(epsilon(s),sdobs(s));
		fvec(8)+=dgamma(invSDobs(s),1.01,1.01);
		fvec(9)+=dgamma(invSDpro(s),1.01,1.01);
	}
	
	/* Random effects */
	for(s=1;s<=3;s++)
	{
		fvec(7) += norm2(nu(s));
	}
	f= sum(fvec);

GLOBALS_SECTION
	#include <statslib/statsLib.h> // '#' preprocess command - put this in before compiling starts
	#undef REPORT
	#define REPORT(object) report << "#" #object "\n" << object << endl;
	
REPORT_SECTION
	REPORT(q);
	//REPORT(sdobs);
	//REPORT(sdpro);
	
	REPORT(Iobs);
	REPORT(N);
	REPORT(Ipred);
	//REPORT(epsilon);
	REPORT(z);
	REPORT(theta);
	REPORT(recrate);

TOP_OF_MAIN_SECTION
	arrmblsize = 50000000;
	gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
	gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
	gradient_structure::set_MAX_NVAR_OFFSET(11158);
	gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);


