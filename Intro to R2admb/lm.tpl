PARAMETER_SECTION
  vector pred_y(1,nobs);           // model predictions
PROCEDURE_SECTION
  pred_y=xmat*beta;                // linear model
  f=(norm2(pred_y-y));             // objective function to be minimized
  f=nobs/2.*log(f);    
