
PARAMETER_SECTION  
  vector ypred(1,NData);

PROCEDURE_SECTION 
  ypred = a + b*x;
  cout << "iteration" << "a= " << a << "b = " << b << ypred << endl;
  obj_fun = norm2(y-ypred);

REPORT_SECTION
 report << ypred << endl;
 report << y << endl;    
