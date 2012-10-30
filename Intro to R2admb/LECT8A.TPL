DATA_SECTION
  init_int NData;
  init_vector x(1,NData);
  init_vector y(1,NData);

PARAMETER_SECTION  
  init_number a;
  init_number b;

  vector ypred(1,NData);
  objective_function_value obj_fun;

PROCEDURE_SECTION
  
  ypred = a + b*x;
  cout << "iteration" << "a= " << a << "b = " << b << ypred << endl;
  
  obj_fun = norm2(y-ypred);

REPORT_SECTION
 report << ypred << endl;
 report << y << endl;    
