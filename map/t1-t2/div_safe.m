  function quot = div_safe(dividend, divisor)
%|function quot = div_safe(dividend, divisor)
%|
%|  safe division: 0/0 set to 0
%|  
%|  inputs
%|    dividend  [(dims)]        dividend
%|    divisor   [(dims)]        divisor
%|    
%|  outputs
%|    quot      [(dims)]        safe quotient
%|
%|  version control
%|    1.1       2016-06-05      copied from old code

undef = (dividend == 0) & (divisor == 0);
quot = dividend ./ divisor; 
quot(undef) = 0;
end