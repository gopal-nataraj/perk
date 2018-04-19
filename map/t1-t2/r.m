  function [rval] = r(reg, x, dim)
%|function [rval] = r(reg, x, dim)
%|  
%|  regularization function evaluation
%|
%|  inputs
%|    reg       {L cell}        regularizer objects
%|    x         {L cell}        latent object parameters
%|    dim       [1x1 struct]    object containing dimension info  
%|
%|  outputs
%|    rval      [1]             regularization penalty value
%|
%|  version control
%|    1.1       2016-06-02      original

rval = 0;
for l = 1:dim.L
  rval = rval + reg{l}.R.penal(reg{l}.R, x{l});
end
end