  function [cgrad] = r_cgradx(reg, x, dim)
%|function [cgrad] = r_cgradx(reg, x, dim)
%|  
%|  regularization column gradient evaluation w.r.t. latent object parameter x
%|
%|  inputs
%|    reg       {L cell}        regularizer objects
%|    x         {L cell}        latent object parameters
%|    dim       [1x1 struct]    object containing dimension info  
%|
%|  outputs
%|    cgrad     [LV]            column gradient of regularization penalty
%|
%|  version control
%|    1.1       2016-06-02      original

cgrad = NaN(dim.L, dim.V);
for l = 1:dim.L
  cgrad(l,:) = reg{l}.R.cgrad(reg{l}.R, x{l});
end
cgrad = col(cgrad);                                                     % [LV]
end