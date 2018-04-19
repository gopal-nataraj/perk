  function [hess] = r_hessx(reg, x, dim)
%|function [hess] = r_hessx(reg, x, dim)
%|
%|  evaluates diagonal upper bound of regularization hessian (w.r.t. x) 
%|
%|  inputs
%|    reg       {L cell}        regularizer objects
%|    x         {L cell}        latent object parameters
%|    dim       [1x1 struct]    object containing dimension info  
%|
%|  outputs
%|    hess      [LV LV sparse]  tight upper-bound of regularization penalty hessian
%|
%|  version control
%|    1.1       2016-06-02      original

hess = NaN(dim.L, dim.V);
for l = 1:dim.L
  hess(l,:) = reg{l}.R.denom(reg{l}.R, x{l});                           % [L V]
end

% construct sparse diagonal hessian upper-bound
hess = spdiags(col(hess), 0, dim.L*dim.V, dim.L*dim.V);                 % [LV LV sparse]
end