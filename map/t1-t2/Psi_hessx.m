  function [hess] = Psi_hessx(y, x, nu, P, W, dim, reg, bool)
%|function [hess] = Psi_hessx(y, x, nu, P, W, dim, reg, bool)
%|
%|  upper bound to rls cost function hessian
%|    
%|  version control
%|    1.1       2016-06-06      original
%|    1.2       2016-06-07      now using pure-real hessian
%|    1.3       2016-09-16      make regularization optional
%|    1.4       2016-10-27      added term involving hess(f)(x)
%|    1.5       2018-01-07      added signal hessian nonnegativity constraint

% first term
rgradf = rgradx_f(x, nu, P, dim, bool);                                 % [DV LV sparse]
hess = abs(rgradf' * (W * rgradf));                                     % [LV LV sparse]

% second term
hessf = hessx_f(x, nu, P, dim, bool);                                   % [DV L^2*V sparse]
hessf = abs(hessf);
tmp = f(x, nu, P, dim, bool) - y;                                       % [DV]
tmp = bsxfun(@max, tmp, eps);
tmp = W * tmp;                                                          
tmp = hessf' * tmp;                                                     % [L^2*V]
tmp = reshape(tmp, [dim.L dim.L dim.V]);                                % [L L V]
tmp = spmat(tmp);                                                       % [LV LV sparse]
hess = hess + tmp;                                                      % [LV LV sparse]


if bool.reg
  hess = hess + r_hessx(reg, x, dim);                                   % [LV LV sparse]
end
end