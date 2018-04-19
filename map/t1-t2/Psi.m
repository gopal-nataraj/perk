  function [cost] = Psi(y, x, nu, P, W, dim, reg, bool)
%|function [cost] = Psi(y, x, nu, P, W, dim, reg, bool)
%|
%|  rls cost function evaluation
%|  
%|  version control
%|    1.1       2016-06-06      original
%|    1.2       2016-09-16      make regularization optional

err = f(x, nu, P, dim, bool) - y;
cost = (1/2) * real(err' * W * err);
if bool.reg
  cost = cost + r(reg, x, dim);
end
end