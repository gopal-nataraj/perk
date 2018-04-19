  function [cgrad] = Psi_cgradx(y, x, nu, P, W, dim, reg, bool)
%|function [cgrad] = Psi_cgradx(y, x, nu, P, W, dim, reg, bool)
%|
%|  rls cost function column gradient evaluation
%|    assumes x{1} is complex and x{2:L} are pure real
%|
%|  version control
%|    1.1       2016-06-06      original
%|    1.2       2016-09-16      make regularization optional

% note: sparse * (fatrix2 * vec) = sparse * vec = vec
rgradf = rgradx_f(x, nu, P, dim, bool);                                 % [DV LV sparse]
cgrad = rgradf' * (W * (f(x, nu, P, dim, bool) - y));                   % [LV]

for l = 2:dim.L
  cgrad(l:dim.L:end) = real(cgrad(l:dim.L:end));
end

if bool.reg
  cgrad = cgrad + r_cgradx(reg, x, dim);                                % [LV]
end
end