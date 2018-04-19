  function [n] = wnorm2(x, w, dim)
%|function [n] = wnorm2(x, w, dim)
%|
%|  weighted 2-norm, x'W'Wx
%|
%|  inputs
%|    x         [LV]            latent parameters
%|    w         [L]             weights, in units 1/x
%|    dim       [1x1 struct]    object containing dimension info
%|
%|  outputs
%|    n         [1]             weighted 2-norm, x'W'Wx
%|
%|  version control
%|    1.1       2016-10-12      original

% rescale x
x = reshape(x, [dim.L dim.V]);
x = bsxfun(@div_safe, x, w);
x = col(x);

% norm
n = norm(x,2);
end