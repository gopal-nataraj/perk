  function out = stat(samp, varargin)
%|function out = stat(samp, varargin)
%|  
%|  computes sample statistics
%|
%|  inputs
%|    samp      [N]               sample values
%|    
%|  options
%|    true      [1]               true value                            def: empty
%|  
%|  outputs
%|    out       [1x1 struct]      sample statistics object
%|                                always return:
%|     .mean    [1]                 sample mean
%|     .std     [1]                 sample standard deviation
%|     .merr    [1]                 standard error of sample mean
%|     .serr    [1]                 standard error of sample standard deviation
%|                                relative to true value, if passed:
%|     .bias    [1]                 sample bias
%|     .rmse    [1]                 sample root mean-squared error
%|     .rstd    [1]                 relative std dev, normalized by truth
%|  
%|  copyright 2016, gopal nataraj, univeristy of michigan
%|
%|  version control
%|    2016-10-21                  original
%|    2017-07-18                  true value now optional

% default values
arg.true = [];

% substitute in optional values, if provided
arg = vararg_pair(arg, varargin);

% constants
n = length(samp);

tmp = div0(n-1,2);
tmp2 = gammaln(tmp);
tmp3 = gammaln(div0(n,2));
tmp4 = exp(tmp2-tmp3);
k = sqrt(tmp) * tmp4;

tmp = tmp - tmp4^(-2);
v = 2*tmp;

% sample mean and bias
out.mean  = mean(samp);
if ~isempty(arg.true)
  out.bias = out.mean-arg.true;
end

% sample standard deviation
s = std(samp);
out.std = k*s;
if ~isempty(arg.true)
  out.rstd = div0(out.std,arg.true);
end

% root mean squared error
if ~isempty(arg.true)
  out.rmse = sqrt((out.bias).^2 + (out.std).^2);
end

% standard error of sample mean
out.merr  = div0(out.std, sqrt(n));

% standard error of sample standard deviation
tmp = k * sqrt(div0(v,n-1));
out.serr  = out.std * tmp;
end