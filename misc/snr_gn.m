  function snr = snr_gn(y, n, mask, varargin)
%|function snr_gn(y, n)
%|
%|  computes signal-to-noise ratio
%|    note: for matlab 2013b and later, could use built-in snr(...)
%|
%|  inputs
%|    y       [(odims) (ddims)] noisy data
%|    n       [(odims) (ddims)] background noise estimate
%|    mask    [(odims)]         mask over which to compute snr
%|
%|  options
%|    unit    {1}               {'amp', 'dB', 'log10', 'log2'}          def: 'dB'
%|
%|  outputs
%|    snr     [(ddims)]         snr estimates
%|
%|  copyright 2016, gopal nataraj, university of michigan
%|
%|  version control
%|    1.1     2016-06-09    original

% default options
arg.unit = 'dB';

% substitute varargin values as appropriate
arg = vararg_pair(arg, varargin);

% masking 
y = masker(y, mask);                                                    % [V (ddims)]
n = masker(n, mask);                                                    % [V (ddims)]

% compute snr
num = squeeze(sqrt(sum(abs(y.^2), 1)));                                          % [1 (ddims)]
den = squeeze(sqrt(sum(abs(n.^2), 1)));                                          % [1 (ddims)]

switch arg.unit
  case 'amp'
    snr = num ./ den;
  case 'dB'
    snr = 20 * log10(num ./ den);
  case 'log10'
    snr = log10(num ./ den);
  case 'log2'
    snr = log2(num ./ den);
  otherwise
    error('Unknown unit requested?!');
end
end
