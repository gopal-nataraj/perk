  function [sp, sm] = dess_fun_v2(M0, T1, T2, flip, TR, TEp, TEm, varargin)
%|function [sp, sm] = dess_fun_v2(M0, T1, T2, flip, TR, TEp, TEm, varargin)
%|
%|  dual-echo steady state signal model
%|    allows for spatial variation in flip angle
%|
%|  inputs
%|    M0        [(odims)]       spin density
%|    T1        [(odims)]       spin-lattice relaxation time              ms
%|    T2        [(odims)]       spin-spin relaxation time                 ms
%|    flip      [1]             nominal nutation angle                    rad
%|    TR        [1]             repetition time                           ms
%|    TEp       [1]             'defocusing' echo time                    ms
%|    TEm       [1]             'refocusing' echo time                    ms
%|
%|  options
%|    mask      [(odims)]       object mask                               def: true(odims)
%|    kap       [(odims)]       flip angle scaling                        def: ones(odims)
%|    dw        [(odims)]       off-resonance field map (kHz)             def: zeros(odims)
%|    R2p       [(odims)]       broadening linewidth (kHz)                def: zeros(odims)
%|    mag       false|true      toggle magnitude signal off|on            def: false
%|
%|  outputs
%|    sp        [(odims)]       'defocused' dess signal at t=TEp after RF
%|    sm        [(odims)]       'refocused' dess signal at t=TEm before RF
%|
%|  copyright 2016, gopal nataraj, university of michigan
%|
%|  version control
%|    1.1       2013-07-19      original
%|    1.2       2014-01-24      used m0prime, m0star instead
%|    1.3       2014-09-12      added flip angle scaling
%|    2.1       2016-06-08      added masking, subfunction format
%|    2.2       2018-01-12      now sets infinities to zero for safety

% default values
arg.mask = [];
arg.kap = [];
arg.dw = [];
arg.R2p = [];
arg.mag = false;

% substitute varargin values as appropriate
arg = vararg_pair(arg, varargin);

% dimensions
odims = size(T1);

% if no mask specified, extrapolate to all voxels
if isempty(arg.mask)
  arg.mask = true(odims);
  N = prod(odims);
else
  N = numel(arg.mask(arg.mask));
end

% if no flip angle scaling specified, set to one
if isempty(arg.kap)
  arg.kap = ones(odims);
end

% if no off-resonance field map specified, set to zero
if isempty(arg.dw)
  arg.dw = zeros(odims);
end

% if no broadening linewidth specified, set to zero
if isempty(arg.R2p)
  arg.R2p = zeros(odims);
end

% vectorize inputs
M0 = masker(M0, arg.mask);
T1 = masker(T1, arg.mask);
T2 = masker(T2, arg.mask);
kap = masker(arg.kap, arg.mask);
dw = masker(arg.dw, arg.mask);
R2p = masker(arg.R2p, arg.mask);

% use generated analytical expressions
sp = dess_echo1_gen(M0, T1, T2, kap, dw, R2p, flip, TR, TEp);
sm = dess_echo2_gen(M0, T1, T2, kap, dw, R2p, flip, TR, TEm);

% protect against unsafe division (0/0 set to 0)
sp(isnan(sp)) = 0;
sm(isnan(sm)) = 0;

% set infinities to zero for safety
sp(isinf(sp)) = 0;
sm(isinf(sm)) = 0;

% optional: return magnitude signal
if arg.mag
  sp = abs(sp);
  sm = abs(sm);
end

% embed back into original array size
sp = embed(sp, arg.mask);
sm = embed(sm, arg.mask);
end