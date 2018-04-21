  function [sse] = SE_fun_v4(M0, T1, T2, TR, TE, varargin)
%|function [sse] = SE_fun_v4(M0, T1, T2, TR, TE, varargin) 
%|
%|  spin echo signal model
%|    assumes rf spin echo readout at time TE after excitation
%|    allows for incomplete recovery due to short TR
%|    allows for imperfect excitation and refocusing
%|
%|  inputs
%|    M0        [(odims)]       spin density
%|    T1        [(odims)]       spin-lattice relaxation time              ms
%|    T2        [(odims)]       spin-spin relaxation time                 ms
%|    TR        [1]             repetition time                           ms
%|    TE        [1]             delay after excitation, before readout    ms
%|
%|  options
%|    mask      [(odims)]       object mask                               def: true(odims)
%|    kap       [(odims)]       flip angle scaling                        def: ones(odims)
%|    wf        [(odims)]       off-resonance field map (kHz)             def: zeros(odims)
%|    flip_ex   [1]             nominal excitation flip angle (rad)       def: pi/2
%|    flip_ref  [1]             nominal refocusing flip angle (rad)       def: pi
%|    mag       false|true      toggle magnitude signal off|on            def: false
%| 
%|  outputs
%|    sse       [(odims)]       spin-echo signal at t=TE after RF
%|
%|  copyright 2016, gopal nataraj, university of michigan
%|
%|  version control
%|    1.1       2015-05-08      original (SE_fun_v1.m)
%|    2.1       2016-02-17      accounted for T1 relaxation and finite-TR effects
%|    3.1       2016-02-17      accounted for variation in nominal excitation/refocusing flip
%|    4.1       2016-05-31      repackaged with subfunction format
%|    4.2       2016-06-02      changed kap to an optional argument
%|    4.3       2016-06-29      protect against unsafe division

% default values
arg.mask = [];
arg.kap = [];
arg.wf = [];
arg.flip_ex = pi/2;
arg.flip_ref = pi;
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
if isempty(arg.wf)
  arg.wf = zeros(odims);
end

% vectorize inputs
M0 = masker(M0, arg.mask);
T1 = masker(T1, arg.mask);
T2 = masker(T2, arg.mask);
kap = masker(arg.kap, arg.mask);
wf = masker(arg.wf, arg.mask);

% use generated analytical expressions
sse = SE_shortTR_gen(M0, T1, T2, kap, TR, TE, arg.flip_ex, arg.flip_ref, wf);
if arg.mag
  sse = abs(sse);
end

% protect against unsafe division
sse(isnan(sse)) = 0;

% embed back into original array sizes
sse = embed(sse, arg.mask);
end