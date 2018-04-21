  function [ss_hess] = spgr_hessx(M0, T1, T2, flip, TR, TE, varargin)
%|function [ss_hess] = spgr_hessx(M0, T1, T2, flip, TR, TE, varargin)
%|
%|  spoiled gradient-recalled echo signal hessian
%|    allows for spatial variation in flip angle
%|
%|  inputs
%|    M0        [(odims)]       spin density
%|    T1        [(odims)]       spin-lattice relaxation time              ms
%|    T2        [(odims)]       spin-spin relaxation time                 ms
%|    flip      [1]             nominal nutation angle                    rad
%|    TR        [1]             repetition time                           ms
%|	  TE        [1]             echo time                                 ms	
%|
%|  options
%|    mask      [(odims)]       object mask                               def: true(odims)
%|    kap       [(odims)]       flip angle scaling                        def: ones(odims)
%|    dw        [(odims)]       off-resonance field map (kHz)             def: zeros(odims)
%|    R2p       [(odims)]       broadening linewidth (kHz)                def: zeros(odims)
%|    mag       false|true      toggle magnitude signal off|on            def: false
%|
%|  outputs
%|	  ss_hess   [(odims) L L]   spgr signal hessian at t=TE after RF
%|
%|  copyright 2017, gopal nataraj, university of michigan
%|
%|  version control
%|	  1.1       2017-06-07      original
%|    2.1       2018-01-05      now uses separate subroutine for magnitude signal case
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
if arg.mag
  ss_hess = spgr_hessx_abs_gen(...
    M0, T1, T2, kap, R2p, flip, TR, TE);
else
  ss_hess = spgr_hessx_cmplx_gen(...
    M0, T1, T2, kap, dw, R2p, flip, TR, TE);
end

% protect against unsafe division (0/0 set to 0)
ss_hess(isnan(ss_hess)) = 0;

% set infinities to zero for safety
ss_hess(isinf(ss_hess)) = 0;

% embed back into original array size
ss_hess = embed(ss_hess, arg.mask);                                       % [(odims) L L]
end
