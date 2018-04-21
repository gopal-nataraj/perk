  function [sp_hess, sn_hess] = IR_hessx(M0, T1, T2, TR, TI, TE, varargin)
%|function [sp_hess, sn_hess] = IR_hessx(M0, T1, T2, TR, TI, TE, varargin)
%|
%|  spin-echo inversion recovery signal hessian w.r.t. [M0 T1 T2 inveff]
%|    allows for imperfect inversion, excitation and refocusing
%|    assumes perfect spoiling immediately following inversion
%|    assumes rf spin echo readout at time TI+TE after inversion
%|    assumes rf excitations on x'-axis, with alternating phase
%|
%|  inputs
%|    M0        [(odims)]       spin density
%|    T1        [(odims)]       spin-lattice relaxation time              ms
%|    T2        [(odims)]       spin-spin relaxation time                 ms
%|    TR        [1]             repetition time                           ms
%|    TI        [1]             delay after inversion, before excitation  ms
%|    TE        [1]             delay after excitation, before readout    ms
%|
%|  options
%|    mask      [(odims)]       object mask                               def: true(odims)
%|    inveff    [(odims)]       inversion efficiency                      def: ones(odims)
%|    kap       [(odims)]       ex/ref flip angle scaling                 def: ones(odims)
%|    wf        [(odims)]       off-resonance field map (kHz)             def: zeros(odims)
%|    flip_inv  [1]             nominal inversion flip angle (rad)        def: pi
%|    flip_ex   [1]             nominal excitation flip angle (rad)       def: pi/2
%|    flip_ref  [1]             nominal refocusing flip angle (rad)       def: pi
%|    mag       false|true      toggle magnitude signal off|on            def: false
%| 
%|  outputs
%|    sp_hess   [(odims) L L]   se-ir signal hessian at t=TI+TE after RF about +x'-axis
%|    sn_hess   [(odims) L L]   se-ir signal hessian at t=TI+TE after RF about -x'-axis
%|
%|  copyright 2017, gopal nataraj, university of michigan
%|
%|  version control
%|    1.1       2017-06-07      original

% default values
arg.mask = [];
arg.inveff = [];
arg.kap = [];
arg.wf = [];
arg.flip_inv = pi;
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

% if no inversion efficiency specified, assume perfect (one)
if isempty(arg.inveff)
  arg.inveff = ones(odims);
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
inveff = masker(arg.inveff, arg.mask);
kap = masker(arg.kap, arg.mask);
wf = masker(arg.wf, arg.mask);

% use generated analytical expressions
if nargout>1
  [sp_hess, sn_hess] = IR_hessx_gen(...
    M0, T1, T2, TR, TI, TE, arg.flip_inv, arg.flip_ex, arg.flip_ref, inveff, kap, wf);
else
  [sp_hess] = IR_hessx_gen(...
    M0, T1, T2, TR, TI, TE, arg.flip_inv, arg.flip_ex, arg.flip_ref, inveff, kap, wf);
end

% trick: hess(abs(s*)) w.r.t. [abs(M0 ...] = abs(hess(s*)) w.r.t. [M0 ...]
sp_hess = abs(sp_hess);
if nargout>1
  sn_hess = abs(sn_hess);
end

% protect against unsafe division
sp_hess(isnan(sp_hess)) = 0;
if nargout>1
  sn_hess = abs(sn_hess);
end

% embed back into original array size
sp_hess = embed(sp_hess, arg.mask);       
if nargout>1
  sn_hess = embed(sn_hess, arg.mask);
end
  end