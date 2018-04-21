  function [sp_grad, sn_grad] = IR_gradx_v2(M0, T1, T2, TR, TI, TE, varargin)
%|function [sp_grad, sn_grad] = IR_gradx_v2(M0, T1, T2, TR, TI, TE, varargin)
%|  
%|  spin-echo inversion recovery signal gradient w.r.t. [M0 T1 T2 inveff]
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
%|    sp_grad   [(odims) L]     se-ir signal gradient at t=TI+TE after RF about +x'-axis
%|    sn_grad   [(odims) L]     se-ir signal gradient at t=TI+TE after RF about -x'-axis
%|
%|  copyright 2016, gopal nataraj, university of michigan
%|
%|  version control
%|    1.1       2016-05-31      original
%|    1.2       2016-06-02      changed kap to an optional argument
%|    2.1       2016-06-27      allow for imperfect inversion
%|    2.2       2016-06-29      protect against unsafe division

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
  [sp_grad, sn_grad] = IR_gradx_gen_v2(...
    M0, T1, T2, TR, TI, TE, arg.flip_inv, arg.flip_ex, arg.flip_ref, inveff, kap, wf);
  % trick: grad(abs(s*)) w.r.t. [abs(M0) T1 T2 inveff] = abs(grad(s*)) w.r.t. [M0 T1 T2 inveff]
  if arg.mag
    sp_grad = abs(sp_grad);
    sn_grad = abs(sn_grad);
  end
else
  [sp_grad] = IR_gradx_gen_v2(...
    M0, T1, T2, TR, TI, TE, arg.flip_inv, arg.flip_ex, arg.flip_ref, inveff, kap, wf);
  if arg.mag
    sp_grad = abs(sp_grad);
  end
end

% protect against unsafe division
sp_grad(isnan(sp_grad)) = 0;
if nargout>1
  sn_grad(isnan(sn_grad)) = 0;
end

% embed back into original array sizes
sp_grad = embed(sp_grad, arg.mask);       
if nargout>1
  sn_grad = embed(sn_grad, arg.mask);
end
end
