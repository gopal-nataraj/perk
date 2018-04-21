  function [sir_grad] = IR_gradx(M0, T1, T2, TR, TI, TE, varargin)
%|function [sir_grad] = IR_gradx(M0, T1, T2, TR, TI, TE, varargin)
%|  
%|  spin-echo inversion recovery signal gradient w.r.t. [M0 T1 T2]
%|    assumes rf spin echo readout at time TI+TE after inversion
%|    assumes perfect (e.g., adiabatic) inversion
%|    allows for imperfect excitation and refocusing
%|
%|  inputs
%|    M0        [(odims)]       spin density
%|    T1        [(odims)]       spin-lattice relaxation time              ms
%|    T2        [(odims)]       spin-spin relaxation time                 ms       
%|    TR        [1]             repetition time                           ms
%|    TI        [1]             delay after inverison, before excitation  ms
%|    TE        [1]             delay after excitation, before readout    ms
%|
%|  options
%|    mask      [(odims)]       object mask                               def: true(odims)
%|    kap       [(odims)]       flip angle scaling                        def: ones(odims)
%|    wf        [(odims)]       off-resonance field map (kHz)             def: zeros(odims)
%|    flip_ex   [1]             nominal excitation flip angle (rad)       def: pi/2
%|    flip_ref  [1]             nominal refocusing flip angle (rad)       def: pi
%|    mag       false|true      toggle (mag) grad of mag signal off|on    def: false
%| 
%|  outputs
%|    sir_grad  [(odims) L]     spin-echo inversion recovery signal gradient at t=TI+TE after RF
%|
%|  copyright 2016, gopal nataraj, university of michigan
%|
%|  version control
%|    1.1       2016-05-31      original
%|    1.2       2016-06-02      changed kap to an optional argument

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

% use generate analytical expressions
sir_grad = IR_shortTR_gradx_gen(M0, T1, T2, kap, TR, TI, TE, arg.flip_ex, arg.flip_ref, wf);

% trick: grad(abs(sir)) w.r.t. [abs(M0) T1 T2] = abs(grad(sir)) w.r.t. [M0 T1 T2]
if arg.mag
  sir_grad = abs(sir_grad);
end

% embed back into original array size
sir_grad = embed(sir_grad, arg.mask);                                     % [(odims) L]
end