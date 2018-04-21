function [S] = IR_fun_v1(M02, T1, kap, a_inv, a_ex, TR, TI, TE, wf, is_mag)
%IR_FUN.M Evaluates (RF SE) IR signal responses
%   Input: 
%       M02     (A/m)   M0*exp(-TE/T2), T2-compensated spin density
%       T1      (ms)    Spin-lattice relaxation time
%       kap     ()      Flip angle scaling factor  
%       a_inv   (rad)   Nominal inversion angle (usually pi)
%       a_ex    (rad)   Nominal excitation angle (usually pi/2)
%       TR      (ms)    Repetition time
%       TI      (ms)    Time between inversion and excitation pulses
%       TE      (ms)    Time between excitation pulse and gradient echo
%       wf      (rad)   Median off-resonance frequency (rad/ms)
%       is_mag          Toggle returning magnitude signal on/off
%   Output:
%       S       (A/m)   IR signal
%
% Assumptions:
%   1) Data collection with RF spin echo: thus there is T2 decay
%   2) Adiabatic inversion pulse: thus robust to B1 inhomogeneity
%
% Written by: Gopal Nataraj, Copyright 2015


% Definitions
E1_I = exp(-TI ./ T1);
E1_R = exp(-TR ./ T1);
flip_inv = a_inv;                   % Robust to flip angle variation
flip_ex = kap .* a_ex;              % Affected by flip angle variation
ph = exp(1i * (wf * TE + pi/2));

% Signal
S = M02 .* sin(flip_ex) .* ph .* (1 - (1 - cos(flip_inv)).*E1_I + E1_R);

% Error catch: in case any values are not computable, set to zero
S(isnan(S)) = 0;

% Option to return magnitude signal values only
if (is_mag)
    S = abs(S);
end
end

