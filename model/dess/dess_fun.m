function [S1, S2] = dess_fun(M0s, T1, T2, kap, a, TR, TE, wf, is_mag)
%DESS_FUN.M Evaluates DESS signal responses
%   Input: 
%       M0s     (A/m)   M0*exp(-TE/T2s), T2*-compensated spin density
%       T1      (ms)    Spin-lattice relaxation time
%       T2      (ms)    Spin-spin relaxation time
%       kap     ()      Flip angle scaling factor  
%       a       (rad)   Nominal nutation angle
%       TR      (ms)    Repetition time
%       TE      (ms)    Echo time
%       wf      (rad)   Median off-resonance
%       is_mag          Toggle returning magnitude signal on/off
%   Output:
%       S1      (A/m)   SSFP-FID signal
%       S2      (A/m)   SSFP-Echo signal
%
%   Written by: Gopal Nataraj, Copyright 2014

% Definitions
E1_R = exp(-TR ./ T1);
E2_R = exp(-TR ./ T2);
E2_E = exp(-TE ./ T2);
flip = kap .* a;

% Intermediate variables
v1 = (1 - E1_R.*cos(flip)) ./ (E1_R - cos(flip));
R = sqrt((1-E2_R.^2) ./ (1-E2_R.^2./v1.^2));
ph = exp(1i * (wf * TE + pi/2));

% Signals
S1 = M0s .* tan(flip/2) .* ph .* (1 - R./v1);
S2 = M0s .* tan(flip/2) .* conj(ph) .* (E2_E.^-2) .* (1 - R);

% Error catch: in case any values are not computable, set to zero
S1(isnan(S1)) = 0;
S2(isnan(S2)) = 0;

% Option to return magnitude signal values only
if (is_mag)
    S1 = abs(S1);
    S2 = abs(S2);
end
end