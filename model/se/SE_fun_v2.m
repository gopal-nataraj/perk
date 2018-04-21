function [S] = SE_fun_v2(M0, T1, T2, TE, TR, is_mag)
%SE_FUN_V2.M Evaluates SE signal response
%   Input:
%       M0      (A/m)       spin density
%       T1      (ms)        spin-lattice relaxation time
%       T2      (ms)        spin-spin relaxation time
%       TE      (ms)        echo time
%       TR      (ms)        repetition time
%       is_mag              toggle returning magnitude signal on/off
%   Output:
%       S       (A/m)       spin-echo signal
%
%   Written by: Gopal Nataraj, Copyright 2016
%   
%   Changes from v1:
%       more complete model accounts for T1 relaxation, finite TR
%       variation in nominal excitation lumped into M0
%       variation in nominal inversion neglected due to gradient crushers

% signal
ph = exp(1i * (pi/2));
S = M0 .* ph .* exp(-TE ./ T2) .* ...
    (1 - 2*exp((TE - 2*TR) ./ T1) + exp(-TR ./ T1));

% optional: return magnitude signal only
if (is_mag)
    S = abs(S);
end
end