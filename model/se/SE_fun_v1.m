function [S] = SE_fun_v1(M0, T2, kap, a_ex, TE, is_mag)
%SE_FUN_V1.M Evaluates SE signal response
%   Input:
%       M0      (A/m)   Spin density
%       T2      (ms)    Spin-spin relaxation time
%       kap     ()      Flip angle scaling factor  
%       a_ex    (rad)   Nominal excitation angle (usually pi/2)
%       TE      (ms)    Echo time
%       is_mag          Toggle returning magnitude signal on/off
%   Output:
%       S       (A/m)   Spin-echo signal
%
%   Written by: Gopal Nataraj, Copyright 2014
%   Observe that for spin echo sequences, no off-resonance dephasing

% Signal
flip_ex = kap .* a_ex;
ph = exp(1i * (pi/2));
S = M0 .* sin(flip_ex) .* ph .* exp(-TE ./ T2);

% Option to return magnitude signal values only
if (is_mag) 
    S = abs(S);
end
end

