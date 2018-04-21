function [S1est, S2est] = dess_fun_M0prime(TR, TE, M0p, wf, a, v1, T2est)
%DESS_FUN_M0PRIME.M Evaluates DESS signals, given a T2 estimate
%   Input:  
%       TR and TE are predefined constants
%       M0p = M0 * exp(-TE/T2') is compensated for broadening
%       wf is the phase accrual from off-resonance
%       a is the flip angle
%       v1 encapsulates all T1 dependence, (1-E1*cos(a))/(E1-cos(a))
%       T2est is the next test value
%   Output: 
%       S1est and S2est are the corresponding signals to T2est

% Definitions
E_2R = exp(-TR./T2est);
E_2E = exp(-TE./T2est);
R = sqrt((1-E_2R.^2) ./ (1-E_2R.^2./v1.^2));
ph = exp(1i * (wf * TE + pi/2));

S1est = M0p .* tan(a/2) .* ph .* E_2E .* (1 - R./v1);
S2est = M0p .* tan(a/2) .* conj(ph) .* (E_2E.^-1) .* (1 - R);

S1est(isnan(S1est)) = 0;
S2est(isnan(S2est)) = 0;