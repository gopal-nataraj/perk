  function [F, sig_M0, sig_T2] = norm_crlb_SE_2parm(T2, flip_ex, TR, TE, Sig_inv, time_comp)
% function [F, sig_M0, sig_T2] = norm_crlb_SE_2parm(T2, flip_ex, TR, TE, Sig_inv, time_comp)
% 
% Inputs:
%   T2          [1]         Spin-spin relaxation constant (ms)
%   flip_ex     [D]         (Compensated) excitation angles (rad)
%   TR          [D]         Repetition times (ms)
%   TE          [D]         Echo times (ms)
%   Sig_inv     [D D]       Inverse of noise covariance matrix
%   time_comp   [1]         Toggle time compensation on/off
% 
% Outputs: 
%   F           [P P]       Fisher info matrix (assume unit noise variance)
%   sig_M0      [1]         (Time-compensated?) CRLB std. dev. for (unity) M0
%   sig_T1      [1]         (Time-compensated?) CRLB std. dev. for T2
% 
% This CRLB assumes a model that has TR>>T1 for simplicity.
% This CRLB assumes a model that has a perfect 180-deg inversion pulse.
%
% Written by: Gopal Nataraj
% Date created: 2015-05-08
% Date last modified: 2015-05-08

% Gather constant declarations
D = length(TE);
P = 2;                      % Number of parameters

% Add columns corresponding to SE data
grad_se = NaN(P, D);
for d = 1:D
    % Append a SE measurmeent model gradient vector
    grad_se(:,d) = ...
        [se(    1, T2, flip_ex(d), TE(d));...
         dse_T2(1, T2, flip_ex(d), TE(d))];
end

% Construct Fisher information matrix (should be pure real)
F = real(grad_se * Sig_inv * grad_se');

% % % Method one: standard deviations using pinv()
% if (time_comp) 
%     scan_time = sum(TR);
%     sig_all = sqrt(scan_time) * abs(sqrt(diag(pinv(F))));
% else
%     sig_all = abs(sqrt(diag(pinv(F))));
% end
% sig_M0 = sig_all(1);
% sig_T2 = sig_all(2);

% Method two: standard deviations using condition number
max_cond_num = 1e50;
if (time_comp)
    scan_time = sum(TR);
    sig_M0 = sqrt(scan_time) * abs(sqrt(diag_pinv(F, 1, max_cond_num)));
    sig_T2 = sqrt(scan_time) * abs(sqrt(diag_pinv(F, 2, max_cond_num)));
else
    sig_M0 = abs(sqrt(diag_pinv(F, 1, max_cond_num)));
    sig_T2 = abs(sqrt(diag_pinv(F, 2, max_cond_num)));
end
  end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SE signal model
function s_se = se(M0, T2, a_ex, TE)
s_se = M0.*exp(-TE./T2).*sin(a_ex).*1i;
s_se(isnan(s_se)) = 0;
end

%% SE first derivative w.r.t. M0
function sse_M0 = dse_M0(~, T2, a_ex, TE)
sse_M0 = exp(-TE./T2).*sin(a_ex).*1i;
sse_M0(isnan(sse_M0)) = 0;
end

%% SE first derivative w.r.t. T2
function sse_T2 = dse_T2(M0, T2, a_ex, TE)
sse_T2 = M0.*1.0./T2.^2.*TE.*exp(-TE./T2).*sin(a_ex).*1i;
sse_T2(isnan(sse_T2)) = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%