  function [F, sig_M0, sig_T2] = norm_crlb_SE_2parm_v2(T1, T2, TR, TE, Sig_inv, time_comp)
% function [F, sig_M0, sig_T2] = norm_crlb_SE_2parm_v2(T1, T2, TR, TE, Sig_inv, time_comp)
%
% inputs
%   T1          [1]         spin-lattice relaxation time (ms)
%   T2          [1]         spin-spin relaxation time (ms)
%   TR          [D]         repetition times (ms)
%   TE          [D]         echo times (ms)
%   Sig_inv     [D D]       inverse noise covariance matrix
%   time_comp   0|1         toggle time compensation on/off
%
% outputs
%   F           [P P]       fisher info matrix (assume unit noise variance)
%   sig_M0      [1]         (time-compensated?) CRLB std. dev. for (unity) M0
%   sig_T1      [1]         (time-compensated?) CRLB std. dev. for T2
% 
% Written by: Gopal Nataraj, Copyright 2016
%
% Version control
%   v1.1        2015-05-08  original
%   v2.1        2016-02-24  accounts for finite-TR, incomplete recovery effects

% constant declarations
D = length(TE);
P = 2;                      % number of parameters to estimate

% construct SE signal model gradient columnwise
grad_se = NaN(P, D);
for d = 1:D
    % append SE model gradient for dth parameters
    grad_se(:,d) = ...
        [se(    1, T1, T2, TR(d), TE(d));...
         dse_T2(1, T1, T2, TR(d), TE(d))];
end

% construct fisher information matrix (should be pure real)
F = real(grad_se * Sig_inv * grad_se');

% % Method one: standard deviations using pinv()
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
function s_se = se(M0, T1, T2, TR, TE)
    s_se = M0.*exp(-TE./T2).*(exp(-TR./T1)-exp((TE-TR.*2.0)./T1).*2.0+1.0).*1i;
    s_se(isnan(s_se)) = 0;
end

%% SE first derivative w.r.t. M0
function sse_M0 = dse_M0(~, T1, T2, TR, TE)
    sse_M0 = exp(-TE./T2).*(exp(-TR./T1)-exp((TE-TR.*2.0)./T1).*2.0+1.0).*1i;
    sse_M0(isnan(sse_M0)) = 0;
end

%% SE first derivative w.r.t. T2
function sse_T2 = dse_T2(M0, T1, T2, TR, TE)
    sse_T2 = M0.*1.0./T2.^2.*TE.*exp(-TE./T2).*(exp(-TR./T1)-exp((TE-TR.*2.0)./T1).*2.0+1.0).*1i;
    sse_T2(isnan(sse_T2)) = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%