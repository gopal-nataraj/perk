  function [F, sig_M02, sig_T1] = norm_crlb_IR_2parm(T1, wf, flip_inv,...
     flip_ex, TR, TI, TE, Sig_inv, time_comp)
% function [F, sig_M02, sig_T1] = norm_crlb_IR_2parm(T1, wf, flip_inv,...
%    flip_ex, TR, TI, TE, Sig_inv, time_comp)
% 
% Inputs:
%   T1          [1]         Spin-lattice relaxation constant
%   wf          [1]         Median off-resonance frequency (rad/ms)
%   flip_inv    [D]         (Compensated) inversion angles (rad)
%   flip_ex     [D]         (Compensated) excitation angles (rad)
%   TR          [D]         Repetition times (ms)
%   TI          [D]         Time between inversion and excitation pulses (ms)
%   TE          [1]         (Fixed) time between excitation pulse and gradient echo (f(ms)
%   Sig_inv     [D D]       Inverse of noise covariance matrix
%   time_comp   [1]         Toggle time compensation on/off
% 
% Outputs:
%   F           [P P]       Fisher info matrix (assume unit noise variance)
%   sig_M02     [1]         (Time-compensated?) CRLB std. dev. for (unity) M02
%   sig_T1      [1]         (Time-compensated?) CRLB std. dev. for T1
%
% Version 1.1:  2015-05-07  Original
%               2015-05-08  Signal model modified
%               2015-06-25  Modified M02 to M02 --> specific for RF SE IR
%
% Written by: Gopal Nataraj

% Gather constant declarations 
D = length(TR);
P = 2;                      % Number of parameters

% Add columns corresponding to IR data
grad_ir = NaN(P, D);
for d = 1:D
    % Append an IR measurement model gradient vector
    grad_ir(:,d) = ...
        [ir(    1, T1, wf, flip_inv(d), flip_ex(d), TR(d), TI(d), TE);...
         dir_T1(1, T1, wf, flip_inv(d), flip_ex(d), TR(d), TI(d), TE)];
end

% Construct Fisher information matrix (should be pure real)
F = real(grad_ir * Sig_inv * grad_ir');

% % % Method one: standard deviations using pinv()
% if (time_comp) 
%     scan_time = sum(TR);
%     sig_all = sqrt(scan_time) * abs(sqrt(diag(pinv(F))));
% else
%     sig_all = abs(sqrt(diag(pinv(F))));
% end
% sig_M02 = sig_all(1);
% sig_T1 =  sig_all(2);

% Method two: standard deviations using condition number
max_cond_num = 1e50;
if (time_comp)
    scan_time = sum(TR);
    sig_M02 = sqrt(scan_time) * abs(sqrt(diag_pinv(F, 1, max_cond_num)));
    sig_T1  = sqrt(scan_time) * abs(sqrt(diag_pinv(F, 2, max_cond_num)));
else
    sig_M02 = abs(sqrt(diag_pinv(F, 1, max_cond_num)));
    sig_T1  = abs(sqrt(diag_pinv(F, 2, max_cond_num)));
end
  end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% IR signal model
function s_ir = ir(M02, T1, wf, a_inv, a_ex, TR, TI, TE)
s_ir = M02.*exp(pi.*5.0e-1i+TE.*wf.*1i).*sin(a_ex).*(exp(-TR./T1)+exp(-TI./T1).*...
    (cos(a_inv)-1.0)+1.0);
s_ir(isnan(s_ir)) = 0;
end

%% IR first derivative w.r.t. M02
function sir_M02 = dir_M02(~, T1, wf, a_inv, a_ex, TR, TI, TE)
sir_M02 = exp(TE.*wf.*1i).*sin(a_ex).*(exp(-TR./T1)+exp(-TI./T1).*(cos(a_inv)-1.0)+1.0).*1i;
sir_M02(isnan(sir_M02)) = 0;
end

%% IR first derivative w.r.t. T1
function sir_T1 = dir_T1(M02, T1, wf, a_inv, a_ex, TR, TI, TE)
sir_T1 = M02.*exp(TE.*wf.*1i).*sin(a_ex).*(1.0./T1.^2.*TR.*exp(-TR./T1)+1.0./T1.^2.*...
    TI.*exp(-TI./T1).*(cos(a_inv)-1.0)).*1i;
sir_T1(isnan(sir_T1)) = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%