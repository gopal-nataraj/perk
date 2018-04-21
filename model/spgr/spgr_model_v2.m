% script spgr_model_v2.m
% deriving a careful spgr signal model
%
% copyright 2016, gopal nataraj, university of michigan
%
% version control
%   1.1     2013-12-26      original
%   2.1     2016-06-08      saving subroutines
%   2.2     2017-06-07      added hessian evaluation

% variable declarations
syms M0 T1 T2 positive;
syms kap dw R2p real;
syms flip TR TE positive;

% assumptions
assume(0 < kap);
assume(0 < R2p);
assume(0 < TE & TE < TR);

% add spatial variation to flip
a = flip * kap;

% explicit calculation, avoiding integrals
E1      = exp(-TR / T1);
tmp     = (1-E1) / (1-E1*cos(a));

ssxy_0  = +1i*M0*sin(a) * tmp;
ssxy_te = ssxy_0 * exp(-TE * (R2p + 1/T2)) * exp(+1i * TE * dw);

% generate signal function
if ~exist('spgr_gen.m', 'file')
  fsxy_te = matlabFunction(ssxy_te,...
    'file', 'spgr_gen.m',...
    'vars', [M0 T1 T2 kap dw R2p flip TR TE]);
end

% row gradient of received signal w.r.t. independent variables, x
ssxy_te_gradx = simplify(jacobian(ssxy_te, [M0 T1 T2]));            % [1 L]

% generate gradient function
if ~exist('spgr_gradx_gen.m', 'file')
  fsxy_te_gradx = matlabFunction(ssxy_te_gradx,...
    'file', 'spgr_gradx_gen.m',...
    'vars', [M0 T1 T2 kap dw R2p flip TR TE]);
end

% hessian of received signal w.r.t. independent variables, x
ssxy_te_hessx = simplify(hessian(ssxy_te, [M0 T1 T2]));             % [L L]
ssxy_te_hessx = reshape(ssxy_te_hessx, [1 3 3]);                    % [1 L L]

% generate hessian function
if ~exist('spgr_hessx_gen.m', 'file')
  fsxy_te_hessx = matlabFunction(ssxy_te_hessx,...
    'file', 'spgr_hessx_gen.m',...
    'vars', [M0 T1 T2 kap dw R2p flip TR TE]);
end
