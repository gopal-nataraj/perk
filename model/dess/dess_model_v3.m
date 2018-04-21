% script dess_model_v3.m
% deriving a more thorough dess signal model
%
% copyright 2016, gopal nataraj, university of michigan
%
% version control
%   1.1     2014-04-29      original
%   2.1     2016-06-08      saving subroutines
%   2.2     2017-06-07      added hessian evaluation
%   3.1     2018-01-05      now separately considers real vs. complex cases

% variable declarations
syms M0 T1 T2 positive;
syms kap dw R2p real;
syms flip TR TEp TEm positive;

% assumptions
assume(0 < kap);
assume(0 < R2p);
assume(0 < TEp & TEp < TR);
assume(0 < TEm & TEm < TR);

% add spatial variation to flip
a = flip * kap;

% explicit calculation, avoiding integrals
E1    = exp(-TR / T1);
E2    = exp(-TR / T2);
xi    = (1 - E1*cos(a)) / (E1 - cos(a));
eta   = sqrt((1-E2^2) / (1-E2^2/xi^2));

% complex signal functions
spxy_0 = +1i*M0*tan(a/2) * (1 - eta/xi);
smxy_0 = -1i*M0*tan(a/2) * (1 - eta);
spxy_tep  = spxy_0 * exp(-TEp * (R2p + 1/T2)) * exp(+1i * TEp * dw);
smxy_tem  = smxy_0 * exp(-TEm * (R2p - 1/T2)) * exp(-1i * TEm * dw);
if ~exist('dess_echo1_gen.m', 'file')
  fpxy_tep = matlabFunction(spxy_tep,...
    'file', 'dess_echo1_gen.m',...
    'vars', [M0 T1 T2 kap dw R2p flip TR TEp]);
end
if ~exist('dess_echo2_gen.m', 'file')
  fmxy_tem = matlabFunction(smxy_tem,...
    'file', 'dess_echo2_gen.m',...
    'vars', [M0 T1 T2 kap dw R2p flip TR TEm]);
end

% row gradient of received complex signals w.r.t. independent variables, x
spxy_tep_gradx = simplify(jacobian(spxy_tep, [M0 T1 T2]));          % [1 L]
smxy_tem_gradx = simplify(jacobian(smxy_tem, [M0 T1 T2]));          % [1 L]
if ~exist('dess_echo1_gradx_cmplx_gen.m', 'file')
  fpxy_tep_gradx = matlabFunction(spxy_tep_gradx,...
    'file', 'dess_echo1_gradx_cmplx_gen.m',...
    'vars', [M0 T1 T2 kap dw R2p flip TR TEp]);
end
if ~exist('dess_echo2_gradx_cmplx_gen.m', 'file')
  fmxy_tem_gradx = matlabFunction(smxy_tem_gradx,...
    'file', 'dess_echo2_gradx_cmplx_gen.m',...
    'vars', [M0 T1 T2 kap dw R2p flip TR TEm]);
end

% row gradient of received magnitude signals w.r.t. independent variables, x
mag_spxy_tep_gradx = simplify(jacobian(abs(spxy_tep), [M0 T1 T2])); % [1 L]
mag_smxy_tem_gradx = simplify(jacobian(abs(smxy_tem), [M0 T1 T2])); % [1 L]
if ~exist('dess_echo1_gradx_abs_gen.m', 'file')
  fpxy_tep_gradx = matlabFunction(mag_spxy_tep_gradx,...
    'file', 'dess_echo1_gradx_abs_gen.m',...
    'vars', [M0 T1 T2 kap R2p flip TR TEp]);
end
if ~exist('dess_echo2_gradx_abs_gen.m', 'file')
  fmxy_tem_gradx = matlabFunction(mag_smxy_tem_gradx,...
    'file', 'dess_echo2_gradx_abs_gen.m',...
    'vars', [M0 T1 T2 kap R2p flip TR TEm]);
end

% hessian of received complex signals w.r.t. independent variables, x
spxy_tep_hessx = simplify(hessian(spxy_tep, [M0 T1 T2]));           % [L L]
spxy_tep_hessx = reshape(spxy_tep_hessx, [1 3 3]);                  % [1 L L]
smxy_tem_hessx = simplify(hessian(smxy_tem, [M0 T1 T2]));           % [L L]
smxy_tem_hessx = reshape(smxy_tem_hessx, [1 3 3]);                  % [1 L L]
if ~exist('dess_echo1_hessx_cmplx_gen.m', 'file')
  fpxy_tep_hessx = matlabFunction(spxy_tep_hessx,...
    'file', 'dess_echo1_hessx_cmplx_gen.m',...
    'vars', [M0 T1 T2 kap dw R2p flip TR TEp]);
end
if ~exist('dess_echo2_hessx_cmplx_gen.m', 'file')
  fmxy_tem_hessx = matlabFunction(smxy_tem_hessx,...
    'file', 'dess_echo2_hessx_cmplx_gen.m',...
    'vars', [M0 T1 T2 kap dw R2p flip TR TEm]);
end

% hessian of received magnitude signals w.r.t. independent variables, x
mag_spxy_tep_hessx = simplify(hessian(abs(spxy_tep), [M0 T1 T2]));  % [L L]
mag_spxy_tep_hessx = reshape(mag_spxy_tep_hessx, [1 3 3]);          % [1 L L]
mag_smxy_tem_hessx = simplify(hessian(abs(smxy_tem), [M0 T1 T2]));  % [L L]
mag_smxy_tem_hessx = reshape(mag_smxy_tem_hessx, [1 3 3]);          % [1 L L]
if ~exist('dess_echo1_hessx_abs_gen.m', 'file')
  fpxy_tep_hessx = matlabFunction(mag_spxy_tep_hessx,...
    'file', 'dess_echo1_hessx_abs_gen.m',...
    'vars', [M0 T1 T2 kap R2p flip TR TEp]);
end
if ~exist('dess_echo2_hessx_abs_gen.m', 'file')
  fmxy_tem_hessx = matlabFunction(mag_smxy_tem_hessx,...
    'file', 'dess_echo2_hessx_abs_gen.m',...
    'vars', [M0 T1 T2 kap R2p flip TR TEm]);
end
