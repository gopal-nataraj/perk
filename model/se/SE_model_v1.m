% Symbolic SE calculations
% Permitted: possibly imperfect pi/2 excitation pulse.
% Assumption: perfect refocusing pulse (perhaps adiabatic)
% Assumption: TR >> T1 (otherwise really complicated)
% 
% Written by: Gopal Nataraj
% Date created: 2015-05-08
% Date last modified: 2015-05-08

syms M0;
syms T2 TE a_ex;

% Signal model 
s_se = 1i .* M0 .* sin(a_ex) .* exp(-TE ./ T2);
f_se = matlabFunction(s_se);

% First derivative w.r.t. M0
ds_se_M0 = simplify(diff(s_se, 'M0'));
df_se_M0 = matlabFunction(ds_se_M0);

% First derivative w.r.t. T2
ds_se_T2 = simplify(diff(s_se, 'T2'));
df_se_T2 = matlabFunction(ds_se_T2);

% Second derivative w.r.t. T2
dds_se_T2 = simplify(diff(ds_se_T2, 'T2'));
ddf_se_T2 = matlabFunction(dds_se_T2);