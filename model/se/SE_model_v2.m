% script SE_calc_v2.m
% symbolic SE calculations
% 
% Written by: Gopal Nataraj, Copyright 2016
%
% Version control
%   v1.1        2015-05-08      assumes perfect refocusing, TR >> T1
%   v2.1        2016-02-24      now accounts for T1 relaxation, finite TR
%                               variation in nominal excitation lumped into M0
%                               variation in nominal inversion neglected due to gradient crushers

syms M0;
syms T1 T2 TR TE positive;

s_se = 1i .* M0 .* exp(-TE ./ T2) .* (1 - 2*exp((TE - 2*TR) ./ T1) + exp(-TR ./ T1));
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