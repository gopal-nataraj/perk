% Symbolic IR calculations
% Assumption: gradient spoiling on z following inversion to kill any
% remaining transverse magnetization.
%
% Written by: Gopal Nataraj 
% Date created: 2015-05-07
% Date last modified: 2015-06-25    

syms M02 wf;
syms T1 wf a_inv a_ex TR TI TE positive;
wf = 0;

% Definitions
E1_I = exp(-TI ./ T1);
E1_R = exp(-TR ./ T1);
ph = exp(1i * (wf * TE + pi/2));

% Signal model
s_ir = M02 .* sin(a_ex) .* ph .* (1 - (1 - cos(a_inv)).*E1_I + E1_R);
f_ir = matlabFunction(s_ir);

% First derivative w.r.t. M02
ds_ir_M02 = simplify(diff(s_ir, 'M02'));
df_ir_M02 = matlabFunction(ds_ir_M02);

% First derivative w.r.t. T1
ds_ir_T1  = simplify(diff(s_ir, 'T1'));
df_ir_T1  = matlabFunction(ds_ir_T1);

% Second derivative w.r.t. T1
dds_ir_T1 = simplify(diff(ds_ir_T1, 'T1'));
ddf_ir_T1 = matlabFunction(dds_ir_T1);