% script IR_model_v2.m
% deriving a more thorough IR SE model
%
% written by: gopal nataraj
% copyright 2016, university of michigan
%
% version control
%   1.1     2016-02-12      original
%   1.2     2016-05-27      added gradient; generate signal function

% variable definitions
syms M0 T1 T2 kap TR TI TE t positive;
syms a flip_ex flip_ref dw real;

% assumptions
assume((TR > TI) & (TI > TE) & (TE > 0) & (T1 > T2));

% add spatial variation to flip angle
a_ex  = kap * flip_ex;
a_ref = kap * flip_ref;

%% define matrix operations
% excitation about x'
Rx(a) =     [1,             0,              0;...
             0,             cos(a),         sin(a);...
             0,             -sin(a),        cos(a)];
         
% excitation about y'
Ry(a) =     [cos(a),        0,              sin(a);...
             0,             1,              0;...
             -sin(a),       0,              cos(a)]; 
                
% precession and relaxation
E(t) =      [exp(-t/T2),    0,              0;...
             0,             exp(-t/T2),     0;...
             0,             0,              exp(-t/T1)];
Rz(t) =     [cos(dw*t),     sin(dw*t),      0;...
             -sin(dw*t),    cos(dw*t),      0;...
             0,             0,              1];
A(t) = Rz(t) * E(t);

% ideal spoiling
S =         [0,             0,              0;...
             0,             0,              0;...
             0,             0,              1];
         
% mz recovery
b(t) =      [0;             0;              M0*(1-exp(-t/T1))];


%% matrix inverse solution
m0 = (eye(3) + S*A(TR - (TI + TE/2))*Ry(a_ref)*A(TE/2)*Rx(a_ex)*A(TI)) \...
  (S * (A(TR-(TI+TE/2))*Ry(a_ref) * (A(TE/2)*Rx(a_ex)*b(TI) + b(TE/2))...
  + b(TR-(TI+TE/2))));
% pretty(simplify(m0))

mTI_m = -A(TI)*m0 + b(TI);
mTI_p = Rx(a_ex) * mTI_m;
mTI_hTE_m = A(TE/2)*mTI_p + b(TE/2);
mTI_hTE_p = Ry(a_ref) * mTI_hTE_m;
mTI_TE = A(TE/2)*mTI_hTE_p + b(TE/2);
% pretty(simplify(mTI_TE))

% % simple case: ideal excitation/refocusing
% nTI_TE_ideal = simplify(subs(mTI_TE, [a_ex, a_ref], [pi/2, pi]));
% pretty(nTI_TE_ideal)
% 
% 
% %% sanity check: alternate solution
% syms m0_m_x m0_m_y m0_m_z real;
% m0_m = [m0_m_x; m0_m_y; m0_m_z];
% 
% m0_p = -m0_m;
% mTI_m = A(TI)*m0_p + b(TI);
% mTI_p = Rx(a_ex) * mTI_m;
% mTI_hTE_m = A(TE/2)*mTI_p + b(TE/2);
% mTI_hTE_p = Ry(a_ref) * mTI_hTE_m;
% mTI_TE_m = A(TE/2)*mTI_hTE_p + b(TE/2);
% mTI_TE_p = S*mTI_TE_m;
% mTR_m = A(TR-(TI+TE))*mTI_TE_p + b(TR-(TI+TE));
% 
% [B, y] = equationsToMatrix(mTR_m == m0_m, m0_m);
% m0_alt = simplify(linsolve(B, y));
% mTI_TE_alt = simplify(subs(mTI_TE_m, [m0_m_x, m0_m_y, m0_m_z],...
%   [m0_alt(1), m0_alt(2), m0_alt(3)]));
% % pretty(mTI_TE_alt)
% 
% % simple case: ideal excitation/refocusing
% mTI_TE_ideal_alt = simplify(subs(mTI_TE_alt, [a_ex, a_ref], [pi/2, pi]));
% pretty(mTI_TE_ideal_alt)

% received signal
sTI_TE_xy = mTI_TE(1) + 1i*mTI_TE(2);

% generate signal functions
if ~exist('IR_shortTR_gen.m', 'file')
  fTI_TE_xy = matlabFunction(sTI_TE_xy,...
    'file', 'IR_shortTR_gen.m',...
    'vars', [M0 T1 T2 kap TR TI TE flip_ex flip_ref dw]);
end

% row gradient of received signals w.r.t. independent variables, x
sTI_TE_xy_gradx = simplify(jacobian(sTI_TE_xy, [M0 T1 T2]));        % [1 L]

% generate gradient function
if ~exist('IR_shortTR_gradx_gen.m', 'file')
  fTI_TE_xy_gradx = matlabFunction(sTI_TE_xy_gradx,...
    'file', 'IR_shortTR_gradx_gen.m',...
    'vars', [M0 T1 T2 kap TR TI TE flip_ex flip_ref dw]);
end

