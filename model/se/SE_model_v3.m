% script SE_model_v3.m
% deriving a more thorough spin-echo model
%
% written by: gopal nataraj
% copyright 2016, university of michigan
%
% version control
%   1.1     2016-02-10      original
%   1.2     2016-02-11      added spoiler
%   1.3     2016-02-17      added matlabFunction()
%   1.4     2016-05-31      added gradient, generate signal function
%   1.5     2017-06-07      added hessian evaluation

% variable definitions
syms M0 T1 T2 kap TR TE t positive;
syms a flip_ex flip_ref dw real;

% assumptions
assume((TR > TE) & (T1 > T2));

% add spatial variation to flip angle
a_ex = kap * flip_ex;
a_ref = kap * flip_ref;

%% define matrix operators
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
m0 = (eye(3) - S*A(TR - TE/2)*Ry(a_ref)*A(TE/2)*Rx(a_ex)) \ ...
  (S * (A(TR - TE/2)*Ry(a_ref)*b(TE/2) + b(TR - TE/2)));
% pretty(simplify(m0))

mTE = A(TE/2) * (Ry(a_ref) * (A(TE/2)*Rx(a_ex)*m0 + b(TE/2))) + b(TE/2);
% pretty(simplify(mTE))

% % simple case: ideal excitation/refocusing
% mTE_ideal = simplify(subs(mTE, [a_ex, a_ref], [pi/2, pi]));
% pretty(mTE_ideal)
% 
% 
% %% sanity check: alternate solution
% syms m0_m_x m0_m_y m0_m_z real;
% m0_m = [m0_m_x; m0_m_y; m0_m_z];
% 
% m0_p = Rx(a_ex)*m0_m;
% mhTE_m = A(TE/2)*m0_p + b(TE/2);
% mhTE_p = Ry(a_ref)*mhTE_m;
% mTE_m = A(TE/2)*mhTE_p + b(TE/2);
% mTE_p = S*mTE_m;
% mTR_m = A(TR-TE)*mTE_p + b(TR-TE);
% 
% [B, y] = equationsToMatrix(mTR_m == m0_m, m0_m);
% m0_alt = simplify(linsolve(B, y));
% mTE_alt = simplify(subs(mTE_m, [m0_m_x, m0_m_y, m0_m_z],...
%   [m0_alt(1), m0_alt(2), m0_alt(3)]));
% % pretty(mTE_alt)
% 
% % simple case: ideal excitation/refocusing
% mTE_ideal_alt = simplify(subs(mTE_alt, [a_ex, a_ref], [pi/2, pi]));
% pretty(mTE_ideal_alt)

% received signal
sTE_xy = mTE(1) + 1i*mTE(2);

% generate signal functions
if ~exist('SE_shortTR_gen.m', 'file')
  fTE_xy = matlabFunction(sTE_xy,...
    'file', 'SE_shortTR_gen.m',...
    'vars', [M0 T1 T2 kap TR TE flip_ex flip_ref dw]);
end

% row gradient of received signals w.r.t. independent variables, x
sTE_xy_gradx = simplify(jacobian(sTE_xy, [M0 T1 T2]));              % [1 L]

% generate gradient function
if ~exist('SE_shortTR_gradx_gen.m', 'file')
  fTE_xy_gradx = matlabFunction(sTE_xy_gradx,...
    'file', 'SE_shortTR_gradx_gen.m',...
    'vars', [M0 T1 T2 kap TR TE flip_ex flip_ref dw]);
end

% hessian of received signals w.r.t. independent variables, x
sTE_xy_hessx = simplify(hessian(sTE_xy, [M0 T1 T2]));               % [L L]
sTE_xy_hessx = reshape(sTE_xy_hessx, [1 3 3]);                      % [1 L L]

% generate hessian function
if ~exist('SE_shortTR_hessx_gen.m', 'file')
  fTE_xy_hessx = matlabFunction(sTE_xy_hessx,...
    'file', 'SE_shortTR_hessx_gen.m',...
    'vars', [M0 T1 T2 kap TR TE flip_ex flip_ref dw]);
end
