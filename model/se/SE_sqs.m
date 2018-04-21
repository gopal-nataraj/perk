%% Computations to find a quadratic surrogate for SE datafit term

syms M0 y;
syms T2 TE a_ex real;

% Signal model 
s_se = 1i .* M0 .* sin(a_ex) .* exp(-TE ./ T2);
f_se = matlabFunction(s_se);

% Datafit term
psi_se = (1/2) .* abs((y - s_se).^2);
psi_se_mf = matlabFunction(psi_se);

% First derivative
dpsi_se = simplify(diff(psi_se, 'T2'));
dpsi_se_mf = matlabFunction(dpsi_se);

% Second derivative
ddpsi_se = simplify(diff(dpsi_se, 'T2'));
ddpsi_se_mf = matlabFunction(ddpsi_se);

% Third derivative
dddpsi_se = simplify(diff(ddpsi_se, 'T2'));
dddpsi_se_mf = matlabFunction(dddpsi_se);

% Fourth derivative
ddddpsi_se = simplify(diff(dddpsi_se, 'T2'));
ddddpsi_se_mf = matlabFunction(ddddpsi_se);


% Plot fourth derivative for various TE values
% If always positive, then Hessian is convex
nTE = 10;
TE = logspace(0, 3, nTE);
T2 = linspace(1, 1000, 1000);
t = round(length(T2)/10);

M0 = 1;
a_ex = pi/2;

for i = 1:nTE
    y = f_se(M0, T2(t), TE(i), a_ex);
    d2psi_se_dat = ddpsi_se_mf(M0, T2, TE(i), a_ex, y);
    d4psi_se_dat = ddddpsi_se_mf(M0, T2, TE(i), a_ex, y);
    
    figure(1); plot(T2, d2psi_se_dat); title('2nd deriv of Psi');
    figure(2); plot(T2, d4psi_se_dat); title('4th deriv of Psi');
    pause(3);
end