%% Diffusion Coefficients for SSFP Signal 
% Written by: Gopal Nataraj, 12/11/2014

% Symbolic variables
syms M0 E1m E2m a positive
syms theta m real

% Helper variables
E2m_cos = E2m*cos(m*theta);
E2m_sin = E2m*sin(m*theta);
E1m_cos = E1m*cos(m*theta);
E1m_sin = E1m*sin(m*theta);

% System matrix, when m is nonzero.
A = [E2m_cos - 1,           -E2m_sin,           E2m_cos - 1,        E2m_sin,            0,                      0;...
     E2m_sin*cos(a),        E2m_cos*cos(a) - 1, -E2m_sin*cos(a),    E2m_cos*cos(a) - 1, -2*E1m_cos*sin(a),      2*E1m_sin*sin(a);...
     E2m_sin*sin(a),        E2m_cos*sin(a),     -E2m_sin*sin(a),    E2m_cos*sin(a),     2*(E1m_cos*cos(a) - 1), -2*E1m_sin*cos(a);...
     -E2m_sin,              1 - E2m_cos,        -E2m_sin,           -1+E2m_cos,         0,                      0;...
     -1 + E2m_cos*cos(a),   -E2m_sin*cos(a),    1 - E2m_cos*cos(a), -E2m_sin*cos(a),    2*E1m_sin*sin(a),       2*E1m_cos*sin(a);...
     E2m_cos*sin(a),        -E2m_sin*sin(a),    -E2m_cos*sin(a),    -E2m_sin*sin(a),    -2*E1m_sin*cos(a),      2*(1-E1m_cos*cos(a))];
    
% Solve Ax = b, where b = 0 when m is nonzero.
b = zeros(6, 1);
x = linsolve(A, b);

% Extract coefficients
aR_pm = x(1); % aR(m-1)
aJ_pm = x(2); % aJ(m-1)
aR_mm = x(3); % aR(-m-1)
aJ_mm = x(4); % aJ(-m-1)
bR_pm = x(5); % bR(m)
bR_mm = x(6); % bJ(m)

% Simpler system matrix when m = 0
% However, nonzero b when m = 0
m = 0;
A0 = [E2m,  0,              0;
      0,    E2m*cos(a),     -E1m*sin(a);
      0,    E2m*sin(a),     E1m*cos(a)];
b0 = [0;    M0*sin(a);      M0*(1-cos(a))];
x0 = (A0 - eye(3)) \ b0;
    