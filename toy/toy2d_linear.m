% script toy2d_linear.m
% krr demo in 2d ill-conditioned linear toy example
%
% gopal nataraj
% university of michigan
%
% version control
%   2017-05-23      original

% irt setup
if ~exist('irtdir', 'var')
  curdir = cd('../../../../irt');
  irtdir = pwd;
  setup();
  cd(curdir);
end

% seed rng
rng(0);

% hyperparameters
n = 1000;
eps = 0.01;
lam0 = 2^5;
t = 100;

% regressands
x1 = -3 + 6*rand(n,1);
x2 = -3 + 6*rand(n,1);
x = transpose([x1 x2]);                   % [2 N]

% regressors
del = 1e-1;
A = [1  0; 0.5 del];
q = A*x + eps*randn(2,n);                 % [2 N]

% fisher and crb
Sig = eps^2*speye(2);
F = A'*(Sig\A);
crb = inv(F);
std.x1 = sqrt(crb(1,1));
std.x2 = sqrt(crb(2,2));

% more hyperparameters
Lam = lam0*diag(mean(q,2));               % [2 2]
reg.x1 = norm(x1-mean(x1)).^2 / n^2;      % [1]
reg.x2 = norm(x2-mean(x2)).^2 / n^2;      % [1]

% krr coefficients
K = NaN(n,n);
for i = 1:n
  dq = bsxfun(@minus, q(:,i), q);
  dq = Lam \ dq;
  K(i,:) = exp(-(sum(abs(dq.^2),1))/2);
end
a.x1 = (K + reg.x1*speye(n)) \ x1;        % [N]
a.x2 = (K + reg.x2*speye(n)) \ x2;        % [N]

% simulate noisy test data
test.q1 = linspace(-3,3,t)';
test.q2 = linspace(-3,3,t)';
[test.qq1, test.qq2] = ndgrid(test.q1, test.q2);
test.qq1 = test.qq1 + eps*randn(t);
test.qq2 = test.qq2 + eps*randn(t);
test.q = [col(test.qq1) col(test.qq2)]';  % [2 T^2]

% apply to test data
Kt = NaN(n,t^2);
for i = 1:n
  dq = bsxfun(@minus, q(:,i), test.q);
  dq = Lam \ dq;
  Kt(i,:) = exp(-(sum(abs(dq.^2),1))/2);
end
test.x1 = reshape(Kt'*a.x1, [t t]);
test.x2 = reshape(Kt'*a.x2, [t t]);

% true function values
tmp = A \ test.q;                         % [2 T^2]
true.x1 = reshape(col(tmp(1,:)), [t t]);  
true.x2 = reshape(col(tmp(2,:)), [t t]);  

% fig-x1
figure;...
  im(test.q1, test.q2, cat(3, true.x1, test.x1), 'cbar', ' ');
  text(col([0 6]), col([3 3]), col({'$x_1(y_1,y_2)$', '$\hat{x}_1(y_1,y_2)$'}),...
    'interpreter', 'latex',...
    'FontSize', 16,...
    'VerticalAlignment', 'bottom',...
    'HorizontalAlignment', 'center',...
    'Color', 'k');
figure;...
  hold on;...
  im(test.q1, test.q2, abs(test.x1-true.x1), [0 0.5], 'cbar', ' ');...
  scatter(q(1,:)',q(2,:)',80,'m');...
  hold off;...
  text(0, 3, {'$|\hat{x}_1(y_1,y_2)-x_1(y_1,y_2)|$'},...
    'interpreter', 'latex',...
    'FontSize', 16,...
    'VerticalAlignment', 'bottom',...
    'HorizontalAlignment', 'center',...
    'Color', 'k');

% fig-x2
tmp = minmax(-0.5*test.qq1+test.qq2)/del;
figure;...
  im(test.q1, test.q2, cat(3, true.x2, test.x2), 'cbar', ' ');
  text(col([0 6]), col([3 3]), col({'$x_2(y_1,y_2)$', '$\hat{x}_2(y_1,y_2)$'}),...
    'interpreter', 'latex',...
    'FontSize', 16,...
    'VerticalAlignment', 'bottom',...
    'HorizontalAlignment', 'center',...
    'Color', 'k');
figure;...
  hold on;...
  im(test.q1, test.q2, abs(test.x2-true.x2), [0 10], 'cbar', ' ');...
  scatter(q(1,:)',q(2,:)',80,'m');...
  hold off;...
  text(0, 3, {'$|\hat{x}_2(y_1,y_2)-x_2(y_1,y_2)|$'},...
    'interpreter', 'latex',...
    'FontSize', 16,...
    'VerticalAlignment', 'bottom',...
    'HorizontalAlignment', 'center',...
    'Color', 'k');
  