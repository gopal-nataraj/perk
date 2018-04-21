% script toy2d_nonlin.m
% krr demo in 2d nonlinear toy example containing branch cut
% 2d toy example of krr
%
% gopal nataraj
% university of michigan
%
% version control
%   2017-05-22      original
%   2017-05-23      true.phi now injective 

% irt setup
if ~exist('irtdir', 'var')
  curdir = cd('~/Box/work/irt');
  irtdir = pwd;
  setup();
  cd(curdir);
end

% seed rng
rng(0);

% hyperparameters
n = 1000;
eps = 0.05;
lam0 = 2^5;
t = 100;

% regressands
rho = 3*rand(n,1);
phi = -pi + 2*pi*rand(n,1);

% regressors
x = rho.*cos(phi) + eps*randn(n,1);     % [N] 
y = rho.*sin(phi) + eps*randn(n,1);     % [N]
q = transpose([x y]);                   % [2 N]

% more hyperparameters
Lam = lam0*diag(mean(q,2));             % [2 2]
reg.rho = norm(rho-mean(rho)).^2 / n^2; % [1]
reg.phi = norm(phi-mean(phi)).^2 / n^2; % [1]

% krr coefficients
K = NaN(n,n);
for i = 1:n
  dq = bsxfun(@minus, q(:,i), q);
  dq = Lam \ dq;
  K(i,:) = exp(-(sum(abs(dq.^2),1))/2);
end
a.rho = (K + reg.rho*speye(n)) \ rho;   % [N]
a.phi = (K + reg.phi*speye(n)) \ phi;   % [N]

% simulate noisy test data
test.x = linspace(-3,3,t)';
test.y = linspace(-3,3,t)';
[test.xx, test.yy] = ndgrid(test.x, test.y);
test.xx = test.xx + eps*randn(t);
test.yy = test.yy + eps*randn(t);
test.q = [col(test.xx) col(test.yy)]';  % [2 T^2]

% apply to test data
Kt = NaN(n,t^2);
for i = 1:n
  dq = bsxfun(@minus, q(:,i), test.q);  
  dq = Lam \ dq;                        % [2 T^2]
  Kt(i,:) = exp(-(sum(abs(dq.^2),1))/2);
end
test.rho = reshape(Kt'*a.rho, [t t]);   % [T T]
test.phi = reshape(Kt'*a.phi, [t t]);   % [T T]

% true function values
true.rho = sqrt(test.xx .^2 + test.yy .^2);
true.phi = atan(test.yy./test.xx) + pi*sign(test.yy).*(test.xx<0);

% fig-rho
figure;...
  im(test.x, test.y, cat(3, true.rho, test.rho), [0 3*sqrt(2)], 'cbar', ' ');...
  text(col([0 6]), col([3 3]), col({'$\rho(x,y)$', '$\hat{\rho}(x,y)$'}),...
    'interpreter', 'latex',...
    'FontSize', 16,...
    'VerticalAlignment', 'bottom',...
    'HorizontalAlignment', 'center',...
    'Color', 'k');
figure;...
  hold on;...
  im(test.x, test.y, abs(test.rho-true.rho), [0 0.1], 'cbar', ' ');...
%   scatter(x,y,80,'m');...
  hold off;...
  text(0, 3, {'$|\hat{\rho}(x,y)-\rho(x,y)|$'},...
    'interpreter', 'latex',...
    'FontSize', 16,...
    'VerticalAlignment', 'bottom',...
    'HorizontalAlignment', 'center',...
    'Color', 'k');

% fig-phi
figure;...
  im(test.x, test.y, cat(3, true.phi, test.phi), [-pi pi], 'cbar', ' ');...
  text(col([0 6]), col([3 3]), col({'$\phi(x,y)$', '$\hat{\phi}(x,y)$'}),...
    'interpreter', 'latex',...
    'FontSize', 16,...
    'VerticalAlignment', 'bottom',...
    'HorizontalAlignment', 'center',...
    'Color', 'k');
  
figure;...
  hold on;...
  im(test.x, test.y, abs(test.phi-true.phi), [0 pi/4], 'cbar', ' ');...
%   scatter(x,y,80,'m');...
  hold off;
  text(0, 3, {'$|\hat{\phi}(x,y)-\phi(x,y)|$'},...
    'interpreter', 'latex',...
    'FontSize', 16,...
    'VerticalAlignment', 'bottom',...
    'HorizontalAlignment', 'center',...
    'Color', 'k');
  