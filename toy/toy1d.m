% script toy1d.m
% krr demo in 1d
% 
% gopal nataraj
% university of michigan
%
% version control
%   2017-05-22      copied from ~/Box\ Sync/work/proj/mwf/conf/17-ismrm/mat/toy2.m  

% seed rng
rng(0);

% define kernel
sig = 0.2;
kern = @(q1,q2) exp(-(minus(q1,q2).^2) ./ (2*sig^2));

% data
n = 50;
eps = 0.05;
x = -pi/2 + pi*rand(n,1);
y = tan(x)/10 + eps*randn(n,1);

% krr coefficients
rho = 0.2;
[tmp,tmp2] = ndgrid(y,y);
K = kern(tmp,tmp2);
a = (K + rho*speye(n)) \ x;

% apply to test data
t = 1000;
yt = linspace(-2,2,t)';
[tmp,tmp2] = ndgrid(y,yt);
Kt = kern(tmp,tmp2);
xt = Kt' * a;

% figure
figure('color', 'k');

hold on;
plot(yt, atan(10*yt),...
  'color', 'c',...
  'linestyle', '--',...
  'linewidth', 2);
plot(yt, xt,...
  'Color', 'm',...
  'linewidth', 2);
scatter(y,x,80,'w');
for p = 2:10:n
  plot(yt,Kt(p,:)'*a(p),...
  'color','g');
end
hold off;

set(gca,...
  'color', 'k',...
  'xcolor', 'y',...
  'ycolor', 'y',...
  'fontsize', 16);
xlabel('$y$',...
  'interpreter', 'latex',...
  'fontsize', 20);
legend({'$x(y)$', '$\hat{x}(y)$'},...
  'interpreter', 'latex',...
  'textcolor', 'y',...
  'edgecolor', 'y',...
  'location', 'se');

% set(gcf, 'InvertHardcopy', 'off');
% print('-depsc', 'toy-1d.eps');

