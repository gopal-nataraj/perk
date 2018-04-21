% script toy1d_t2est_se.m
% krr comparison to linear regression on log-transformed data in 1d
% 
% gopal nataraj
% university of michigan
%
% version control
%   2017-01-31      original

% irt setup
if ~exist('irtdir', 'var')
  curdir = cd('../../../../irt');
  irtdir = pwd;
  setup();
  cd(curdir);
end

% seed rng
rng(0);

% echo time, in ms
te = 30;
D = length(te);

% bounds
bound.train = col([10 500]);
bound.test = col([eps 700]);
bound.plot.x = [0 1];
bound.plot.y = [0 700];

% perk training parameters
N = 50;
noise.train = 0.01;
lam0 = 2^-1.5;
rho = 2^-20;

% perk testing parameters
T = 100000;
noise.test = 0.01;

% regressands, in ms
x = random('unif', log(bound.train(1)), log(bound.train(2)), [N, 1]);
x = exp(x);                               % [N]

% noisy regressors
tmp = bsxfun(@rdivide, 1, x.');           % [1 N]
s = exp(-te * tmp);                       % [D N]
n = randn(D,N) + 1i*randn(D,N);
n = bsxfun(@times, n, noise.train);
q = bsxfun(@plus, s, n);
q = abs(q);

% perk coefficients
Lam = lam0*diag(mean(q,2));               % [D D]
K = nan(N,N);
for i = 1:N
  dq = bsxfun(@minus, q(:,i), q);
  dq = Lam \ dq;
  K(i,:) = exp(-(sum(abs(dq.^2),1))/2);
end
a = (K + rho*speye(N)) \ x;               % [N]

% noisy test data
xt = linspace(bound.test(1), bound.test(2), T).';
% xt = random('unif', log(bound.test(1)), log(bound.test(2)), [T, 1]);
% xt = exp(xt);
tmp = bsxfun(@rdivide, 1, xt).'; 
st = exp(-te * tmp);                      % [D T]
nt = randn(D,T) + 1i*randn(D,T);         
nt = bsxfun(@times, nt, noise.test);   
qt = st + nt;
qt = abs(qt);                             % [D T]
[qt,idx] = sort(qt);
xt = xt(idx);

% apply perk to test data
Kt = nan(N,T);
for i = 1:N
  dq = bsxfun(@minus, q(:,i), qt);
  dq = Lam \ dq;
  Kt(i,:) = exp(-(sum(abs(dq.^2),1))/2);
end
xhat.perk = Kt'*a;                        % [T]

% apply log-transformed linear regression to test data
tmp = -log(qt); 
num = tmp' .* te;
den = col(sum(abs(tmp.^2),1));
xhat.ltlr = div0(num,den);

% plot 
figure;
hold on;
scatter(qt, xt, '.');
plot(qt, xhat.perk,...
  'linewidth', 2);
plot(qt, xhat.ltlr,...
  'linewidth', 2);
scatter(q, x, 200, 'k');
plot(bound.plot.x, ones(2,1)*bound.train(1), 'k--');
plot(bound.plot.x, ones(2,1)*bound.train(2), 'k--');
hold off;
axis([bound.plot.x bound.plot.y]);
xlabel('$y$',...
  'interpreter', 'latex',...
  'fontsize', 24);
ylabel('ms',...
  'fontsize', 16);
legend({'$T_2$', '$\hat{T}_2^{\mathrm{PERK}}(y)$', '$\hat{T}_2^{\mathrm{MOM}}(y)$'},...
  'interpreter', 'latex',...
  'location', 'nw',...
  'fontsize', 24);
set(gca,...
  'fontsize', 24);
tmp = sprintf('toy-1d,t2-se,n-%u.eps', N);
% print('-depsc', tmp);

% rmse
rmse.perk = sqrt(sum(abs(xhat.perk-xt).^2));
rmse.ltlr = sqrt(sum(abs(xhat.ltlr-xt).^2));
