% script inv_fn.m
% exploring identifiability of inverse spgr/dess signal model
%
% copyright 2017, gopal nataraj, university of michigan
%
% version control
%   2017-07-14        original

% irt
if (~exist('irtdir', 'var'))
  curdir = pwd;
  cd('~/Box/work/irt');
  irtdir = pwd;
  setup();
  cd(curdir);
end

% add relevant directories
addpath('../exp/sim');
addpath('../model/ir');
addpath('../model/se');
addpath('../model/spgr');
addpath('../model/dess');
addpath('../map/sense');
addpath('../map/t1-t2');
addpath('../misc');

% constant declarations
c.sig.sp = 3.8607e-4;
c.sig.de = 3.8607e-4;

% parameter ranges
vec.m0 = linspace(eps,1,20);
vec.t1 = logspace(log10(100),log10(3500),20);
vec.t2 = logspace(log10(10),log10(700),20);

% true parameters
[x.m0, x.t1, x.t2] = ndgrid(vec.m0, vec.t1, vec.t2);
dim.odims = size(x.m0);

% variable acquisition parameters
P.sp.aex  = [15 5]' * pi/180;           % rad
P.sp.tr   = [12.2 12.2]';               % ms
P.de.aex  = (30)' * pi/180;             % rad
P.de.tr   = (17.5)';                    % ms

% fixed acquisition parameters
S.sp      = length(P.sp.aex);
E.sp      = 1;
P.sp.te   = ones(S.sp,E.sp) * 4.67;     % ms

S.de      = length(P.de.aex);
E.de      = 2;
P.de.te   = ones(S.de,E.de) * 4.67;     % ms

% recon options
wght.sp = [1 1]'      * ones(1,E.sp);    
wght.de = (1)'        * ones(1,E.de);
meth.init = {'vpm','perk'};
meth.iter = 'pgpm';
dist.x.t1.prior = 'logunif';
dist.x.t1.supp = col(minmax(vec.t1));
dist.x.t1.nsamp = 300;
dist.x.t2.prior = 'logunif';
dist.x.t2.supp = col(minmax(vec.t2));
dist.x.t2.nsamp = 300;
kmean.C = 1;
rff.snr = 0.15;
rff.std = [];
rff.len = [];
rff.c = 2^0.6;
rff.H = 10^3;
rff.K = 10^5;
inv.perk = 2^-41;
stop.iter = 0;
bool.mag.sp = 1;
bool.mag.de = 1;
bool.chat = 1;
bool.norm = 0;
bool.m0dess = 0;
bool.reset = 1;
bool.rfftst = 0;
bool.nuclip = 1;
bool.reg = 0;
bool.precon = 1;
bool.disp = 0;

% noiseless spgr data
ybar.sp = NaN([dim.odims S.sp E.sp]);
for s = 1:S.sp
  ybar.sp(:,:,:,s,1) = spgr_fun_v2(...
    x.m0, x.t1, x.t2,...
    P.sp.aex(s), P.sp.tr(s), P.sp.te(s,1));
end

% noiseless dess data
ybar.de = NaN([dim.odims S.de E.de]);
for s = 1:S.de
  [ybar.de(:,:,:,s,1), ybar.de(:,:,:,s,2)] = dess_fun_v2(...
    x.m0, x.t1, x.t2,...
    P.de.aex(s), P.de.tr(s), P.de.te(s,1), P.de.te(s,2));
end

% add complex gaussian noise
rng(0);
n.sp = c.sig.sp * (randn(size(ybar.sp)) + 1i*randn(size(ybar.sp)));
n.de = c.sig.de * (randn(size(ybar.de)) + 1i*randn(size(ybar.de)));
y.sp = ybar.sp + n.sp;
y.de = ybar.de + n.de;

% use magnitude data as appropriate
if bool.mag.sp, y.sp = abs(y.sp); end
if bool.mag.de, y.de = abs(y.de); end

% parameter estimation
xhat = struct('init', [], 'iter', []);
t = struct('init', [], 'iter', []);
for p = 1:length(meth.init)
  opt.map = {...
    'wght.sp', wght.sp,...
    'wght.de', wght.de,...
    'meth.init', meth.init{p},...
    'meth.iter', meth.iter,...
    'dist.x.t1', dist.x.t1,...
    'dist.x.t2', dist.x.t2,...
    'kmean.C', kmean.C,...
    'rff', rff,...
    'inv.perk', inv.perk,...
    'stop.iter', stop.iter,...
    'bool', bool...
  };
  [xhat(p), t(p)] = mri_m0t1t2inveff_map(y, P, opt.map{:});
end

% indices of mean values
xf = fieldnames(vec);
for i = 1:numel(xf)
  tmp = minmax(vec.(xf{i}));
  if i==1
    avg.(xf{i}) = mean(tmp);
  else
    avg.(xf{i}) = div0(tmp(2)-tmp(1),log(tmp(2))-log(tmp(1)));
  end
  tmp = abs(vec.(xf{i})-avg.(xf{i}));
  [~,idx.(xf{i})] = min(tmp);
end

% plots
k = 0;
figure;
for c = 1:4
  switch c
    case 1
      tmp = abs(ybar.sp(:,:,:,1,1));
    case 2
      tmp = abs(ybar.sp(:,:,:,2,1));
    case 3
      tmp = abs(ybar.de(:,:,:,1,1));
    case 4
      tmp = abs(ybar.de(:,:,:,1,2));
  end
  
  for r = 1:numel(xf)
    switch r
      case 1
        horiz     = tmp(:,idx.t1,idx.t2);
        vert.ml   = col(xhat(1).init.(xf{r})(:,idx.t1,idx.t2));
        vert.perk  = col(xhat(2).init.(xf{r})(:,idx.t1,idx.t2));
      case 2
        horiz     = tmp(idx.m0,:,idx.t2);
        vert.ml   = col(xhat(1).init.(xf{r})(idx.m0,:,idx.t2));
        vert.perk  = col(xhat(2).init.(xf{r})(idx.m0,:,idx.t2));
      case 3
        horiz     = tmp(idx.m0,idx.t1,:);
        vert.ml   = col(xhat(1).init.(xf{r})(idx.m0,idx.t1,:));
        vert.perk  = col(xhat(2).init.(xf{r})(idx.m0,idx.t2,:));
    end
    vert.tru = col(vec.(xf{r}));
    
    k = k+1;
    subplot(numel(xf),4,c+4*(r-1));
    hold on;
    scatter(horiz,vert.tru);
    scatter(horiz,vert.ml);
    scatter(horiz,vert.perk);
    set(gca, 'ylim', minmax(vert.tru));
    if r==2
      loc = 'ne';
    else
      loc = 'nw';
    end
    legend('truth','grid search','perk','location',loc);
    hold off;
    switch r
      case 1
        ylabel('$\widehat{M}_0$ (a.u.)', 'interpreter', 'latex');
        str = '(M_0, \bar{T}_1, \bar{T}_2)$ (a.u.)';
      case 2
        ylabel('$\widehat{T}_1$ (ms)', 'interpreter', 'latex');
        str = '(\bar{M}_0, T_1, \bar{T}_2)$ (a.u.)';
      case 3
        ylabel('$\widehat{T}_2$ (ms)', 'interpreter', 'latex');
        str = '(\bar{M}_0, \bar{T}_1, T_2)$ (a.u.)';
    end
    xlabel(strcat(sprintf('$\\bar{y}_%u', c), str), 'interpreter', 'latex');
  end
end
% print('-depsc', sprintf('inv-fn,t2min-%ums.eps', min(vec.t2)));