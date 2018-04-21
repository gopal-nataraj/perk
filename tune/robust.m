% script robust.m
% exploring sensitivity to training noise std dev rff.std
%
% copyright 2017, gopal nataraj, university of michigan
%
% version control
%   2018-02-06      adapted from holdout.m

% irt 
if (~exist('irtdir', 'var'))
  curdir = cd('~/Box/work/irt'); 
  irtdir = pwd;
  setup();
  cd(curdir);
end

% mapping
addpath('../model/ir');
addpath('../model/se');
addpath('../model/spgr');
addpath('../model/dess');
addpath('../map/sense');
addpath('../map/t1-t2');
addpath('../misc');

% noise standard deviation
c.sig = 3.8607e-4; 

% latent parameter distribution
dist.x.m0.supp    = [0.01 1];
dist.x.m0.prior   = 'unif';
dist.x.t1.supp    = [400 2000].';
dist.x.t1.prior   = 'logunif';
dist.x.t2.supp    = [40 200].';
dist.x.t2.prior   = 'logunif';

% known parameter distribution
dist.nu.kap.supp  = [0.5 2];
dist.nu.kap.prior = 'lognormal';      % log(nu) is normal
dist.nu.kap.logmu = 0;                % log(nu) mean
dist.nu.kap.logsd = 0.2;              % log(nu) std dev
dist.nu.b0.supp   = [0 0];
dist.nu.b0.prior  = 'unif';
dist.nu.r2p.supp  = [0 0];
dist.nu.r2p.prior = 'unif';

% header options
header.br = 0;                        % use broad range of parameters
header.tru = 0;                       % initialize with truth (to test)
header.sv = 0;                        % save holdout results
header.im = 1;                        % show holdout plots
header.pr = 1;                        % print holdout plots

% holdout options
if header.br
  hld.sig = 2.^col(linspace(-40,10,51));
else
  hld.sig = 2.^col(linspace(-15,-5,51));
end
hld.V = 10^5;

% sample latent parameters
field = fieldnames(dist.x);
for l = 1:length(field)
  tmp = dist.x.(field{l}).supp;
  switch dist.x.(field{l}).prior
    case 'unif'
      x.(field{l}) = random('unif', tmp(1), tmp(2), [hld.V 1]);
      
    case 'logunif'
      tmp = log(tmp);
      tmp2 = random('unif', tmp(1), tmp(2), [hld.V 1]);
      x.(field{l}) = exp(tmp2);
      
    otherwise
      error('unknown dist for latent parameter %u.', l);
  end
  x.(field{l}) = reshape(x.(field{l}), 100, []);
end

% sample known parameters
field = fieldnames(dist.nu);
for k = 1:length(field)
  tmp = dist.nu.(field{k}).supp;
  switch dist.nu.(field{k}).prior
    case 'unif'
      nu.(field{k}) = random('unif', tmp(1), tmp(2), [hld.V 1]);
      
    case 'logunif'
      tmp = log(tmp);
      tmp2 = random('unif', tmp(1), tmp(2), [hld.V 1]);
      nu.(field{k}) = exp(tmp2);
      
    case 'lognormal'
      tmp2 = makedist('normal',...
        'mu', dist.nu.(field{k}).logmu,...
        'sigma', dist.nu.(field{k}).logsd);
      tmp2 = truncate(tmp2, log(tmp(1)), log(tmp(2)));
      tmp2 = random(tmp2, [hld.V 1]);
      nu.(field{k}) = exp(tmp2);
      
    otherwise
      error('unknown dist for known parameter %u.', k);
  end
  nu.(field{k}) = reshape(nu.(field{k}), 100, []);
end

P.sp.aex  = [15 5]' * pi/180;           % rad
P.sp.tr   = [12.2 12.2]';               % ms
P.de.aex  = (30)' * pi/180;             % rad
P.de.tr   = (17.5)';                    % ms

S.sp      = length(P.sp.aex);
E.sp      = 1;
P.sp.te   = ones(S.sp,E.sp) * 4.67;     % ms

S.de      = length(P.de.aex);
E.de      = 2;
P.de.te   = ones(S.de,E.de) * 4.67;     % ms

% recon options
wght.sp = [1 1]'  * ones(1,E.sp);    
wght.de = (1)'    * ones(1,E.de);
meth.init = 'perk';
meth.iter = 'pgpm';
rff.snr = 0.15;
rff.std = [];
rff.len = [];
rff.c = 2^0.6;
rff.H = 10^3;
rff.K = 10^6;
inv.perk = 2^-41;
stop.iter = 0;
bool.mag.ir = 0;
bool.mag.se = 0;
bool.mag.sp = 1;
bool.mag.de = 1;
bool.chat = 1;
bool.norm = 1;%0
bool.reset = 1;
bool.rfftst = 0;
bool.nuclip = 1;
bool.reg = 0;
bool.precon = 0;
bool.disp = 0;

% acquire noiseless spgr data
y.sp = NaN([size(x.m0), S.sp, E.sp]);
for s = 1:S.sp
  y.sp(:,:,s,1) = spgr_fun_v2(...
    x.m0, x.t1, x.t2,...
    P.sp.aex(s), P.sp.tr(s), P.sp.te(s,1),...
    'kap', nu.kap, 'dw', nu.b0, 'R2p', nu.r2p);
end

% acquire noiseless dess data
y.de = NaN([size(x.m0), S.de, E.de]);
for s = 1:S.de
  [y.de(:,:,s,1), y.de(:,:,s,2)] = dess_fun_v2(...
    x.m0, x.t1, x.t2,...
    P.de.aex(s), P.de.tr(s), P.de.te(s,1), P.de.te(s,2),...
    'kap', nu.kap, 'dw', nu.b0, 'R2p', nu.r2p);
end

% add complex gaussian noise
n.sp = c.sig * (randn(size(y.sp)) + 1i*randn(size(y.sp)));
n.de = c.sig * (randn(size(y.de)) + 1i*randn(size(y.de)));
ycc.sp = y.sp + n.sp;
ycc.de = y.de + n.de;

% use magnitude data as appropriate
if bool.mag.sp, ycc.sp = abs(ycc.sp); end
if bool.mag.de, ycc.de = abs(ycc.de); end

% set wmse dist.m0.supp to match map setting
tmp = cat(3,...
  reshape(ycc.sp(:,:,wght.sp(:,1)>0,:), [size(x.m0) sum(wght.sp>0)]),...
  reshape(ycc.de(:,:,wght.de(:,1)>0,:), [size(x.m0) sum(wght.de>0)]));
tmp = col(tmp);
dist.m0.supp = [eps div0(max(tmp),rff.snr)].';

% initial estimate
if header.tru
  x0 = x;
else
  x0.m0 = [];
  x0.t1 = [];
  x0.t2 = [];
  x0.inveff = [];
end

% perk hyperparameter optimization
try
  if header.br
    load('robust,broad.mat');
  else
    load('robust,tight.mat');
  end
catch
  % weighted nrmse measures
  nrmse.t12 = nan(length(hld.sig),1);
  nrmse.all = nan(length(hld.sig),1);
  nrmse.m0  = nan(length(hld.sig),1);
  nrmse.t1  = nan(length(hld.sig),1);
  nrmse.t2  = nan(length(hld.sig),1);
  
  % parameter weighting
  w.t12 = [0 1/2 1/2].';
  w.all = [1/3 1/3 1/3].';
  w.m0  = [1 0 0].';
  w.t1  = [0 1 0].';
  w.t2  = [0 0 1].';

  for h1=1:length(hld.sig)
    rff.std = hld.sig(h1);

    % parameter estimation
    fprintf('\n\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
    fprintf(...
      'Sensitivity: rff.std = 2^%.2f...\n',...
      log2(rff.std));
    fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
    opt.map = {...
      'nu', nu,...
      'x0', x0,...
      'wght.sp', wght.sp,...
      'wght.de', wght.de,...
      'dist.x.t1', dist.x.t1,...
      'dist.x.t2', dist.x.t2,...
      'meth.init', meth.init,...
      'meth.iter', meth.iter,...
      'rff', rff,...
      'inv.perk', inv.perk,...
      'stop.iter', stop.iter,...
      'bool', bool
    };
    [xhat, ~] = mri_m0t1t2inveff_map(ycc, P, opt.map{:});

    % holdout nrmse of perk
    tmp = fieldnames(nrmse);
    for i = 1:length(tmp)
      opt.wnrmse = {...
        'wght', w.(tmp{i})
      };
      nrmse.(tmp{i})(h1) = wnrmse(xhat.init, x, opt.wnrmse{:});
    end
  end
  
  % save wmse
  if header.sv
    if header.br
      save('robust,broad.mat', 'nrmse', 'w', 'hld');
    else
      save('robust,tight.mat', 'nrmse', 'w', 'hld');
    end
  end
end

% find minima
tmp = fieldnames(nrmse);
idx = cell(length(tmp),1);
for i = 1:length(tmp)
	[tmp2, idx{i}] = min(nrmse.(tmp{i}));
  fprintf('\n%3s-weighted nrmse minimized at sig = 2^%.2f with value %0.4f.\n',...
    tmp{i},...
    log2(hld.sig(idx{i})),...
    tmp2);
end

% weighted nrmse plots
if header.im
  % dynamic range
  if header.br
    dyn.t12 = col([0 2]);
  else
    dyn.t12 = col([0.1 0.3]);
  end
  
  % reduce broad plot range
  if header.br
    hld.sig = hld.sig(11:end);
    tmp = fieldnames(nrmse);
    for i = 1:length(tmp)
      nrmse.(tmp{i}) = nrmse.(tmp{i})(11:end);
      idx{i} = idx{i} - 10;
    end
  end
  
  % tick labels
  tick.x = log2(hld.sig(1:10:end));
  for i = 1:length(tick.x)
    if header.br
      ticklab.x{i} = sprintf('$2^{%d}$', tick.x(i));
    else
      ticklab.x{i} = sprintf('$2^{%0.1f}$', tick.x(i));
    end
  end
  
  % plots
  tmp = fieldnames(nrmse);
  for i = 1%1:length(tmp)
    figure;
    hold on;
    plot(log2(hld.sig), nrmse.(tmp{i}),...
      'linewidth', 2);
    plot(log2(hld.sig(idx{i})), nrmse.(tmp{i})(idx{i}), 'p',...
      'markeredgecolor', 'k',...
      'markerfacecolor', 'k',...
      'markersize', 16);
    %  normalization constant needed if bool.norm=1
    if bool.norm
      tmp2 = 0.16460719;
    else
      tmp2 = 1;
    end
    plot(log2(ones(2,1)*c.sig/tmp2), col(ylim), 'r--',...
      'linewidth', 2);
    hold off;
    xlabel('$\sigma$',...
      'interpreter', 'latex',...
      'fontsize', 24);
    ylabel('$\Psi(\sigma)$',...
      'interpreter', 'latex',...
      'fontsize', 24);
    set(gca,...
      'xtick', tick.x,...
      'ylim', dyn.(tmp{i}),...
      'ticklabelinterpreter', 'latex',...
      'xticklabel', ticklab.x,...
      'fontsize', 24);
    if header.pr && strcmp(tmp{i}, 't12')
      if header.br
        tmp2 = sprintf('robust,broad,w-%s.eps', tmp{i});
      else
        tmp2 = sprintf('robust,tight,w-%s.eps', tmp{i});
      end
      print('-depsc', tmp2);
    end
  end
end

