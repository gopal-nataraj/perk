% conditional bias and covariance plots
%
% copyright 2017, gopal nataraj, university of michigan
%
% version control
%   2017-09-04      original

% irt
if (~exist('irtdir', 'var'))
  curdir = cd('~/Box/work/irt');
  irtdir = pwd;
  setup();
  cd(curdir);
end

% signal models
addpath('../model/spgr');
addpath('../model/dess');
addpath('../model/ir');
addpath('../model/se');

% header options
header.prof = [1 0 0 0 0];              % specify profiles(s) to process
header.pr = 0;                          % print figures to file

% constant declarations
c.p = length(header.prof);        
c.sig = 3.8607e-4;

% acquisition parameters
P.ir.ti   = [50 150 450 1350]';         % ms
S.ir      = length(P.ir.ti);
E.ir      = 1;
P.ir.tr   = ones(S.ir,1) * 1400;        % ms
P.ir.te   = ones(S.ir,E.ir) * 14;       % ms
P.ir.ainv = ones(S.ir,1) * pi;          % rad
P.ir.aex  = ones(S.ir,1) * pi/2;        % rad
P.ir.aref = ones(S.ir,1) * pi;          % rad

P.se.te   = [10 30 60 150]';            % ms
S.se      = length(P.se.te);
E.se      = 1;         
P.se.tr   = ones(S.se,1) * 1000;        % ms
P.se.aex  = ones(S.se,1) * pi/2;        % rad
P.se.aref = ones(S.se,1) * pi;          % rad

P.sp.aex  = [15 5 15]' * pi/180;        % rad
P.sp.tr   = [12.2 12.2 13.9]';          % ms
S.sp      = length(P.sp.aex);          
E.sp      = 1;                          
P.sp.te   = ones(S.sp,E.sp) * 4.67;     % ms

P.de.aex  = [35 30 10 10]' * pi/180;    % rad
P.de.tr   = [24.4 17.5 28 17.5]';       % ms
S.de      = length(P.de.aex);           
E.de      = 2;
P.de.te   = ones(S.de,E.de) * 4.67;     % ms

% set data weights based on profile
wght = struct([]);
for p = 1:c.p
  switch p
    case 1                              % (2 sp, 1 de)
      wght(p).ir = [0 0 0 0]'   * ones(1,E.ir);                   
      wght(p).se = [0 0 0 0]'   * ones(1,E.se);   
      wght(p).sp = [1 1 0]'     * ones(1,E.sp);
      wght(p).de = [0 1 0 0]'   * ones(1,E.de);
    case 2                              % (1 sp, 1 de)
      wght(p).ir = [0 0 0 0]'   * ones(1,E.ir);                   
      wght(p).se = [0 0 0 0]'   * ones(1,E.se);   
      wght(p).sp = [0 0 1]'     * ones(1,E.sp);
      wght(p).de = [0 0 1 0]'   * ones(1,E.de);
    case 3                              % (0 sp, 2 de)
      wght(p).ir = [0 0 0 0]'   * ones(1,E.ir);                   
      wght(p).se = [0 0 0 0]'   * ones(1,E.se);   
      wght(p).sp = [0 0 0]'     * ones(1,E.sp);
      wght(p).de = [1 0 0 1]'   * ones(1,E.de);
    case 4                              % (4 ir, 0 se)
      wght(p).ir = [1 1 1 1]'   * ones(1,E.ir);
      wght(p).se = [0 0 0 0]'   * ones(1,E.se);
      wght(p).sp = [0 0 0]'     * ones(1,E.sp);
      wght(p).de = [0 0 0 0]'   * ones(1,E.de);
    case 5                              % (0 ir, 4 se)
      wght(p).ir = [0 0 0 0]'   * ones(1,E.ir);
      wght(p).se = [1 1 1 1]'   * ones(1,E.se);   
      wght(p).sp = [0 0 0]'     * ones(1,E.sp);
      wght(p).de = [0 0 0 0]'   * ones(1,E.de);
  end
end

% define grid on which to evaluate
eval.x.m0     = col(linspace(1,1,1));
eval.x.t1     = col(logspace(log10(400), log10(2000), 50));
eval.x.t2     = col(logspace(log10(40), log10(200), 50));
eval.x.inveff = col(linspace(1,1,1));
eval.nu.kap   = col(linspace(1,1,1));
eval.nu.b0    = col(linspace(0,0,1));
eval.nu.r2p   = col(linspace(0,0,1));
[x.m0, x.t1, x.t2, x.inveff, nu.kap, nu.b0, nu.r2p] = ndgrid(...
  eval.x.m0, eval.x.t1, eval.x.t2, eval.x.inveff,...
  eval.nu.kap, eval.nu.b0, eval.nu.r2p...
);

% define masks
[~,idx.m0] = min(abs(eval.x.m0-1));
[~,idx.t1] = min(abs(eval.x.t1-1000));
[~,idx.t2] = min(abs(eval.x.t2-80));
mask.est = false(size(x.t1));
mask.est(:,idx.t1,idx.t2) = true;
mask.est(idx.m0,:,idx.t2) = true;
mask.est(idx.m0,idx.t1,:) = true;
mask.disp = mask.est;

% bias-covariance options
dist.m0.supp = [0.99 1.01].';
dist.m0.prior = 'unif';
dist.t1.supp = [400 2000].';
dist.t1.prior = 'logunif';
dist.t2.supp = [40 200].';
dist.t2.prior = 'logunif';
dist.inveff.supp = [0.51 1].';
dist.inveff.prior = 'unif';
rff.snr = 0.15;
rff.std = c.sig;
rff.len.y = [];
rff.len.nu = [];
rff.c = 2^0.6;
rff.H = 10^2; 
% rff.H = 10^3;
rff.K = 10^4; 
% rff.K = 10^5;
inv.krr = 2^-30; 
% inv.krr = 2^-41;
bool.mag.ir = 0;
bool.mag.se = 0;
bool.mag.sp = 1;
bool.mag.de = 1;
bool.chat = 1;
bool.norm = 1;
bool.reset = 1;
bool.rfftst = 0;
bool.disp = 0;

% bias-covariance calculation
stat = struct('bias', [], 'cov', []);
for p = 1:c.p
  if header.prof(p)
    fprintf('Testing profile %d: (%u,%u,%u,%u) ir/se/sp/de scans...\n', p,...
      sum(wght(p).ir(:,1)), sum(wght(p).se(:,1)),...
      sum(wght(p).sp(:,1)), sum(wght(p).de(:,1)));
    
    % set rff.len.y from simulation
    if p==1
      rff.len.y = [...
        0.145193780995431;...
        0.158289232720290;...
        0.281746574027058;...
        0.208023018335361...
      ];
    else
      warn('rff.len.y not set!');
    end
    
    opt.map = {...
      'mask', mask,...
      'nu', nu,...
      'dist', dist,...
      'wght', wght(p),...
      'rff', rff,...
      'inv', inv,...
      'bool', bool...
    };
    [stat(p)] = bias_cov(x, P, opt.map{:});
  end
end

% figures
opt.xhat = {'m0', 't1', 't2'};
opt.unit = {'a.u.', 'ms', 'ms'};
for p = 1:c.p
  if header.prof(p)
    figure;
    n = 0;
    for i = 2:length(opt.xhat)
      switch i
        case 1
          tmp = stat(p).bias(:,idx.t1,idx.t2,:);
        case 2
          tmp = stat(p).bias(idx.m0,:,idx.t2,:);
        case 3
          tmp = stat(p).bias(idx.m0,idx.t1,:,:);
      end
      tmp = squeeze(tmp);
      
      for j = 2:length(opt.xhat)
        n = n+1;
        subplot(2,2,n);
        plot(eval.x.(opt.xhat{i}), col(tmp(:,j)));
        set(gca, 'xlim', minmax(eval.x.(opt.xhat{i})));
        xlabel(sprintf('%s (%s)', opt.xhat{i}, opt.unit{i}));
        ylabel(sprintf('%s estimation bias (%s)', opt.xhat{j}, opt.unit{j}));
      end
    end
    
    figure; 
    n = 0;
    for i = 2:length(opt.xhat)
      switch i
        case 1
          tmp = stat(p).cov(:,idx.t1,idx.t2,:,:);
        case 2
          tmp = stat(p).cov(idx.m0,:,idx.t2,:,:);
        case 3
          tmp = stat(p).cov(idx.m0,idx.t1,:,:,:);
      end
      tmp = sqrt(squeeze(tmp));
      
      for j = 2:length(opt.xhat)
        n = n+1;
        subplot(2,2,n);
        plot(eval.x.(opt.xhat{i}), col(tmp(:,j,j)));
        set(gca, 'xlim', minmax(eval.x.(opt.xhat{i})));
        xlabel(sprintf('%s (%s)', opt.xhat{i}, opt.unit{i}));
        ylabel(sprintf('%s estimation std dev (%s)', opt.xhat{j}, opt.unit{j}));
      end
    end   
  end
end
