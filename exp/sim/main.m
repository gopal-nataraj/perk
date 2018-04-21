% script main.m
% comparison of perk versus local optimization and varpro in simulation
%   used simulated brainweb data
%   from http://mouldy.bic.mni.mcgill.ca/cgi/brainweb2?alias=phantom_1.0mm_msles2_crisp&download=1
%
% copyright 2017, gopal nataraj, university of michigan
%
% version control
%   2017-07-21      copied from mri_m0t1t2inveff_map_test.m
%   2017-08-04      modified to produce figures similar to ../hpd/main.m
%   2017-12-19      added local optimization

% irt
if (~exist('irtdir', 'var'))
  curdir = cd('~/Box/work/irt'); 
  irtdir = pwd;
  setup();
  cd(curdir);
end

% mapping
addpath('../../model/spgr');
addpath('../../model/dess');
addpath('../../model/ir');
addpath('../../model/se');
addpath('../../map/sense');
addpath('../../map/t1-t2');
addpath('../../misc');

% header options
header.prof = [1 0 0 0 0];              % specify profiles(s) to process
header.sv = 0;                          % save xhat
header.im = 1;                          % show xhat
header.pr = 0;                          % print images to file
header.stat = 0;                        % print sample statistics to file     

% constant declarations
c.b0 = 3.0;                             % field strength
c.sl = 81;                              % slice number
c.coil = 1;                             % number of rx coils to simulate
c.p = length(header.prof);        

c.sig.ir = 3.8607e-4;
c.sig.se = 3.8607e-4;
c.sig.sp = 3.8607e-4;
c.sig.de = 3.8607e-4;

c.rng.inveff.tru = 0.05;
c.rng.kap.tru = 0.2;
c.rng.kap.est = 0.2;
c.rng.b0.tru  = 0.00;                   % kHz
c.rng.b0.est  = 0.00;                   % kHz
c.rng.r2p.tru = 0.00;                   % kHz
c.rng.r2p.est = 0.00;                   % kHz    

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

% load digital phantom file
f.name = 'phantom_1.0mm_msles2_crisp.fld';
if exist('flip', 'builtin')
  f.label = flip(fld_read(f.name, 'slice', c.sl), 2);
else
  f.label = flipdim(fld_read(f.name, 'slice', c.sl), 2);
end
[c.n.x, c.n.y] = size(f.label);

% true latent parameter maps
[xx, yy] = ndgrid(linspace(-1,1,c.n.x), linspace(-1,1,c.n.y));

xtrue.m0 = NaN(c.n.x, c.n.y);
xtrue.t1 = NaN(c.n.x, c.n.y);
xtrue.t2 = NaN(c.n.x, c.n.y);
for i = 0:10
  f.param = mri_brainweb_params(i, 'b0', c.b0);
  xtrue.m0(f.label == i) = f.param.pd;
  xtrue.t1(f.label == i) = f.param.t1;
  xtrue.t2(f.label == i) = f.param.t2;
end
xtrue.inveff = (0.8+c.rng.inveff.tru/2) - c.rng.inveff.tru*(xx.^2 + yy.^2);

% true known parameter maps
nu.tru.kap = (1+c.rng.kap.tru/2) - c.rng.kap.tru*(xx.^2 + yy.^2);
nu.tru.b0  = (0+c.rng.b0.tru/2)  - c.rng.b0.tru *(xx.^2 + yy.^2);
nu.tru.r2p = (0+c.rng.r2p.tru/2) - c.rng.r2p.tru*(xx.^2 + yy.^2);

% estimated known parameter maps
nu.est.kap = (1+c.rng.kap.est/2) - c.rng.kap.est*(xx.^2 + yy.^2);
nu.est.b0  = (0+c.rng.b0.est/2)  - c.rng.b0.est *(xx.^2 + yy.^2);
nu.est.r2p = (0+c.rng.r2p.est/2) - c.rng.r2p.est*(xx.^2 + yy.^2);

% estimation options
meth.init = {'vpm','mom','perk'};
meth.iter = 'pgpm';
dist.x.t1.nsamp = 5; 
% dist.x.t1.nsamp = 500;
dist.x.t1.prior = 'logunif';
dist.x.t2.nsamp = 5; 
% dist.x.t2.nsamp = 500;
dist.x.t2.prior = 'logunif';
kmean.C = [20, 50, NaN];
rff.snr = 0.15;
rff.std = [];
rff.len = [];
rff.c = 2^0.6;
rff.H = 10^3;
rff.K = 10^5;
inv.perk = 2^-41;
boxcon.t1 = [100 3000'; 
stop.iter = [0, 0, 0]; 
% stop.iter = [0, 1000, 0];
bool.mag.ir = 0;
bool.mag.se = 0;
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

% coil sensitivities
if c.coil > 1
  smap.true = mri_sensemap_sim(...
    'nx', c.n.x,...
    'ny', c.n.y,...
    'dx', 1,...
    'dy', 1,...
    'ncoil', c.coil);
else
  smap.true = ones(c.n.x,c.n.y);
end

% acquire noiseless ir data
y.ir = NaN(c.n.x, c.n.y, S.ir, E.ir, c.coil);
for s = 1:S.ir
  tmp = IR_fun_v4(...
    xtrue.m0, xtrue.t1, xtrue.t2,...
    P.ir.tr(s), P.ir.ti(s), P.ir.te(s,1),...
    'inveff', xtrue.inveff, 'kap', nu.tru.kap, 'wf', nu.tru.b0,...
    'flip_inv', P.ir.ainv(s), 'flip_ex', P.ir.aex(s), 'flip_ref', P.ir.aref(s));
  for coil = 1:c.coil
    y.ir(:,:,s,1,coil) = smap.true(:,:,coil) .* tmp;
  end
end

% acquire noiseless se data
y.se = NaN(c.n.x, c.n.y, S.se, E.se, c.coil);
for s = 1:S.se
  tmp = SE_fun_v4(...
    xtrue.m0, xtrue.t1, xtrue.t2,...
    P.se.tr(s), P.se.te(s,1),...
    'kap', nu.tru.kap, 'wf', nu.tru.b0, 'flip_ex', P.se.aex(s), 'flip_ref', P.se.aref(s));
  for coil = 1:c.coil
    y.se(:,:,s,1,coil) = smap.true(:,:,coil) .* tmp;
  end
end

% acquire noiseless spgr data
y.sp = NaN(c.n.x, c.n.y, S.sp, E.sp, c.coil);
for s = 1:S.sp
  tmp = spgr_fun_v2(...
    xtrue.m0, xtrue.t1, xtrue.t2,...
    P.sp.aex(s), P.sp.tr(s), P.sp.te(s,1),...
    'kap', nu.tru.kap, 'dw', nu.tru.b0, 'R2p', nu.tru.r2p);
  for coil = 1:c.coil
    y.sp(:,:,s,1,coil) = smap.true(:,:,coil) .* tmp;
  end
end

% acquire noiseless dess data
y.de = NaN(c.n.x, c.n.y, S.de, E.de, c.coil);
for s = 1:S.de
  [tmp, tmp2] = dess_fun_v2(...
    xtrue.m0, xtrue.t1, xtrue.t2,...
    P.de.aex(s), P.de.tr(s), P.de.te(s,1), P.de.te(s,2),...
    'kap', nu.tru.kap, 'dw', nu.tru.b0, 'R2p', nu.tru.r2p);
  for coil = 1:c.coil
    y.de(:,:,s,1,coil) = smap.true(:,:,coil) .* tmp;
    y.de(:,:,s,2,coil) = smap.true(:,:,coil) .* tmp2;
  end
end

% random number generator seed
rng(0);

% add complex gaussian noise
n.ir = c.sig.ir * (randn(size(y.ir)) + 1i*randn(size(y.ir)));
n.se = c.sig.se * (randn(size(y.se)) + 1i*randn(size(y.se)));
n.sp = c.sig.sp * (randn(size(y.sp)) + 1i*randn(size(y.sp)));
n.de = c.sig.de * (randn(size(y.de)) + 1i*randn(size(y.de)));
y.ir = y.ir + n.ir;
y.se = y.se + n.se;
y.sp = y.sp + n.sp;
y.de = y.de + n.de;

% masks
mask.wm = f.label == 3;
mask.gm = f.label == 2;
mask.disp = f.label ~= 0;
mask.est = imdilate(mask.disp, strel('disk', 10));
mask.noise = ~imdilate(mask.disp, strel('disk', 20));

% compute snr
tmp = 'amp';
snr.ir.wm = snr_gn(y.ir, n.ir, mask.wm, 'unit', tmp);
snr.ir.gm = snr_gn(y.ir, n.ir, mask.gm, 'unit', tmp);
snr.se.wm = snr_gn(y.se, n.se, mask.wm, 'unit', tmp);
snr.se.gm = snr_gn(y.se, n.se, mask.gm, 'unit', tmp);
snr.sp.wm = snr_gn(y.sp, n.sp, mask.wm, 'unit', tmp);
snr.sp.gm = snr_gn(y.sp, n.sp, mask.gm, 'unit', tmp);
snr.de.wm = snr_gn(y.de, n.de, mask.wm, 'unit', tmp);
snr.de.gm = snr_gn(y.de, n.de, mask.gm, 'unit', tmp);

% coil-combine
ycc.ir = NaN(c.n.x, c.n.y, S.ir, E.ir);
ycc.se = NaN(c.n.x, c.n.y, S.se, E.se);
ycc.sp = NaN(c.n.x, c.n.y, S.sp, E.sp);
ycc.de = NaN(c.n.x, c.n.y, S.de, E.de);
if c.coil > 1
  opt.coil = {...
    'nouter', 10,...
    'thresh', 0.1,...
    'disp', 0};
  tmp = cat(3,...
    reshape(y.ir, [c.n.x c.n.y S.ir*E.ir c.coil]),...
    reshape(y.se, [c.n.x c.n.y S.se*E.se c.coil]),...
    reshape(y.sp, [c.n.x c.n.y S.sp*E.sp c.coil]),...
    reshape(y.de, [c.n.x c.n.y S.de*E.de c.coil]));
  [tmp, smap.est] = mri_multidata_coil_combine(tmp, opt.coil{:});
  
  yf = fieldnames(ycc);
  tmp2 = 0;
  for i = 1:length(yf)
    ycc.(yf{i}) = reshape(tmp(:,:,tmp2+1:tmp2+S.(yf{i})*E.(yf{i}),:),...
      [c.n.x c.n.y S.(yf{i}) E.(yf{i})]);
    tmp2 = tmp2 + S.(yf{i})*E.(yf{i});
  end
else
  ycc.ir = y.ir;
  ycc.se = y.se;
  ycc.sp = y.sp;
  ycc.de = y.de;
end

% use magnitude data as appropriate
if bool.mag.ir, ycc.ir = abs(ycc.ir); end
if bool.mag.se, ycc.se = abs(ycc.se); end
if bool.mag.sp, ycc.sp = abs(ycc.sp); end
if bool.mag.de, ycc.de = abs(ycc.de); end

% parameter estimation
try
  load(sprintf('im,xhat,sl-%u,nc-%u.mat', c.sl, c.coil));
catch
  xhat = struct('init', [], 'iter', []);
  t = struct('init', [], 'iter', []);
  for p = 1:c.p
    if header.prof(p) 
      fprintf('Testing profile %d: (%u,%u,%u,%u) ir/se/sp/de scans...\n', p,...
        sum(wght(p).ir(:,1)), sum(wght(p).se(:,1)),...
        sum(wght(p).sp(:,1)), sum(wght(p).de(:,1)));
      
      % set x0
      if p==5
        % for 4-se, set t1 from 4-ir ml est
        x0.t1 = x(4).t1.ml;
      else
        x0.t1 = [];
      end
      
      for i = 1:length(meth.init)
        % set t1,t2 dist support
        if strcmp(meth.init{i},'vpm')
          dist.x.t1.supp = [10^1.5 10^3.5].';
          if p==4
            dist.x.t2.supp = [10^6 10^6].';
          else
            dist.x.t2.supp = [10^0.5 10^3].';
          end
        elseif strcmp(meth.init{i},'perk')
          dist.x.t1.supp = [400 2000].';
          if p==4
            dist.x.t2.supp = [10^6 10^6].';
          else
            dist.x.t2.supp = [40 200].';
          end
        elseif strcmp(meth.init{i},'mom')
          dist.x.t1.supp = [];
          dist.x.t2.supp = [];
        else
          error('unknown method?');
        end
        
        % set options
        opt.map = {...
          'mask.disp', mask.disp,...
          'mask.est', mask.est,...
          'nu', nu.est,...
          'wght', wght(p),...
          'dist.x.t1', dist.x.t1,...
          'dist.x.t2', dist.x.t2,...
          'x0.t1', x0.t1,...
          'meth.init', meth.init{i},...
          'meth.iter', meth.iter,...
          'kmean.C', kmean.C(i),...
          'rff', rff,...
          'inv.perk', inv.perk,...
          'boxcon.t1', boxcon.t1,...
          'boxcon.t2', boxcon.t2,...
          'stop.iter', stop.iter(i),...
          'bool', bool...
        };
      
        % map
        [xhat(p,i), t(p,i)] = mri_m0t1t2inveff_map(ycc, P, opt.map{:});
      end
    else
      fprintf('Skipping profile %d: (%u,%u,%u,%u) ir/se/sp/de scans.\n', p,...
        sum(wght(p).ir(:,1)), sum(wght(p).se(:,1)),...
        sum(wght(p).sp(:,1)), sum(wght(p).de(:,1)));
    end
    
    % save xhat maps
    if header.sv
      save(sprintf('im,xhat,sl-%u,nc-%u.mat', c.sl, c.coil), 'xhat', 't');
    end
  end
end
  
% combine profiles 4 and 5
if header.prof(4) && header.prof(5)
  x(5,:).init.inveff = x(4,:).init.inveff;
  x(5,:).iter.inveff = x(4,:).iter.inveff;
  x(4,:) = [];
  c.p = 4;
end

% results display options
opt.prof = {'sp2de1', 'sp1de1', 'sp0de2', 'ir4se4'};
opt.type = {'im','err'};
opt.cmap = {'jet','gray'};
opt.xhat = {'m0','t1','t2'};
opt.meth = {'Truth', 'VPM', 'PGPM', 'PERK'};
opt.unit = {'a.u.','ms','ms'};
opt.roi = {'wm','gm'};

% colors
color.o = [1 0.6 0.2];                  % orange
color.f = [34 139 34]/256;              % forest green
color.m = [102 0 102]/256;              % maroon
color.l = [229 204 255]/256;            % lavendar
color.s = [255 0 127]/256;              % salmon
color.p = [255 153 255]/256;            % pink

% images
if header.im
  % dynamic ranges
  dyn.im.m0 = [0 1];
  dyn.im.t1 = [0 1500];
  dyn.im.t2 = [0 150];
  dyn.im.inveff = [0 1];
  
  dyn.err.m0 = [0 0.1];
  dyn.err.t1 = [0 150];
  dyn.err.t2 = [0 15];
  dyn.err.inveff = [0 0.1];
  
  for p = 1:c.p
    if header.prof(p)
      for m = 1:length(opt.cmap)
        for t = 1:length(opt.type)
          for i = 1:length(opt.xhat)
            tmp = cat(3,...
              xtrue.(opt.xhat{i}),...
              xhat(p,1).init.(opt.xhat{i}),...
              xhat(p,2).iter.(opt.xhat{i}),...
              xhat(p,3).init.(opt.xhat{i}));
            tmp = embed(masker(tmp, mask.wm|mask.gm), mask.wm|mask.gm);
            if strcmp(opt.type{t}, 'err')
              tmp = bsxfun(@minus, tmp, xtrue.(opt.xhat{i}));
              tmp = abs(tmp);
              tmp = embed(masker(tmp, mask.wm|mask.gm), mask.wm|mask.gm);
            end
            figure; 
            im('notick', 'row', 1, tmp, dyn.(opt.type{t}).(opt.xhat{i}), 'cbar' ,' ');
            colormap(gca, opt.cmap{m});...
            tmp = colorbar;
            ylabel(tmp, opt.unit{i});
            set(tmp, 'ytick',...
              linspace(dyn.(opt.type{t}).(opt.xhat{i})(1), dyn.(opt.type{t}).(opt.xhat{i})(2), 6));
            hold on;
            tmp = length(opt.meth);
            if strcmp(opt.type{t}, 'im')
              text(col(c.n.x/2-1:c.n.x:(2*tmp-1)*c.n.x/2), zeros(tmp,1), col(opt.meth),...
                'VerticalAlignment', 'bottom',...
                'HorizontalAlignment', 'center',...
                'FontSize', 16,...
                'Color', 'k');
            end
            if strcmp(opt.type{t}, 'im')
              tmp = upper(opt.xhat{i});
            else
              tmp = [upper(opt.xhat{i}), ' Magnitude Error'];
              text(2, 2, '(x10 magnified)',...
                'VerticalAlignment', 'top',...
                'HorizontalAlignment', 'left',...
                'FontSize', 12,...
                'Color', 'w');
            end
            text(0, c.n.y/2-1, tmp,...
              'VerticalAlignment', 'bottom',...
              'HorizontalAlignment', 'center',...
              'FontSize', 14,...
              'Color', 'k',...
              'Rotation', 90);
            hold off;
            if header.pr
              tmp = sprintf('%s,sl-%u,%s,%s,%s',...
                opt.prof{p}, c.sl, opt.xhat{i}, opt.type{t}, opt.cmap{m});
              print('-depsc', tmp);
            end
          end
        end
      end
    end
  end
end

% summary statistics
for p = 1:c.p
  if header.prof(p)
    if header.stat
      tmp = sprintf('%s,sl-%u,stat', opt.prof{p}, c.sl);
      fid = fopen(tmp,'w');
      fprintf(fid, 'estimator statistics for profile %s, slice %u\n', opt.prof{p}, c.sl);
      fprintf(fid, '\trows distinguish regions of interest\n');
      fprintf(fid, '\tcolumns distinguish estimation methods\n');
    end
    for i = 1:length(opt.xhat)
      if header.stat
        fprintf(fid, '\n%s\n', opt.xhat{i});
        tmp = '';
        for z = 1:length(opt.meth)
          tmp = strcat(tmp, '%28s');
        end
        tmp = strcat(tmp, '\n');
        fprintf(fid, tmp, opt.meth{:});
      end
      for r = 1:length(opt.roi)
        for m = 1:length(opt.meth)
          switch opt.meth{m}
            case 'Truth'
              tmp = xtrue;
            case 'VPM'
              tmp = xhat(p,1).init;
            case 'MOM'
              tmp = xhat(p,2).init;
            case 'PGPM'
              tmp = xhat(p,2).iter;
            case 'PERK'
              tmp = xhat(p,3).init;
          end
          summ.(opt.prof{p}).(opt.xhat{i})(m,r) = stat(...
            masker(tmp.(opt.xhat{i}), mask.(opt.roi{r})),...
            'true', mean(masker(xtrue.(opt.xhat{i}), mask.(opt.roi{r}))));
          if header.stat
            fprintf(fid, '\t%12.4f\t%c%8.4f',...
              summ.(opt.prof{p}).(opt.xhat{i})(m,r).mean,...
              char(177),...
              summ.(opt.prof{p}).(opt.xhat{i})(m,r).std);
          end
        end
        if header.stat 
          fprintf(fid, '\n');
        end
      end
    end
    if header.stat
      fclose(fid);
    end
  end
end
