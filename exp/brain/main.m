% script main.m
% comparison of perk versus local optimization and varpro in brain
%   brain data acquired on 2016-05-31
%
% copyright 2017, gopal nataraj, university of michigan
%
% version control
%   2017-07-19      original
%   2017-08-11      coil-combine data here, with stronger regularization
%   2018-01-23      added local optimization

% irt
if ~exist('irtdir', 'var')
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
header.sc   = 1;                        % scale kap using prior calibration
header.reg  = 1;                        % register coil-combined images
header.prof = [1 0 0 0 0];              % specify profile(s) to process
header.sv   = 0;                        % save xhat
header.rs   = 0;                        % roi manual selection
header.im   = 1;                        % show xhat
header.roi  = 1;                        % show roi labels in images
header.pr   = 0;                        % print images to file
header.stat = 0;                        % print sample statistics to file

% constant declarations
n.x = 256;
n.y = 256;
n.z = 8;
n.c = 32;
n.sl = 5;
n.p = length(header.prof);

fov.x = 240;                            % mm
fov.y = 240;                            % mm
fov.z = 5;                              % mm

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
for p = 1:n.p
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

% coil combination options
opt.coil = {...
  'log2b', -2,...
  'thresh', 0.1,...
  'fwhm', 0,...
  'nouter', 10,...
  'disp', 1 ...
};

opt.b1 = {...
  'coilOpt', opt.coil,...
  'scale', 1,...
  'reg.beta', 2^2,...
  'stop.iter', 5000,...
  'bool.chat', 1,...
  'bool.disp', 1 ...
};

% load b1 maps
load(sprintf('im,nu,sl-%u.mat', n.sl));

% uniform b1/kap scaling
if header.sc
  scale = 1.1;                          % from external calibration
  nu.b1.mom = nu.b1.mom * scale;
  nu.b1.rls = nu.b1.rls * scale;
  nu.kap.mom = nu.kap.mom * scale;
  nu.kap.rls = nu.kap.rls * scale;
end

% load coil-combined image data
load('im,irse-se,coil-comb.mat');
tmp = ycc;
load(sprintf('im,spgr-dess,coil-comb,sl-%u.mat', n.sl));
ycc = catstruct(tmp, ycc);

% register coil-combined images and nu.kap.rls to first ycc.ir image
if header.reg && ~exist(sprintf('im,xhat,sl-%u.mat', n.sl), 'file')
  fprintf('\nCoil-combined image registration to reg.fixed:\n');
  reg.fixed = abs(ycc.ir(:,:,1));
  reg.ttype = 'rigid';
  [reg.optim, reg.metric] = imregconfig('multimodal');
  
  tmp = fieldnames(ycc);
  for i = 1:length(tmp)
    for s = 1:S.(tmp{i})
      fprintf('  registering %s dataset %u...', upper(tmp{i}), s);
      tmp2 = imregtform(abs(ycc.(tmp{i})(:,:,s,1)),...
        reg.fixed, reg.ttype, reg.optim, reg.metric);
      for j = 1:E.(tmp{i})
        ycc.(tmp{i})(:,:,s,j) = imwarp(ycc.(tmp{i})(:,:,s,j),...
          tmp2, 'OutputView', imref2d(size(reg.fixed)));
      end
      fprintf('done.\n');
    end
  end
  
  fprintf('  registering flip angle scaling map...');
  tmp = imregtform(abs(nu.kap.rls), reg.fixed, reg.ttype, reg.optim, reg.metric);
  nu.kap.rls = imwarp(nu.kap.rls, tmp, 'OutputView', imref2d(size(reg.fixed)));
  fprintf('done.\n\n');
end

% create masks from short te image
tmp = squeeze(abs(ycc.se(:,:,1)));
mask.thresh = 0.05;
mask.t = imfill(tmp > mask.thresh * max(col(tmp)), 'holes');
mask.t = imdilate(~imdilate(~mask.t, strel('disk', 5)), strel('disk', 5));
mask.b = imdilate(mask.t, strel('disk', 10));

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
boxcon.t1 = [100 3000].';
boxcon.t2 = [10 700].';
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

% use magnitude data as appropriate
tmp = fieldnames(bool.mag);
for i = 1:length(tmp)
  if bool.mag.(tmp{i})
    ycc.(tmp{i}) = abs(ycc.(tmp{i}));
  end
end

% parameter estimation
try
  load(sprintf('im,xhat,sl-%u.mat', n.sl));
catch
  xhat = struct('init', [], 'iter', []);
  t = struct('init', [], 'iter', []);
  for p = 1:n.p
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
          'mask.disp', mask.t,...
          'mask.est', mask.b,...
          'nu.kap', double(nu.kap.rls),...
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
          'bool',  bool...
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
      save(sprintf('im,xhat,sl-%u.mat', n.sl), 'xhat', 't');
    end
  end
end

% roi selection
if header.rs
  tmp = 'ROI selection: [m]anual or [s]aved? ';
  roi.sel = input(tmp, 's');
else
  roi.sel = 's';
end
roi.label = {'wm-ra'; 'wm-la'; 'wm-rp'; 'wm-lp'; 'gm-a'};
n.roi = length(roi.label);

switch roi.sel
  case 'm'
    roi.mask = false(n.x, n.y, n.roi);
    for r = 1:n.roi
      while true
        fprintf('Select %s polygon.\n', roi.label{r});
        tmp = roipoly(x(n.p).t1.rls/2000);
        roi.mask(:,:,r) = roi.mask(:,:,r) | tmp;
        tmp = input('Add another polygon to this ROI [*/n]? ', 's');
        if strcmp(tmp,'n'), break; end
      end
    end

    % save mask
    dir = cd('../../data/Brain_05,31,16');
    if exist('roi-mask.mat', 'file')
      tmp = input('Overwrite previously-defined ROI mask [y/n]? ', 's');
      switch tmp
        case 'y'
          save('roi-mask.mat', 'roi'); 
          fprintf('roi-mask.mat overwritten.\n');
        case 'n'
          fprintf('New ROI mask not saved.\n');
        otherwise
          error('Unrecognized input.');
      end
    else
      save('roi-mask.mat', 'roi');
    end
    cd(dir);

  case 's'
    try
      load('roi-mask.mat');
    catch
      error('ROI masks not found!');
    end

  otherwise
    error('Unrecognized input!');
end

% roi boundaries
roi.bound = cell(n.roi,1);
for r = 1:n.roi
  roi.bound{r} = bwboundaries(roi.mask(:,:,r), 'noholes');
end

% combine profiles 4 and 5
if header.prof(4) && header.prof(5)
  x(5,:).init.inveff = x(4,:).init.inveff;
  x(5,:).iter.inveff = x(4,:).iter.inveff;
  x(4,:) = [];
  n.p = 4;
end

% results display options
opt.prof = {'sp2de1', 'sp1de1', 'sp0de2', 'ir4se4'};
opt.cmap = {'jet','gray'};
opt.xhat = {'m0','t1','t2'};
opt.meth = {'VPM', 'PGPM', 'PERK'};
opt.unit = {'a.u.','ms','ms'};

% colors
roi.color = {'y';'m';'g';'b';'c'};

% images
if header.im
  dyn.m0 = [0 30];
  dyn.t1 = [600 1600];
  dyn.t2 = [20 120];
  dyn.inveff = [0.5 1];
  
  crop.x = 35:222;
  crop.y = 11:256;
  n.xcrop = length(crop.x);
  n.ycrop = length(crop.y);

  for p = 1:n.p
    if header.prof(p)
      for m = 1:length(opt.cmap)
        for i = 1:length(opt.xhat)
          tmp = cat(3,...
            xhat(p,1).init.(opt.xhat{i})(crop.x,crop.y),...
            xhat(p,2).iter.(opt.xhat{i})(crop.x,crop.y),...
            xhat(p,3).init.(opt.xhat{i})(crop.x,crop.y));
          figure; 
          im('notick', 'row', 1, tmp, dyn.(opt.xhat{i}), 'cbar' ,' ');
          colormap(gca, opt.cmap{m});...
          tmp = colorbar;
          ylabel(tmp, opt.unit{i});
          set(tmp, 'ytick', linspace(dyn.(opt.xhat{i})(1), dyn.(opt.xhat{i})(2), 6));
          hold on;
          tmp = length(opt.meth);
          if strcmp(opt.xhat{i}, 'm0')
            text(col(n.xcrop/2-1:n.xcrop:(2*tmp-1)*n.xcrop/2), zeros(tmp,1), col(opt.meth),...
              'VerticalAlignment', 'bottom',...
              'HorizontalAlignment', 'center',...
              'FontSize', 16,...
              'Color', 'k');
          end
          text(0, n.ycrop/2-1, upper(opt.xhat{i}),...
            'VerticalAlignment', 'bottom',...
            'HorizontalAlignment', 'center',...
            'FontSize', 16,...
            'Color', 'k',...
            'Rotation', 90);
          if header.roi && strcmp(opt.xhat{i},'m0') && strcmp(opt.cmap{m},'gray')
            for r = 1:n.roi
              tmp = roi.bound{r};
              for b = 1:length(roi.bound{r})
                tmp2 = tmp{b};
                plot(tmp2(:,1)-crop.x(1)+1, tmp2(:,2)-crop.y(1)+1, roi.color{r},...
                  'LineWidth', 0.5);
              end
            end
            text(n.xcrop/2, 0, 'A',...
              'VerticalAlignment', 'top',...
              'HorizontalAlignment', 'center',...
              'FontSize', 16,...
              'Color', 'w');
            text(n.xcrop, n.ycrop/2, 'L',...
              'VerticalAlignment', 'middle',...
              'HorizontalAlignment', 'right',...
              'FontSize', 16,...
              'Color', 'w');
            text(n.xcrop/2, n.ycrop, 'P',...
              'VerticalAlignment', 'bottom',...
              'HorizontalAlignment', 'center',...
              'FontSize', 16,...
              'Color', 'w');
            text(1, n.ycrop/2, 'R',...
              'VerticalAlignment', 'middle',...
              'HorizontalAlignment', 'left',...
              'FontSize', 16,...
              'Color', 'w');
          end
          hold off;
          if header.pr
            tmp = sprintf('%s,sl-%u,%s,im-%s',...
              opt.prof{p}, n.sl, opt.xhat{i}, opt.cmap{m});
            print('-depsc', tmp);
          end
        end
      end
    end
  end
end
     
% summary statistics
for p = 1:n.p
  if header.prof(p)
    if header.stat
      tmp = sprintf('%s,sl-%u,stat', opt.prof{p}, n.sl);
      fid = fopen(tmp,'w');
      fprintf(fid, 'estimator statistics for profile %s, slice %u\n', opt.prof{p}, n.sl);
      fprintf(fid, '\trows distinguish regions of interest\n');
      fprintf(fid, '\tcolumns distinguish estimation methods\n');
    end
    for i = 1:length(opt.xhat)
      if header.stat && ~strcmp(opt.xhat{i},'m0')
        fprintf(fid, '\n%s\n', opt.xhat{i});
        tmp = '';
        for z = 1:length(opt.meth)
          tmp = strcat(tmp, '%28s');
        end
        tmp = strcat(tmp, '\n');
        fprintf(fid, tmp, opt.meth{:});
      end
      for r = 1:n.roi
        for m = 1:length(opt.meth)
          switch opt.meth{m}
            case 'VPM'
              tmp = xhat(p,1).init;
            case 'MOM'
              tmp = xhat(p,2).init;
            case 'PGPM'
              tmp = xhat(p,2).iter;
            case 'PERK'
              tmp = xhat(p,3).init;
          end
          summ.(opt.prof{p}).(opt.xhat{i})(m,r) = ...
            stat(masker(tmp.(opt.xhat{i}), stackpick(roi.mask,r)));
          summlog10.(opt.prof{p}).(opt.xhat{i})(m,r) = ...
            stat(masker(log10(tmp.(opt.xhat{i})), stackpick(roi.mask,r)));
          if header.stat && ~strcmp(opt.xhat{i},'m0')
            fprintf(fid, '\t%12.2f\t%c%8.2f',...
              summ.(opt.prof{p}).(opt.xhat{i})(m,r).mean,...
              char(177),...
              summ.(opt.prof{p}).(opt.xhat{i})(m,r).std);
          end
        end
        if header.stat && ~strcmp(opt.xhat{i},'m0')
          fprintf(fid, '\n');
        end
      end
    end
    if header.stat
      fclose(fid);
    end
  end
end
