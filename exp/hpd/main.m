% script main.m
% comparison of perk versus local optimization and varpro in phantom
%   High Precision Devices (HPD) phantom data acquired on 2016-06-20
%
% copyright 2017, gopal nataraj, university of michigan
%
% version control
%   2017-07-17      original
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
header.prof = [1 0 0 0 0];              % specify profile(s) to process
header.sv   = 0;                        % save xhat
header.im   = 1;                        % show xhat
header.roi  = 1;                        % show roi labels in images
header.pr   = 0;                        % print images to file
header.stat = 0;                        % print sample statistics to file
header.cmp  = 0;                        % print t1/t2 comparison plots

% constant declarations
n.x = 256;
n.y = 256;
n.z = 8;
n.c = 8;
n.sl = 6;
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

% load images
load('im,irse-se,coil-comb.mat');
tmp = ycc;
load(sprintf('im,spgr-dess,coil-comb,sl-%u.mat', n.sl));
ycc = catstruct(tmp, ycc);

% create masks from short te image
tmp = squeeze(abs(ycc.se(:,:,1)));
mask.thresh = 0.05;
mask.t = imfill(tmp > mask.thresh * max(col(tmp)), 'holes');
mask.t = imdilate(~imdilate(~mask.t, strel('disk', 5)), strel('disk', 5));
mask.b = imdilate(mask.t, strel('disk', 10));

% estimation options
meth.init = {'vpm','mom','krr'};
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
inv.krr = 2^-41;
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
  load(sprintf('dist-broad/im,xhat,sl-%u.mat', n.sl));
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
        elseif strcmp(meth.init{i},'krr')
          dist.x.t1.supp = [400 2000].';
%           dist.x.t1.supp = [10^1.5 10^3.5].';
          if p==4
            dist.x.t2.supp = [10^6 10^6].';
          else
            dist.x.t2.supp = [40 200].';
%             dist.x.t2.supp = [10^0.5 10^3].';
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
          'inv.krr', inv.krr,...
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

% roi locations
roi.x = [129 160 180 179 160 128  97  77  77  98 107 150 150 107].';
roi.y = [ 75  85 112 145 171 181 171 143 111  84 106 107 149 149].';
n.roi = length(roi.x);
roi.r = ones(n.roi,1) * 3;

% roi masks
mask.roi = false(n.x,n.y,n.roi);
[tmp,tmp2] = ndgrid(1:n.x, 1:n.y);
for r = 1:n.roi
  mask.roi(:,:,r) = (tmp-roi.x(r)).^2 + (tmp2-roi.y(r)).^2 <= roi.r(r)^2;
  roi.label{r} = sprintf('Vial%2u',r);
end
mask.all = sum(mask.roi,3)>0;

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
color.o = [1 0.6 0.2];                  % orange
color.f = [34 139 34]/256;              % forest green
color.m = [102 0 102]/256;              % maroon
color.l = [229 204 255]/256;            % lavendar
color.s = [255 0 127]/256;              % salmon
color.p = [255 153 255]/256;            % pink

% tight and loose rois
roi.vial.t.t1 = 5:7;
roi.vial.t.t2 = 6:7;
roi.vial.b.t1 = 3:9;
roi.vial.b.t2 = 4:8;

% images
if header.im
  dyn.m0 = [0 30];
  dyn.t1 = [0 3000];
  dyn.t2 = [0 500];
  dyn.inveff = [0.5 1];
  for p = 1:n.p
    if header.prof(p)
      for m = 1:length(opt.cmap)
        for i = 1:length(opt.xhat)
          tmp = cat(3,...
            xhat(p,1).init.(opt.xhat{i}),...
            xhat(p,2).iter.(opt.xhat{i}),...
            xhat(p,3).init.(opt.xhat{i}));
          figure; 
          im('notick', 'row', 1, tmp, dyn.(opt.xhat{i}), 'cbar' ,' ');
          colormap(gca, opt.cmap{m});...
          tmp = colorbar;
          ylabel(tmp, opt.unit{i});
          set(tmp, 'ytick', linspace(dyn.(opt.xhat{i})(1), dyn.(opt.xhat{i})(2), 6));
          hold on;
          tmp = length(opt.meth);
          if strcmp(opt.xhat{i}, 'm0')
            text(col(n.x/2-1:n.x:(2*tmp-1)*n.x/2), zeros(tmp,1), col(opt.meth),...
              'VerticalAlignment', 'bottom',...
              'HorizontalAlignment', 'center',...
              'FontSize', 16,...
              'Color', 'k');
          end
          text(0, n.y/2-1, upper(opt.xhat{i}),...
            'VerticalAlignment', 'bottom',...
            'HorizontalAlignment', 'center',...
            'FontSize', 16,...
            'Color', 'k',...
            'Rotation', 90);
          if header.roi && strcmp(opt.xhat{i},'m0') && strcmp(opt.cmap{m},'gray')
            for r = 1:n.roi
              if min(roi.vial.b.t2)<=r && r<=max(roi.vial.b.t2)
                tmp = 'y';
              else
                tmp = 'k';
              end
              text(roi.x(r), roi.y(r), num2str(r),...
                'VerticalAlignment', 'middle',...
                'HorizontalAlignment', 'center',...
                'FontSize', 6,...
                'Color', tmp,... 
                'FontWeight', 'bold');
            end
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

% nist statistics
nist.t1.mean = ...
  [2480   2173   1907   1604    1332    1044    801.7...
    608.6  458.4  336.5  244.2   176.6   126.9   90.9]';
nist.t1.std = ...
  [  10.8   14.7   10.3    7.2     0.8     3.2    1.70...
      1.03   0.33   0.18   0.09    0.09   0.03    0.05]'; 
nist.t2.mean = ...
  [ 581.3  403.5  278.1  190.94  133.27   96.89  64.07...
     46.42  31.97  22.56  15.813  11.237   7.911  5.592]';
nist.t2.std = ...
  [   0.39   0.55   0.28   0.011   0.073   0.049  0.034...
      0.014  0.083  0.012  0.0061  0.0057  0.0037 0.0055]';
          
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
        tmp = strcat(tmp, '%28s\n');
        fprintf(fid, tmp, opt.meth{:}, 'NIST');
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
            stat(masker(tmp.(opt.xhat{i}), stackpick(mask.roi,r)));
          summlog10.(opt.prof{p}).(opt.xhat{i})(m,r) = ...
            stat(masker(log10(tmp.(opt.xhat{i})), stackpick(mask.roi,r)));
          if header.stat && ~strcmp(opt.xhat{i},'m0')
            fprintf(fid, '\t%12.2f\t%c%8.2f',...
              summ.(opt.prof{p}).(opt.xhat{i})(m,r).mean,...
              char(177),...
              summ.(opt.prof{p}).(opt.xhat{i})(m,r).std);
          end
        end
        if header.stat && ~strcmp(opt.xhat{i},'m0')
          fprintf(fid, '\t%12.2f\t%c%8.2f\n',...
            nist.(opt.xhat{i}).mean(r),...
            char(177),...
            nist.(opt.xhat{i}).std(r));
        end
      end
    end
    if header.stat
      fclose(fid);
    end
  end
end

% t1/t2 roi fill areas
roi.dyn.t.t1  = [800 1400];
roi.dyn.t.t2  = [50 120];
roi.dyn.b.t1  = [400 2000];
roi.dyn.b.t2  = [40 200];

% fig options
fig.c = {'b', 'r', color.f, color.s, color.m};
fig.m = {'o', 'v', 'd', 's', '^'};
fig.f = {'c', 'm', 'g', color.p, color.l};
fig.tick.t1 = log10([80:10:100 200:100:1000 2000:1000:3000]);
fig.tick.t2 = log10([5:1:10 20:10:100 200:100:600]);

% comparison plots
for p = 1:n.p
  if header.prof(p)
    for i = 1:length(opt.xhat)
      if ~strcmp(opt.xhat{i},'m0')
        % roi boxes
        roi.box.t.(opt.xhat{i}) = {...
          log10(kron(roi.dyn.t.(opt.xhat{i}), [1 1])),...
          log10([roi.dyn.t.(opt.xhat{i}) fliplr(roi.dyn.t.(opt.xhat{i}))])};
        roi.box.b.(opt.xhat{i}) = {...
          log10(kron(roi.dyn.b.(opt.xhat{i}), [1 1])),...
          log10([roi.dyn.b.(opt.xhat{i}) fliplr(roi.dyn.b.(opt.xhat{i}))])};

        % tick labels
        fig.ticklab.(opt.xhat{i}) = cell(length(fig.tick.(opt.xhat{i})),1);
        for j = 1:length(fig.tick.(opt.xhat{i}))
          if rem(fig.tick.(opt.xhat{i})(j),1)==0
            fig.ticklab.(opt.xhat{i}){j} = sprintf('%u', 10^fig.tick.(opt.xhat{i})(j));
          else
            fig.ticklab.(opt.xhat{i}){j} = ' ';
          end
        end

        % figure
        figure;
        hold on;
        fill(roi.box.b.(opt.xhat{i}){:}, 'y');
        handles = cell(length(opt.meth)+1,1);
        tmp = col(minmax(fig.tick.(opt.xhat{i})));
        tmp = log10(logspace(tmp(1), tmp(2), 1000));
        handles{1} = plot(tmp, tmp, 'k--', 'LineWidth', 1);
        tmp = log10(nist.(opt.xhat{i}).mean);
        tmp2 = nist.(opt.xhat{i}).mean ./ nist.(opt.xhat{i}).std;
        tmp2 = tmp ./ tmp2;
        for m = 1:length(opt.meth)
          handles{m+1} = errorbarxy(...
            tmp, col([summlog10.(opt.prof{p}).(opt.xhat{i})(m,:).mean]),...
            tmp2, col([summlog10.(opt.prof{p}).(opt.xhat{i})(m,:).std]),...
            'Color', fig.c{m},...
            'LineStyle', 'none',...
            'Marker', fig.m{m},...
            'LineWidth', 1.0,...
            'MarkerSize', 8,...
            'MarkerFaceColor', fig.f{m});
          handles{m+1} = col(handles{m+1});
          if strcmp(opt.meth{m},'VPM')
            text(tmp, col([summlog10.(opt.prof{p}).(opt.xhat{i})(m,:).mean]),...
              num2str(col(1:n.roi)),...
              'VerticalAlignment', 'top',...
              'HorizontalAlignment', 'left',...
              'FontSize', 16);
          end
        end
        
        hold off;
        axis tight;
        axis square;
        grid on;
        set(gca,...
          'XTick', fig.tick.(opt.xhat{i}),...
          'YTick', fig.tick.(opt.xhat{i}),...
          'XLim', minmax(fig.tick.(opt.xhat{i})),...
          'YLim', minmax(fig.tick.(opt.xhat{i})),...
          'XTickLabel', fig.ticklab.(opt.xhat{i}),...
          'YTickLabel', fig.ticklab.(opt.xhat{i}),...
          'Layer', 'top');
        tmp = ['NIST ', upper(opt.xhat{i}), ' NMR estimates (ms)'];
        xlabel(tmp, 'FontSize', 16);
        tmp = ['Candidate ', upper(opt.xhat{i}), ' estimates (ms)'];
        ylabel(tmp, 'FontSize', 16);
        tmp = [handles{2:end}];
        tmp = legend([handles{1} tmp(1,:)],...
          'Ideal', opt.meth{:},...
          'Location', 'SE');
        set(tmp, 'FontSize', 18);
        if header.cmp
          tmp = sprintf('%s,sl-%u,%s,plot',...
            opt.prof{p}, n.sl, opt.xhat{i});
          print('-depsc', tmp);
        end
      end
    end
  end
end
