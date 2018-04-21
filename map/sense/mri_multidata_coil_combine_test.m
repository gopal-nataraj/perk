% script mri_multidata_coil_combine_test.m
% test script for multi-dataset coil combination via mri_multidata_coil_combine.m
%
% copyright 2015, gopal nataraj, university of michigan
%
% version control
%   1.1     2015-03-19    original
%   1.2     2016-06-12    testing noise amplification issues
%   1.3     2018-04-20    now generating data using spgr_fun_v2(...) and dess_fun_v2(...)

%% Load true parameter maps
% irt
if (~exist('irtdir', 'var'))
  curdir = cd('~/Box/work/irt'); 
  irtdir = pwd;
  setup();
  cd(curdir);
end

% mapping
addpath('../../model/spgr/');
addpath('../../model/dess/');
addpath('../../exp/sim/');

% Load phantom file
f.filename = 'phantom_1.0mm_msles2_crisp.fld';
labels = fld_read(f.filename); 

% Get true values
slice = 81; [nx, ny, nz] = size(labels(2:end, 2:end, :));
M0_true = double(squeeze(labels(2:end, 2:end, slice)));
T1_true = double(squeeze(labels(2:end, 2:end, slice)));
T2_true = double(squeeze(labels(2:end, 2:end, slice)));
T2s_true = double(squeeze(labels(2:end, 2:end, slice)));
for idx = 0:10
    f = mri_brainweb_params(idx);
    M0_true(M0_true == idx) = f.pd;
    T1_true(T1_true == idx) = f.t1; 
    T2_true(T2_true == idx) = f.t2;
    T2s_true(T2s_true == idx) = f.t2s;
end
wf_true = zeros(nx, ny);

% Create a true excitation profile map
% [xx, yy] = ndgrid(linspace(-0.5, 0.5, nx), linspace(-0.5, 0.5, ny));
% kap_true = 1.5 - 2*(xx.^2 + yy.^2); 
% kap_true(M0_true == 0) = 0; 
kap_true = ones(nx, ny);


%% Simulate a scan protocol
% Global values
TE = 4;                                 % ms
nc = 8;                                 % Number of SENSE coils

% SPGR imaging parameters
flip_s = []' * pi/180;                  % rad
nfs = length(flip_s);  
TRs = []';                              % ms
TEs = TE * ones(nfs, 1);                % ms

% DESS imaging parameters
flip_d = [15 45]' * pi/180;
nfd = length(flip_d);                   
TRd = [10 10]';
TEd = TE * ones(nfd, 1);                % ms (symmetric echoes)

% Total number of datasets 
M = length(flip_s) + 2*length(flip_d);

% Make a M0s_true map, where M0s = M0 * exp(-TE/T2s)
% This is to be used for generating the forward model 
T2s_msk = T2_true ~= T2s_true; 
M0s_true = M0_true;
M0s_true(T2s_msk) = M0_true(T2s_msk) .* exp(-TE ./ T2s_true(T2s_msk));


%% Synthesize Phantom Coil Data
% Simulate true coil sensitivities 
smap = mri_sensemap_sim('nx', nx, 'ny', ny, 'dx', 1, 'dy', 1, 'ncoil', nc, 'chat', 0);

% Forward model: make noiseless SPGR data
ys_true = NaN(nx, ny, nfs, nc);
xs_true = NaN(nx, ny, nfs);
for a = 1:nfs
    xs_true(:,:,a) = spgr_fun_v2(...
      M0s_true, T1_true, T2_true,...
      flip_s(a), TRs(a), TEs(a),...
      'kap', kap_true,...
      'dw', wf_true,...
      'mag', 0);
    for c = 1:nc
        ys_true(:,:,a,c) = fft2(smap(:,:,c) .* xs_true(:,:,a));
    end
end

% Forward model: make noiseless DESS data
yp_true = NaN(nx, ny, nfd, nc);
ym_true = NaN(nx, ny, nfd, nc);
xp_true = NaN(nx, ny, nfd);
xm_true = NaN(nx, ny, nfd);
for a = 1:nfd
    [xp_true(:,:,a), xm_true(:,:,a)] = dess_fun_v2(...
      M0s_true, T1_true, T2_true,...
      flip_d(a), TRd(a), TEd(a), TEd(a),...
      'kap', kap_true,...
      'dw', wf_true,...
      'mag', 0);
    for c = 1:nc
        yp_true(:,:,a,c) = fft2(smap(:,:,c) .* xp_true(:,:,a));
        ym_true(:,:,a,c) = fft2(smap(:,:,c) .* xm_true(:,:,a));
    end
end

% Concatenate all object data
x_true = cat(3, xs_true, xp_true, xm_true);

% Add complex white gaussian noise 
var_im = 1.31e-7; 
% var_im = 1e-5;
sig_k = sqrt(var_im * (nx*ny));
ys = ys_true + sig_k * (randn(size(ys_true)) + 1i * randn(size(ys_true)));
yp = yp_true + sig_k * (randn(size(yp_true)) + 1i * randn(size(yp_true)));
ym = ym_true + sig_k * (randn(size(ym_true)) + 1i * randn(size(ym_true)));
y = cat(3, ys, yp, ym);

% Since fully-sampled, can IFT back to image domain
ys_im = ifft2(ys);
yp_im = ifft2(yp);
ym_im = ifft2(ym);
y_im = cat(3, ys_im, yp_im, ym_im);


%% Estimate coil sensitivities and underlying object
[x_reg, s_reg, x_sos, s_ml, cost] = ...
    mri_multidata_coil_combine(y_im, 'normalize', 0, 'nouter', 10);

% Rotate true smap to have zero phase in 1st coil, relative phase after
% Compensate phase of x_true to reflect this modification
phase1 = angle(squeeze(smap(:,:,1)));
smap = smap .* exp(-1i * repmat(phase1, [1 1 nc]));
x_true = x_true .* exp(+1i * repmat(phase1, [1 1 M]));

% Compute differences
x_reg_dif = x_reg - x_true;
x_sos_dif = x_sos - x_true;
s_reg_dif = s_reg - smap;
s_ml_dif  = s_ml  - smap;

% Display estimates and absolute differences
figure(1); im('col', M, cat(4, abs(x_sos), abs(x_reg)), 'cbar', 'x-sos/x-reg');
figure(2); im('col', M, cat(4, abs(x_sos_dif), abs(x_reg_dif)),...
    'cbar', 'x-sos-diff/x-reg-diff');

figure(3); im('col', nc, cat(4, abs(s_ml),  abs(s_reg)), 'cbar', 's-ml/s-reg');
figure(4); im('col', nc, cat(4, abs(s_ml_dif), abs(s_reg_dif)),...
    'cbar', 's-ml-diff/s-reg-diff');
