% script mri_bs_kapb1_map_test.m
% test script for regularized b1 estimation
%
% copyright 2015, gopal nataraj, university of michigan
% 
% version control
%   1.1       2015-07-17      original
%   1.2       2016-06-10      modified for v1.3
%   1.3       2016-06-13      strong gauss filtering of init sense maps 
%   1.4       2018-04-20      now generating data using spgr_fun_v2(...)

%% Load true parameter maps and set scan protocol
% irt
if (~exist('irtdir', 'var'))
  curdir = cd('~/Box/work/irt');
  irtdir = pwd;
  setup();
  cd(curdir);
end

% mapping
addpath('../../model/spgr/');
addpath('../../map/sense/');
addpath('../../exp/sim/');

% Load phantom file
f.filename = 'phantom_1.0mm_msles2_crisp.fld';
labels = fld_read(f.filename); 

% Scan parameters
nx = 256; ny = 256; nz = 6; nc = 32;
FOVx = 240; FOVy = 240; FOVz = 30;          % mm
dx = FOVx/nx; dy = FOVy/ny; dz = FOVz/nz;   % mm
voxelVol = dx*dy*dz;                        % mm^3

b0 = 3.0;                                   % Tesla
flip_s = 10 * pi/180;                       % radians
nfs = length(flip_s);                       % ms
TRs = 39;                                   % ms
TEs = 12;                                   % ms (long, need to check sequence exactly)

% Get true m0, t1, t2, t2s values
slice = 81; [nx, ny, nz] = size(labels);
M0_true = double(squeeze(labels(:,:,slice)));
T1_true = double(squeeze(labels(:,:,slice)));
T2_true = double(squeeze(labels(:,:,slice)));
T2s_true = double(squeeze(labels(:,:,slice)));
for idx = 0:10
    f = mri_brainweb_params(idx, 'b0', b0);
    M0_true(M0_true == idx) = f.pd;
    T1_true(T1_true == idx) = f.t1; 
    T2_true(T2_true == idx) = f.t2;
    T2s_true(T2s_true == idx) = f.t2s;
end
wf_true = zeros(size(M0_true)); 

% Make a M0s_true map, where M0s = M0 * exp(-TE/T2s)
% This is to be used for generating the forward model 
T2s_msk = M0_true ~= 0; 
M0s_true = M0_true;
M0s_true(T2s_msk) = M0_true(T2s_msk) .* exp(-TEs ./ T2s_true(T2s_msk));

% Create excitation profile map
% Cannot allow too much variation in flip because then ML T1/T2 estimate biased
kap_mean = 0.8;                             % flip angle mean
kap_range = 0.1;                            % flip angle variation
[xx, yy] = ndgrid(linspace(-1, 1, nx), linspace(-1, 1, ny));
kap_true = kap_mean * ((1+kap_range/2) - kap_range*(xx.^2 + yy.^2));

% Create true b1 map
b1_peak = 0.075;                            % Gauss
b1_true = b1_peak * kap_true;               % Latent b1 map, magnitude only

% Create true coil sensitivities
nRc = 8;                                % Number of rX coils
smap = mri_sensemap_sim('nx', nx, 'ny', ny, 'ncoil', nRc);


%% Synthesize Phantom BS SPGR Data
% Create Fermi pulse (B1_peak preset)
w_rf = 8000;                                % Hz
dt = 4e-6;                                  % s
tot_dur = 0.009;                            % s
pk_dur = 0.004;                             % s
tr_dur = 0.0001;                            % s

t = 0:dt:tot_dur;                           % s
env = 1 ./ (1 + exp((abs(t-tot_dur/2) - pk_dur)/tr_dur));    
b1 = b1_peak * env .* exp(1i*2*pi*w_rf*t);  % Gauss

% Calculate Kbs
gam = 4.2576e3;                             % Hz/G
w_bs = (gam*abs(b1)).^2/(2*w_rf);             % Hz
phi_bs = trapz(t, (2*pi*w_bs));             % rad
Kbs = phi_bs / max(abs(b1).^2);             % rad/Gauss^2

% Compute phase offset
ph = exp(1i * Kbs * abs(b1_true).^2);

% Forward model: create coil-by-coil data
yp_coil_true = NaN(nx, ny, nRc);
ym_coil_true = NaN(nx, ny, nRc);
for c = 1:nRc
  yp_coil_true(:,:,c) = smap(:,:,c) .* ph .* spgr_fun_v2(...
    M0s_true, T1_true, T2_true,...
    flip_s, TRs, TEs,...
    'kap', kap_true,...
    'dw', wf_true,...
    'mag', 1);
  ym_coil_true(:,:,c) = smap(:,:,c) .* conj(ph) .* spgr_fun_v2(...
    M0s_true, T1_true, T2_true,...
    flip_s, TRs, TEs,...
    'kap', kap_true,...
    'dw', wf_true,...
    'mag', 1);
end

% Add complex white Gaussian noise
var_im = 1.31e-7;                           % Measured for 1mm x 1mm x 3mm voxel
var_im = var_im * ((1*1*3)/voxelVol)^2;     % Rescaled for our voxel size
sig_im = sqrt(var_im);
y.p = yp_coil_true + sig_im * (randn(size(yp_coil_true)) + 1i * randn(size(yp_coil_true)));
y.m = ym_coil_true + sig_im * (randn(size(ym_coil_true)) + 1i * randn(size(ym_coil_true)));


%% Regularized b1 map estimation
% Create tight and loose masks
tight_mask = imfill(T2s_msk, 'holes');
tight_mask = imdilate(~imdilate(~tight_mask, strel('disk', 5)), strel('disk', 5));
loose_mask = imdilate(tight_mask, strel('disk', 3)); 

bs.mag = transpose(abs(b1));
bs.wrf = w_rf / 1000;
bs.dt = dt * 1000;

% Mapping
mapArg = {...
  'mask', loose_mask,...
	'coilOpt', {'nouter', 3, 'fwhm', 5},...
  'stop.iter', 1000,...
  'reg.beta', 2^2,...
  'bool.chat', true,...
  'bool.disp', true};
x = mri_bs_kapb1_map(y, bs, mapArg{:});

% Display output
disp_range = b1_peak * kap_mean * [1-2*kap_range 1+2*kap_range];
err_range  = [0 b1_peak/10];

figure, im(cat(3, x.b1.mom, x.b1.rls, b1_true), disp_range, 'cbar');
figure, im(cat(3, x.kap.mom, x.kap.rls, kap_true), disp_range/b1_peak, 'cbar');

figure, im(cat(3, abs(x.b1.mom - b1_true), abs(x.b1.rls - b1_true)), err_range, 'cbar');
figure, im(cat(3, abs(x.kap.mom - kap_true), abs(x.kap.rls - kap_true)), err_range/b1_peak, 'cbar');

