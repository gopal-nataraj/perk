%% Test Script with Digital Phantom Data as to whether BS B1 mapping works
% This ensures that regularized B1 estimation from multiple-receive coil
% Bloch-Siegert data works in the "inverse crime" case
%
% Written by: Gopal Nataraj
% Originally created 2015-07-17

%% Load true parameter maps and set scan protocol
% Load digital phantom
if (~exist('irtdir', 'var'))
    curdir = pwd; 
    cd ../../../irt; 
    setup(); 
    cd(curdir);
end

addpath('../../data/DigitalPhantom');
addpath('../../model/bs/');
addpath('../../model/spgr/');
addpath('../../map/sense/');

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
wf_true = 0; 

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
w_rf = 2*pi*4000;                           % rad/s
dt = 4e-6;                                  % s
tot_dur = 0.009;                            % s
pk_dur = 0.004;                             % s
tr_dur = 0.0001;                            % s

t = 0:dt:tot_dur;                           % s
env = 1 ./ (1 + exp((abs(t-tot_dur/2) - pk_dur)/tr_dur));    
b1 = b1_peak * env .* exp(1i*w_rf*t);       % Gauss

% Calculate Kbs
gam = 4.2576e3 * (2*pi);                    % rad/Gauss
Kbs = (gam^2/(2*w_rf)) * ...
    trapz(t, env.^2);                       % rad/Gauss^2

% Compute phase offset
ph = exp(1i * Kbs * abs(b1_true).^2);

% Forward model: create coil-by-coil data
yp_coil_true = NaN(nx, ny, nRc);
ym_coil_true = NaN(nx, ny, nRc);
for c = 1:nRc
    yp_coil_true(:,:,c) = smap(:,:,c) .* (spgr_fun(M0s_true, T1_true, kap_true,...
        flip_s, TRs, TEs, wf_true, 1) .* ph);
    ym_coil_true(:,:,c) = smap(:,:,c) .* (spgr_fun(M0s_true, T1_true, kap_true,...
        flip_s, TRs, TEs, wf_true, 1) .* conj(ph));
end

% Add complex white Gaussian noise
var_im = 1.31e-7;                           % Measured for 1mm x 1mm x 3mm voxel
var_im = var_im * ((1*1*3)/voxelVol)^2;     % Rescaled for our voxel size
sig_im = sqrt(var_im);
yp_coil = yp_coil_true + sig_im * (randn(size(yp_coil_true)) + 1i * randn(size(yp_coil_true)));
ym_coil = ym_coil_true + sig_im * (randn(size(ym_coil_true)) + 1i * randn(size(ym_coil_true)));


%% Regularized b1 map estimation

% Create tight and loose masks
tight_mask = imfill(T2s_msk, 'holes');
tight_mask = imdilate(~imdilate(~tight_mask, strel('disk', 5)), strel('disk', 5));
loose_mask = imdilate(tight_mask, strel('disk', 3)); 

coil_iter = 20;
b1_iter = 1000;
disp_range = b1_peak * kap_mean * [1-kap_range/2 1+kap_range/2];
err_range  = [0 b1_peak/10];
log2b_coil = 10;                            % Strong coil reg to avoid ssos noise spikes
log2b_b1 = 5;                               % Higher requires better masking

b1Arg = {...
  'mask', loose_mask,...
  'Kbs', Kbs,...
  'log2b_b1', log2b_b1,...
  'iter_b1', b1_iter,...
  'disp_rng', disp_range,...
  'coilOpt', {'log2b', log2b_coil, 'nouter', coil_iter}};
[b1_init, b1_reg, cost] = mri_bsb1map(yp_coil, ym_coil, b1Arg{:});

% Display output
figure, im(b1_init, disp_range, 'cbar');
figure, im(b1_reg,  disp_range, 'cbar');
figure, im(b1_true, disp_range, 'cbar');

figure, im(abs(b1_init - b1_true), err_range, 'cbar');
figure, im(abs(b1_reg  - b1_true), err_range, 'cbar');

