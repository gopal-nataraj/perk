 function [b1_init, b1_reg, cost] = mri_bsb1map(yp_coil, ym_coil, varargin)
%function [b1_init, b1_reg, cost] = mri_bsb1map(yp_coil, ym_coil, [options])
%|
%| Regularized estimation of (smooth) |b1+| maps from multiple
%| receive-coil Bloch-Siegert data (single transmit channel)
%|
%| This method avoids noise-amplifying "ratios" of conventional methods,
%| and smoothly interpolates over regions with signal voids. 
%|
%| Inputs
%|  yp_coil         [(odims) nRc]   Noisy coil data at +freq offset
%|  ym_coil         [(odims) nRc]   Noisy coil data at -freq offset
%|
%| Options
%|  ph_init         [(odims)]       Initial phase estimate  (def: from ratio)
%|                                  ph_init = Kbs*|B1+|^2
%|  mask            [(odims)]       Image mask              (def: true(odims))
%|  Kbs             [1]             BS constant, rad/G^2    (def: 111.0488)
%|  log2b_coil      [1]             log2(coil reg param)    (def: -5
%|  log2b_b1        [1]             log2(b1 reg parameter)  (def: 3)
%|  iter_b1         [1]             B1-update iterations    (def: 50)
%|  log10tol        [1]             log10(conv criterion)   (def: -5)
%|  disp            0|1             Toggle display on/off   (def: 1)
%|  disp_rng        [2]             Iteration display range (def: [0 0.075])
%|  coilOpt         {1 2*nOpt}      Coil-combine options    (def: {'nouter', 10, 'log2b', -5})
%|
%| Outputs
%|  b1_init         [(odims)]       Initial b1 estimate
%|  b1_reg          [(odims)]       Regularized b1 estimate
%|  cost            [niter+1]       Cost after each update
%| 
%| Other notes:
%|  1)  Since this method depends heavily on initialization, coil data
%|      is first coil-combined using mri_multidata_coil_combine().
%|  2)  This method regularizes phi \propto |B1+|^2, which may be
%|      undesirable since regions with higher |B1+| are regularized more.
%|
%| Related work:
%|  1)  "Regularized Estimation of Bloch-Siegert |B1+| maps in MRI"
%|      H. Sun et al, Proc. Intl. Symp. Im. Proc., 2014
%|
%| Written by: Gopal Nataraj
%| Copyright 2015, University of Michigan
%|
%| Version control
%|  2015-09-09      original
%|  2016-05-17      Added optional arguments to mri_multidata_coil_combine(...)

% Default values
arg.ph_init = [];
arg.mask = [];
arg.Kbs = 111.0488;
arg.log2b_b1 = 3;
arg.iter_b1 = 50;
arg.log10tol = -5;
arg.disp = 1;
arg.disp_rng = [0 0.075];
arg.coilOpt = {'nouter', 10, 'log2b', -5};

% Substitute varargin values as appropriate
arg = vararg_pair(arg, varargin);

% Dimensions
tmp = size(yp_coil);
odims = tmp(1:end-1);                               % Object dimensions
coildim = numel(tmp);

% If no mask specified, extrapolate to all voxels
if isempty(arg.mask)
    arg.mask = true(odims);
    N = prod(odims);                                % Number of voxels of interest
else 
    N = numel(arg.mask(arg.mask));                  % Number of voxels of interest
end

% Combine the coil data 
y_coil = permute(cat(coildim+1, yp_coil, ym_coil), [1:length(odims) coildim+1 coildim]);
[y_im, ~] = mri_multidata_coil_combine(...
    y_coil, arg.coilOpt{:});                        % [(odims) 2]
yp_im = stackpick(y_im, 1);                         % [(odims)]
ym_im = stackpick(y_im, 2);                         % [(odims)]

% If no ph_init specified, compute via conventional estimate
if isempty(arg.ph_init)
    arg.ph_init = angle(yp_im .* conj(ym_im)) / 2;  % rad
    arg.ph_init(arg.ph_init<0) = ...
        arg.ph_init(arg.ph_init<0) + 2*pi;          % ensure positive phase
end

% Columnize data and initial estimate
yp_im = col(yp_im(arg.mask));                       % [N]
ym_im = col(ym_im(arg.mask));                       % [N]
ph_reg = col(arg.ph_init(arg.mask));                % [N]

% Output initial b1 guess
b1_init = sqrt(arg.ph_init / arg.Kbs);              % Gauss
if (arg.disp)
    figure(1); 
    im(b1_init, arg.disp_rng, 'cbar');
    drawnow;
end

% Precompute fixed constants
mag_pp = abs(yp_im).^2;                             % [N]
mag_mm = abs(ym_im).^2;                             % [N]
mag_pm = abs(yp_im .* ym_im);                       % [N]
alpha = angle(yp_im) - angle(ym_im);                % [N]

% Define regularizer object
reg_args = {'beta', 2^arg.log2b_b1, 'distance_power', 2, 'type_penal', 'mat'};
R = Reg1(arg.mask, reg_args{:}, 'type_diff', 'sparse');

% Define fixed preconditioner
d = 2*mag_pm + R.denom(R, ph_reg);                  % [N]
P = Gdiag(1 ./ d, 'mask', arg.mask);                % [N N]   
    
% Compute the initial cost
if (nargout > 2)
    cost = NaN(arg.iter_b1+1,1);
    cost_iter = 1;
    cost(cost_iter) = fullcost(ph_reg, mag_pp, mag_mm, mag_pm, alpha, R);
    cost_iter = cost_iter + 1;
end

%% PGD Iterations to update ph_reg
for iter = 1:arg.iter_b1
    % Display outer iteration
    if(rem(iter,10) == 0), printm('Iteration: %d of %d.', iter, arg.iter_b1); end;
    
    % Compute gradient
    ph_reg_prev = ph_reg;
    grad = mag_pm .* sin(2*ph_reg_prev - alpha) + R.cgrad(R, ph_reg_prev);
    
    % PGD update
    ph_reg = ph_reg_prev - (P * grad);
    
    % Compute new cost
    if (nargout > 2)
        cost(cost_iter) = fullcost(ph_reg, mag_pp, mag_mm, mag_pm, alpha, R);
        cost_iter = cost_iter + 1;
    end
    
    % Update display if desired
    if (arg.disp && (rem(iter,10) == 0))
        figure(1);
        im(embed(sqrt(ph_reg/arg.Kbs), arg.mask), arg.disp_rng, 'cbar');
        drawnow;
    end
    
    % Exit early if tolerance criterion satisfied
    if ((norm(ph_reg - ph_reg_prev)/norm(ph_reg)) < 10^arg.log10tol)
        printm('Exited at iteration %u of %u', iter, arg.iter_b1);
        break;
    end
end

%% Embed regularized estimate before returning to caller
b1_reg = embed(sqrt(ph_reg/arg.Kbs), arg.mask);
 end
 
 % Helper function to compute full cost
 function cost = fullcost(ph_reg, mag_pp, mag_mm, mag_pm, alpha, R)
    cost = (1/4) * sum(mag_pp.^2 + mag_mm.^2 ...
        + 2*mag_pm.*cos(2*ph_reg - alpha)) ...
        + R.penal(R, ph_reg);
 end
    