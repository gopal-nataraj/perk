 function [x_reg, s_reg, x_sos, s_ml, cost] = mri_multidata_coil_combine(y_jmc, varargin)
%function (x_reg, s_reg, x_sos, s_ml, cost] = mri_multidata_coil_combine(y_jmc, [options])
%|
%|  Regularized estimation of (smooth) sensitivity maps for parallel MRI,
%|  particularly when (possibly) more than one dataset is available. 
%|
%|  This method avoids noise-amplifying "ratios" of conventional methods, 
%|  and smoothly interpolates over regions with signal voids. 
%|
%|  Signal model: y_jmc = s_jc f_jm + e_jmc
%|    s_jc    Complex sensitivity map (relative to "bodycoil" reference)
%|    f_jm    Unknown underlying object for mth dataset
%|    e_jmc   Complex-gaussian noise
%| 
%|  Signal model for "body coil" image: z_jm = f_jm + e_jm
%|    If a body coil image is provided, initialize f_jm with it.
%|    Otherwise, initialize f_jm with the square root of the sum-of-squares
%|    across coils of y_jmc * phase of first coil image, for each dataset m.
%|   
%|  Indexing 
%|    j       voxel index, 1,...,N
%|    m       dataset index, 1,...,M
%|    c       coil index, 1,...,C
%| 
%|  Inputs
%|    y_jmc       [(odims) M C]   Noisy complex coil images (2D or 3D)
%|
%|  Options
%|    bodycoil    [(odims) M]     Reference image                             def: SSoS
%|    sinit       [(odims) C]     Initial SENSE maps                          def: from ratio
%|    smask       [(odims)]       Object mask                                 def: true(odims)
%|    var         [1]             Noise variance estimate in coil data        def: 0
%|    coilweights [C]             Coil weights                                def: ones(C,1)
%|    dataweights [M]             Dataset weights                             def: ones(M,1)
%|    log2b       [1]             log_2 of SENSE map reg param                def: -5
%|    order       [1]             Regularization order for smap               def: 2 
%|    thresh      [1]             Fraction of magnitude max used for median
%|                                  initial value in "background"             def: 0.05
%|    fwhm        [1]             FWHM of init smap gauss filter (pixels)     def: 3
%|    nk          [1]             Number of gaussian filter samples           def: 9 
%|    nouter      [1]             Number of outer iterations                  def: 20
%|    ninner      [1]             Max number of smap iterations               def: 50
%|    log10tol    [1]             Stop smap inner iterations if 
%|                                  norm(s_new-s_old)/norm(s_new)<10^log10tol def: -5
%|    relPh       0|1             If on, init bodycoil phase est w/ coil1     def: 0
%|    chol        0|1             Toggle cholesky decomposition off|on        def: 1
%|    normalize   0|1             Toggle normalizing data before estimating
%|                                  smap, to make log2b scale-invariant       def: 1
%|    disp        0|1             Toggle display off|on                       def: 1
%|
%|  Outputs
%|    x_reg       [(odims) M]     Coil-combined image, regularized
%|    s_reg       [(odims) C]     Estimated sensitivity maps
%|    x_sos       [(odims) M]     Sqrt-SoS coil-combined initialization
%|    s_ml        [(odims) C]     Coil sensitivity initialization
%|    cost        [2*nouter+1]    Cost after each inner update
%|
%|  Other Notes:
%|    1)     Non-iterative cholesky decomposition is faster for small J,M
%|    2)     Regularization parameter is "universal" (scale-invariant) 
%| 
%|  Related works:
%|    1)     "Accelerated regularized estimation of MR coil sensitivities 
%|              using augmented Lagrangian methods"
%|              M. J. Allison et al, Mar. 2013 IEEE T-MI
%|    2)     "Smoothing effect of sensitivity map on fMRI data using a novel
%|              self-calibrated estimation method"
%|              Y. Kim et al, 2008 ISMRM abstract, p. 1267
%| 
%|  Copyright 2015, Gopal Nataraj, University of Michigan
%|
%|  Version Control
%|    1.1         2015-03-19      original
%|    1.2         2016-06-13      initial smaps gaussian-filtered
%|                                made optional: initializing body-coil reference w/ crude phase 

% Test case
% todo: package script test case into function with all helpers
if (nargin == 1 && streq(y_jmc, 'test'))
    [x_reg, s_reg] = mri_multidata_coil_combine_test;
    if (~nargout), clear('x_reg','s_reg'), end;
    return
end

% Default values 
arg.bodycoil = [];
arg.sinit = [];
arg.smask = [];
arg.var = 0;
arg.coilweights = [];
arg.dataweights = [];
arg.log2b = -5;
arg.order = 2;
arg.thresh = 0.05;
arg.fwhm = 3;
arg.nk = 9;
arg.nouter = 20;
arg.ninner = 50;
arg.log10tol = -5;
arg.relPh = false;
arg.chol = true;
arg.normalize = true;
arg.disp = true;

% Substitute varargin values as appropriate
arg = vararg_pair(arg, varargin); 

% Dimensions
tmp = size(y_jmc);
odims = tmp(1:end-2);                               % Object dimensions
M = tmp(end-1);                                     % Number of datasets
C = tmp(end);                                       % Number of coils
coildim = numel(tmp);

% If no mask specified, extrapolate to all voxels
if isempty(arg.smask)
    arg.smask = true(odims);
    N = prod(odims);                                % Number of voxels of interest
else
    N = numel(arg.smask(arg.smask));                % Number of voxels of interest
end

% If no coil weights given, set weights to be equal
if isempty(arg.coilweights)
    arg.coilweights = ones(C, 1);
end
Wc = ones(N,1) * arg.coilweights';                  % [N C]

% If no dataset weights given, set weights to be equal
if isempty(arg.dataweights)
    arg.dataweights = ones(M, 1);
end
Wm = Gdiag(col(repmat(arg.dataweights', [N 1])), 'mask',...
    repmat(arg.smask, [ones(1,length(odims)) M]));  % [N*M N*M]

% If no body coil reference given, set default to sqrt(sum-of-squares)
if isempty(arg.bodycoil)
    arg.bodycoil = sqrt(max(sum(abs(y_jmc).^2, coildim) - C*arg.var, 0));
    % optional: add crude, noisy phase estimate from first coil
    if arg.relPh
        tmp = angle(stackpick(y_jmc, 1));           % crude phase estimate
        arg.bodycoil = arg.bodycoil .* exp(1i*tmp); % [(odims) M]
    end
end

% Return the initial reference object if desired
x_sos = arg.bodycoil;                               % [(odims) M]
if (arg.disp)
    dispmax = max(col(abs(x_sos)));                 % Fix display window
    figure(1); 
    im(abs(x_sos), [0 dispmax], 'cbar'); 
    drawnow;
end

% Initialize SENSE maps using ML estimate and body coil reference
if isempty(arg.sinit)
    % Set all uncertain initial map values to median of good values
    % Good values are selected from dataset with highest weighting
    [~, m_best] = max(arg.dataweights);
    body_best = stackpick(arg.bodycoil, m_best);    % [(odims)]
    good = abs(body_best) > arg.thresh * max(abs(body_best(:)));
    
    % Initialize system matrix (same for all coils)
    tmp = cell(M, 1);
    for m = 1:M
        body_m = stackpick(arg.bodycoil, m);        % [(odims)]
        tmp{m} = Gdiag(body_m(arg.smask), 'mask', arg.smask);
    end
    Am = vertcat(tmp{:});                           % [N*M N]
    
    arg.sinit = zeros([N C]);                       % [N C]
    for c = 1:C
        y_jm = stackpick(y_jmc, c);                 % [(odims) M]
        y_jm = col(y_jm(repmat(arg.smask,...
            [ones(1, length(odims)) M])));          % [N*M 1]
        tmp = ((Am'*Wm) * y_jm) ./ ((Am'*Wm*Am) * ones(size(Am,2), 1));
        tmp(~good(arg.smask)) = ...
            median(abs(tmp(good(arg.smask))));      % [N]
        arg.sinit(:,c) = tmp;
    end
    
    % Smooth initial smaps via gaussian filtering along each direction
    tmp = ceil((arg.nk-1)/2);       
    kern = gaussian_kernel(arg.fwhm, tmp);          % [nk]
    for c = 1:C
        tmp = embed(arg.sinit(:,c), arg.smask);     % [(odims)]
        kern = col(kern);
        for d = 1:length(odims)
            tmp = convn(tmp, kern, 'same');
            kern = reshape(kern, [1 size(kern)]);
        end
        arg.sinit(:,c) = tmp(arg.smask);
    end
    arg.sinit = embed(arg.sinit, arg.smask);        % [(odims) C]
end

% Return the initial SENSE maps if desired
s_ml = arg.sinit;                                   % [(odims) C]
if (arg.disp)
    figure(2); 
    im(abs(s_ml), 'cbar'); 
    drawnow;
end
        
% Regularizer object for SENSE maps
reg_args = {'beta', 2^arg.log2b, 'order', arg.order, 'distance_power', 2};
if arg.chol
    R = Reg1(arg.smask, reg_args{:}, 'type_penal', 'mat', 'type_diff', 'spmat');
else
    R = Reg1(arg.smask, reg_args{:}, 'type_penal', 'mat', 'type_diff', 'sparse');
end
Cdiff = R.C;                                        % [K N]

% Trick: normalize data by median of non-background value in bodycoil image
% so that the effective regularization beta is scale-invariant
if (arg.normalize)
    tmp = abs(arg.bodycoil); 
    tmp = median(tmp(tmp > arg.thresh * max(tmp(:))));
    arg.bodycoil = arg.bodycoil / tmp; 
    y_jmc = y_jmc / tmp;
    if (arg.disp), dispmax = dispmax / tmp; end;
end

% Set x_reg initial estimate to body coil reference
x_reg = NaN(N, M);
for m = 1:M
    body_m = stackpick(arg.bodycoil, m);            % [(odims)]
    x_reg(:,m) = body_m(arg.smask);                 % [N]
end

% Set s_reg initial estimate to ML estimate
s_reg = NaN(N, C);
for c = 1:C
    sinit_c = stackpick(arg.sinit, c);              % [(odims)]
    s_reg(:,c) = sinit_c(arg.smask);                % [N]
end

% Reshape and extract data of interest
y_jmc = reshape(y_jmc(repmat(arg.smask,...
    [ones(1, length(odims)) M C])), [N M C]);       % [N M C]

% Compute initial cost
if (nargout > 4)
    cost = NaN(2*arg.nouter+1,1);
    cost_iter = 1;
    cost(cost_iter) = pwls_cost_allcoils(x_reg, s_reg, Wm, arg.coilweights, y_jmc, R); 
    cost_iter = cost_iter + 1; 
end

for outer = 1:arg.nouter
    % Display Outer Iteration 
    printm('Outer Iteration: %d of %d.', outer, arg.nouter);
    
    %% Dataset-by-dataset object update
    % Least-squares solution at each voxel
    den = sum(conj(s_reg) .* (Wc .* s_reg), 2);     % [N]
    for m = 1:M
        y_jc = squeeze(y_jmc(:,m,:));               % [N C]
        num = sum(conj(s_reg) .* (Wc .* y_jc), 2);  % [N]
        x_reg(:,m) = div0(num, den);                % [N]
    end
    
    % Compute cost after object update
    if (nargout > 4)
        cost(cost_iter) = pwls_cost_allcoils(x_reg, s_reg,...
            Wm, arg.coilweights, y_jmc, R);
        cost_iter = cost_iter + 1;
    end
    
    % Update object display if desired
    if (arg.disp)
        figure(1); 
        im(abs(embed(x_reg, arg.smask)), [0 dispmax], 'cbar');
        drawnow;
    end
    
    %% Coil-by-coil SENSE map update
    % Update the system matrices (same for all coils)
    tmp = cell(M, 1);
    for m = 1:M
        tmp{m} = Gdiag(x_reg(:,m), 'mask', arg.smask);
    end
    Am = vertcat(tmp{:});                           % [N*M N]
    
    % For medium-sized problems, fast to use numerical inverse
    if arg.chol
        Wm = sparse(Wm);                            % fatrix2 -> spmat
        Am = sparse(Am);                            % fatrix2 -> spmat
        Hm = (Am' * (Wm * Am)) + (Cdiff' * Cdiff);  % [N N]
        
        for c = 1:C
            y_jm = col(squeeze(y_jmc(:,:,c)));      % [N*M 1]
            tmp = Am' * y_jm;                       % [N]
            s_reg(:,c) = Hm \ tmp;                  % [N]
        end
    % For larger problems, use QPWLS iterative algorithm    
    else
        Dm = (Am' * (Wm * Am));                     % Datafit Hessian
        for c = 1:C
            y_jm = col(squeeze(y_jmc(:,:,c)));      % [N*M 1]
            s_c = squeeze(s_reg(:,c));              % [N]
            
            % Separate preconditioner for each coil
            % To be safe, take magnitude to avoid complex preconditioner
            Dmc = Dm + Gdiag(R.denom(R, s_c), 'mask', arg.smask);
            Pmc = Gdiag(1 ./ abs(Dmc * ones(size(Dmc,1), 1)));
            
            s_reg(:,c) = qpwls_pcg1(s_c, Am, Wm, y_jm, Cdiff,...
                'niter', arg.ninner, 'precon', Pmc,...
                'stop_diff_tol', 10^arg.log10tol, 'chat', arg.disp);
        end
    end
    
    % Compute cost after SENSE map update
    if (nargout > 4)
        cost(cost_iter) = pwls_cost_allcoils(x_reg, s_reg,...
            Wm, arg.coilweights, y_jmc, R);
        cost_iter = cost_iter + 1;
    end
    
    % Update SENSE map display if desired
    if (arg.disp)
        figure(2);
        im(abs(embed(s_reg, arg.smask)), 'cbar');
        drawnow;
    end
end
 
% Embed regularized estimates before returning to caller
x_reg = embed(x_reg, arg.smask);
s_reg = embed(s_reg, arg.smask);
 end
        
% Helper function to compute full cost across all coils
function cost = pwls_cost_allcoils(obj, smap, Wobj, ws, y, Rs)
    [~, M, C] = size(y);

    % Construct system fatrix
    tmp = cell(M, 1);
    for m = 1:M
        tmp{m} = Gdiag(obj(:,m), 'mask', Rs.mask);  % [N N]
    end
    Aobj = vertcat(tmp{:});                         % [N*M N]
    
    % Compute cost separately by coil
    cost_c = NaN(C, 1);
    for c = 1:C
        smap_c = squeeze(smap(:,c));                % [N]
        y_c = col(squeeze(y(:,:,c)));               % [N*M 1]
        cost_c(c) = pwls_cost(smap_c, Aobj, Wobj, y_c, Rs);
    end
    
    % Sum the cost across the coils, weighted by coil weights
    cost = ws.' * cost_c;
end
