  function [x] = mri_bs_kapb1_map(y, bs, varargin)
%|function [x] = mri_bs_kapb1_map(y, bs, varargin)
%|
%|  regularized estimation of kappa (flip scaling) |b1+| field maps
%|    assumes a single transmit channel
%|    smoothly interpolates over regions with signal voids
%|    avoids noise-amplifying 'ratios' of conventional methods
%|
%|  inputs
%|    y         [1x1 struct]      complex coil image data
%|     .p       [(odims) (C)]       at +freq offset
%|     .m       [(odims) (C)]       at -freq offset
%|    bs        [1x1 struct]      contains bloch-siegert pulse info
%|     .mag     [nt]                magnitude of off-resonant pulse                           G
%|     .wrf     [1]                 rf frequency offset                                       kHz
%|     .dt      [1]                 dwell time                                                ms
%|
%|  options
%|    mask      [(odims)]         object mask                           def: true(odims)
%|    x0        [(odims)]         initial estimates                                    
%|     .ph      [(odims)]           phase = Kbs*|B1+|^2                 def: from ratio       G
%|    coilOpt   {1 2*nOpt}        coil-combine options                  def: see below
%|    scale     [1]               global kap scale factor               def: 1
%|    reg       [1x1 struct]      regularizer options                   
%|     .pot     {1 npotarg}         potential function arguments        def: {'quad'}
%|     .beta    [1]                 strength                            def: 10^1
%|    stop      [1x1 struct] 
%|     .iter    [1]               maximum number of iterations          def: 1000
%|     .tolx    [1]               compare vs norm(xnew-xold)/norm(xnew) def: 10^-7
%|    bool      [1x1 struct]      boolean variables
%|     .chat    false|true          verbosity                           def: false
%|     .disp    [(odims)]           show image updates                  def: false
%| 
%|  outputs
%|    x         [1x1 struct]      object estimates (b1,kap)
%|     .*.mom   [(odims)]           method-of-moments, if computed
%|     .*.rls   [(odims)]           regularized least-squares 
%|
%| other notes:
%|  1)  Since this method depends heavily on initialization, coil data
%|      is first coil-combined using mri_multidata_coil_combine().
%|  2)  This method regularizes phi \propto |B1+|^2, which may be
%|      undesirable since regions with higher |B1+| are regularized more.
%|
%| related work:
%|  1)  "Regularized Estimation of Bloch-Siegert |B1+| maps in MRI"
%|      H. Sun et al, Proc. Intl. Symp. Im. Proc., 2014
%|
%| Written by: Gopal Nataraj
%| Copyright 2015, University of Michigan
%|
%| version control
%|    1.1       2015-09-09      original
%|    1.2       2016-05-17      added optional arguments to mri_multidata_coil_combine(...)
%|    1.3       2016-06-10      bloch-siegert constant computation
%|    1.4       2016-06-13      non-negativity constraint on ph -> now pgpm

% default values
arg.mask = [];
arg.x0.ph = [];
arg.coilOpt = {'nouter', 10, 'log2b', -5};
arg.scale = 1;
arg.reg.pot = {'quad'};
arg.reg.beta = 10^1;
arg.stop.iter = 1000;
arg.stop.tolx = 10^-7;
arg.bool.chat = false;
arg.bool.disp = false;

% substitute varargin values as appropriate
arg = vararg_pair(arg, varargin);

% dimensions
tmp = size(y.p);
odims = tmp(1:end-1);
cdim = numel(tmp);

% compute bloch-siegert constant
bs.gam = 4.2576e3;                                                      % Hz/G
bs.wrf = bs.wrf * 1000;                                                 % Hz
bs.t = (bs.dt/1000) * (0:length(bs.mag)-1)';                            % s
bs.wbs = (bs.gam*bs.mag).^2 / (2*bs.wrf);                               % Hz
bs.phi = trapz(bs.t, (2*pi*bs.wbs));                                    % rad
bs.max = max(bs.mag);                                                   % G
bs.kbs = bs.phi / max(bs.max.^2);                                       % rad/G^2
if arg.bool.chat
  printm('computed Bloch-Siegert constant: %0.4f rad/G^2.', bs.kbs);            
end

% if no mask specified, extrapolate to all voxels
if isempty(arg.mask)
  arg.mask = true(odims);
end

% coil combine data
tmp = permute(cat(cdim+1, y.p, y.m), [1:length(odims) cdim+1 cdim]);    % [(odims) 2 C]
[tmp, ~] = mri_multidata_coil_combine(tmp, arg.coilOpt{:});             % [(odims) 2]
ycc.p = stackpick(tmp,1);                                               % [(odims)]
ycc.m = stackpick(tmp,2);                                               % [(odims)]

% initial estimate
if isempty(arg.x0.ph)
  arg.x0.ph = angle(ycc.p .* conj(ycc.m)) / 2;                          % [(odims)] rad
  arg.x0.ph(arg.x0.ph<0) = arg.x0.ph(arg.x0.ph<0) + 2*pi;               % ensure positive phase
 
  % save for output
  x.b1.mom = embed(masker(sqrt(arg.x0.ph / bs.kbs), arg.mask), arg.mask);
  x.kap.mom = embed(masker(x.b1.mom / bs.max * arg.scale, arg.mask), arg.mask); 
  if arg.bool.disp
    figure(1), im(x.b1.mom, [0 bs.max], 'cbar'), drawnow;
    figure(2), im(x.kap.mom, [0 1], 'cbar'), drawnow;
  end
end

% vectorize
ycc.p = ycc.p(arg.mask);                                                % [V]
ycc.m = ycc.m(arg.mask);                                                % [V]
ph = arg.x0.ph(arg.mask);                                               % [V]

% define regularizer
regArg = {...
  'pot_arg', arg.reg.pot,...
  'beta', arg.reg.beta,...
  'type_penal', 'mat',...
  'type_diff', 'sparse'};
R = Reg1(arg.mask, regArg{:});

% precompute some constants
c.pp = abs(ycc.p).^2;                                                   % [V]
c.mm = abs(ycc.m).^2;                                                   % [V]
c.pm = abs(ycc.p .* ycc.m);                                             % [V]
c.ang = angle(ycc.p) - angle(ycc.m);                                    % [V]

% define fixed preconditioner
d = 2*c.pm + R.denom(R, ph);
P = Gdiag(1 ./ d, 'mask', arg.mask);                                    % [V V fatrix2]   

% pgd iterations
for i = 1:arg.stop.iter
  if arg.bool.chat && rem(i,10)==0
    printm('pgd cost %0.6e at iter %u of %u.', Psi(ph, c, R), i, arg.stop.iter);
  end

  % evaluate gradient
  ph_prev = ph;
  g = c.pm .* sin(2*ph_prev - c.ang) + R.cgrad(R, ph_prev);
  
  % projected update
  ph = ph_prev - (P * g);
  ph = max(ph, 0);
  
  % display update
  if arg.bool.disp && rem(i,100)==0
    tmp = sqrt(ph / bs.kbs);
    figure(1), im(embed(tmp, arg.mask), 'cbar'), drawnow;
    figure(2), im(embed(tmp / bs.max * arg.scale, arg.mask), 'cbar'), drawnow;
  end
  
  % exit early if tolerance criterion satisfied
  if norm(ph-ph_prev)/norm(ph) < arg.stop.tolx
    if arg.bool.chat
      printm('return early: norm(xnew-xold)/norm(xnew) < xTol = %0.3e.\n', arg.stop.tolx);
    end
    break;
  end
end

% embed for output
tmp = sqrt(ph / bs.kbs);
x.b1.rls = embed(tmp, arg.mask);
x.kap.rls = embed(tmp / bs.max * arg.scale, arg.mask);
end


% cost function evaluation
function cost = Psi(ph, c, R)
  cost = sum(c.pp.^2 + c.mm.^2 + 2*c.pm.*cos(2*ph - c.ang));
  cost = cost/4 + R.penal(R, ph);
end
