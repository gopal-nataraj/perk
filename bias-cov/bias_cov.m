  function [stat] = bias_cov(x, P, varargin)
%|function [stat] = bias_cov(x, P, varargin)
%|
%|  krr conditional bias-covariance calculations for m0,t1,t2(,inveff) estimation 
%|    handles coil-combined image data from the following pulse sequences:
%|      (single) spin-echo inversion recovery
%|        assume perfect spoiling immediately following (adiabatic) inversion
%|        since adiabatic inversion depends on t1/t2, jointly estimates inveff
%|      (single) spin-echo
%|      spoiled gradient-recalled echo
%|        assume perfect spoiling
%|      dual-echo steady state
%|        neglect diffusion effects
%|    permits b1 and b0 variation, imperfect inversion
%|    assumes single-component t1/t2 models
%|
%|  inputs
%|    x         [1x1 struct]    true latent parameters
%|     .*       [(odims)]         1st field: latent parameter (m0,t1,t2(,inveff))
%|    P         [1x1 struct]    scan parameters
%|                                1st field: data type (ir,se,sp,de)
%|                                2nd field: scan parameter, as appropriate
%|     .*.tr    [S.*]               repetition times                                          ms
%|     .*.ti    [S.ir]              inversion times                                           ms
%|     .*.te    [S.* E.*]           echo times                                                ms
%|     .*.ainv  [S.ir]              nominal effective flip angle of inversion                 rad
%|     .*.aex   [S.*]               nominal flip angle of excitation                          rad
%|     .*.aref  [S.*]               nominal flip angle of refocusing                          rad
%|
%|  options
%|    mask      [1x1 struct]    binary object masks 
%|     .disp    [(odims)]         over which to display images          def: true(odims)
%|     .est     [(odims)]         over which to perform calculations    def: imdilate(mask.disp)
%|    nu        [1x1 struct]    known parameters
%|     .kap     [(odims)]         flip-angle scale map                  def: ones(odims)            
%|     .b0      [(odims)]         off-resonance map                     def: zeros(odims)     kHz
%|     .r2p     [(odims)]         broadening linewidth map              def: zeros(odims)     kHz
%|    wght      [1x1 struct]    dataset weights
%|     .ir      [S.ir 1]          inversion-recovery weights            def: ones(S.ir,1)
%|     .se      [S.se 1]          spin-echo weights                     def: ones(S.se,1)
%|     .sp      [S.sp 1]          spoiled GRE weights                   def: ones(S.sp,1)
%|     .de      [S.de 2]          dual echo steady-state weights        def: ones(S.de,1)
%|    thresh    [1]             frac of max sig deciding 'background' 	def: 0.05 
%|    dist      [1x1 struct]    sampling distribution object (ignored if x0 set!)  
%|                                1st field: latent parameter (m0,t1,t2(,inveff))
%|     .*.supp  [2]               [lb ub] distribution support          def: see below
%|     .*.prior {1}               ('unif', 'logunif') distribution      def: see below 
%|    rff       [1x1 struct]    random fourier features object
%|     .snr     [1]               estimate of max sig for unity m0      def: 0.1
%|     .std     [1]               noise std dev in training data        def: 3.8607e-4
%|     .len.y   [D]               kernel data length scales             def: ones(D,1)
%|     .len.nu  [N]               kernel known param length scales      def: from nu
%|     .c       [1]               global kernel length scale parameter  def: 2^0
%|     .H       [1]               embedding dimension                   def: 10^4
%|     .K       [1]               number of training samples            def: 10^6
%|    inv       [1x1 struct]    matrix inversion reg strength object
%|     .krr     [1]               kernel ridge regression               def: 10^-8
%|    bool      [1x1 struct]    boolean variables
%|     .mag.*   false|true        using magnitude (ir,se,sp,de) data    def: false
%|     .chat    false|true        verbosity                             def: false
%|     .norm    false|true        normalize data for scale-invar reg    def: true
%|     .reset   false|true        reset rng while sampling              def: true
%|     .rfftst  false|true        show kernel approximation             def: false
%|     .disp    false|true        show image updates                    def: false     
%|
%|  outputs
%|    stat      [1x1 struct]    object statistics
%|     .bias    [(odims) L]       conditional bias
%|     .cov     [(odims) L L]     conditional covariance
%|
%|  copyright 2017, gopal nataraj, university of michigan
%|
%|  version control
%|    1.1       2017-09-01      original

% object dimensions
dim.odims = size(x.t1);
  
% number of scans and echoes
tmp = {'ir','se','sp','de'};
tmp2 = isfield(P, tmp);
for i = 1:length(tmp2)
  if tmp2(i)
    dim.S.(tmp{i}) = size(P.(tmp{i}).te,1);
    dim.E.(tmp{i}) = size(P.(tmp{i}).te,2);
  else
    dim.S.(tmp{i}) = 0;
    dim.E.(tmp{i}) = 0;
  end
end

% number of latent parameters
if dim.S.ir>0
  dim.L = 4;
else
  dim.L = 3;
end

% default values
arg.mask.est = [];
arg.mask.disp = [];

arg.nu.kap = [];
arg.nu.b0 = [];
arg.nu.r2p = [];

arg.wght.ir = ones(dim.S.ir, dim.E.ir);
arg.wght.se = ones(dim.S.se, dim.E.se);
arg.wght.sp = ones(dim.S.sp, dim.E.sp);
arg.wght.de = ones(dim.S.de, dim.E.de);

arg.thresh = 0.05;

arg.dist.m0.supp = [eps 1];
arg.dist.m0.prior = 'unif';
arg.dist.t1.supp = [400 2000].';
arg.dist.t1.prior = 'logunif';
arg.dist.t2.supp = [40 200].';
arg.dist.t2.prior = 'logunif';
if dim.L>3
  arg.dist.inveff.supp = [0.51 1].';
  arg.dist.inveff.prior = 'unif';
end

arg.rff.snr = 0.1;
arg.rff.std = [];
arg.rff.len.y = [];
arg.rff.len.nu = [];
arg.rff.c = 2^0;
arg.rff.H = 10^4;
arg.rff.K = 10^6;

arg.inv.krr = 10^-8;

arg.bool.mag.ir = false;
arg.bool.mag.se = false;
arg.bool.mag.sp = false;
arg.bool.mag.de = false;
arg.bool.chat   = false;
arg.bool.norm   = true;
arg.bool.reset  = true;
arg.bool.rfftst = false;
arg.bool.disp   = false;

% substitute varargin values as appropriate
arg = vararg_pair(arg, varargin);

% if no display mask specified, use all voxels
if isempty(arg.mask.disp)
  arg.mask.disp = true(dim.odims);
  arg.mask.est  = true(dim.odims);
% else if no estimation mask specified, dilate mask.disp
elseif isempty(arg.mask.est)
  arg.mask.est = imdilate(arg.mask.disp, strel('disk', 10));
% else make sure that display mask tighter than estimation mask
elseif sum(arg.mask.disp & ~arg.mask.est) > 0
  error('Display mask has voxels outside estimation mask?!');
end

% set V to number of voxels in mask.est
dim.V = numel(arg.mask.est(arg.mask.est));

% omit ir scans with zero weight
keep.ir = arg.wght.ir(:,1)>0;
if ~all(keep.ir)
  P.ir.tr = P.ir.tr(keep.ir);
  P.ir.ti = P.ir.ti(keep.ir);
  P.ir.te = P.ir.te(keep.ir, dim.E.ir);
  P.ir.aex = P.ir.aex(keep.ir);
  P.ir.aref = P.ir.aref(keep.ir);
  arg.wght.ir = arg.wght.ir(keep.ir,:);
  dim.S.ir = size(arg.wght.ir,1);
  if arg.bool.chat
    warn('\nDetected zero-weight IR data: omitting %u datasets.',...
      sum(~keep.ir)*dim.E.ir);
  end
  
  % if no ir scans with nonzero weight remain, remove inveff fields
  if dim.S.ir==0
    dim.L = 3;
    x = rmfield(x, 'inveff');
    arg.dist = rmfield(arg.dist, 'inveff');
    if arg.bool.chat
      warn('\nNo nonzero-weight IR scans remain; omitting inveff estimation.');
    end
  end
end

% omit se scans with zero weight
keep.se = arg.wght.se(:,1)>0;
if ~all(keep.se)
  P.se.tr = P.se.tr(keep.se);
  P.se.te = P.se.te(keep.se,:);
  P.se.aex = P.se.aex(keep.se);
  P.se.aref = P.se.aref(keep.se);
  arg.wght.se = arg.wght.se(keep.se,:);
  dim.S.se = size(arg.wght.se,1);
  if arg.bool.chat
    warn('\nDetected zero-weight SE data: omitting %u datasets.',...
      sum(~keep.se)*dim.E.se);
  end
end

% omit spgr scans with zero weight
keep.sp = arg.wght.sp(:,1)>0;
if ~all(keep.sp)
  P.sp.tr = P.sp.tr(keep.sp);
  P.sp.te = P.sp.te(keep.sp,:);
  P.sp.aex = P.sp.aex(keep.sp);
  arg.wght.sp = arg.wght.sp(keep.sp,:);
  dim.S.sp = size(arg.wght.sp,1);
  if arg.bool.chat
    warn('\nDetected zero-weight SPGR data: omitting %u datasets.',...
      sum(~keep.sp)*dim.E.sp);
  end
end

% omit dess scans with zero weight
keep.de = arg.wght.de(:,1)>0;
if ~all(keep.de)
  P.de.tr = P.de.tr(keep.de);
  P.de.te = P.de.te(keep.de,:);
  P.de.aex = P.de.aex(keep.de);
  arg.wght.de = arg.wght.de(keep.de,:);
  dim.S.de = size(arg.wght.de,1);
  if arg.bool.chat
    warn('\nDetected zero-weight DESS data: omitting %u datasets.',...
      sum(~keep.de)*dim.E.de);
  end
end

% total number of datasets
dim.Dir = dim.S.ir * dim.E.ir;
dim.Dse = dim.S.se * dim.E.se;
dim.Dsp = dim.S.sp * dim.E.sp;
dim.Dde = dim.S.de * dim.E.de;
dim.D   = sum([dim.Dir dim.Dse dim.Dsp dim.Dde]);   

% vectorize test x points
x = struct2cell(x);                                                     % {L} cell of [odims] arrays
for l = 1:dim.L
  x{l} = masker(x{l}, arg.mask.est);                                    % [V]
end

% if nu fields unspecified, set to default values
if isempty(arg.nu.kap)
  arg.nu.kap = ones(dim.odims);
end
if isempty(arg.nu.b0)
  if dim.S.de>0 && ~(arg.bool.mag.sp && arg.bool.mag.de) 
    warn('Using complex dess data w/o b0 map will likely induce bias due to phase-mismatch.');
  end
  arg.nu.b0 = zeros(dim.odims);
end
if isempty(arg.nu.r2p)
  if dim.S.de>0 && nnz(P.de.te(:,1)-P.de.te(:,2))>0
    warn('Using asymmetric dess echo times w/o r2p map will induce bias.');
  end
  arg.nu.r2p = zeros(dim.odims);
end

% vectorize nu fields
arg.nu = struct2cell(arg.nu);                                           % {N} cell of [odims] arrays
dim.N = length(arg.nu);
for n = 1:dim.N
  arg.nu{n} = masker(arg.nu{n}, arg.mask.est);                          % [V]
end

% vectorize data weights
w = [...
  col(arg.wght.ir);...
  col(arg.wght.se);...
  col(arg.wght.sp);...
  col(arg.wght.de)];                                                    % [D]

% to enumerate conveniently, switch from struct to cell
arg.dist = struct2cell(arg.dist);                                       % {L} cell of [1] structs

% check max sig for unity m0
if arg.rff.snr>1
  error('max unity-m0 signal cannot exceed 1!');
elseif arg.rff.snr<=eps
  error('max unity-m0 signal must be positive and less than 1!');
elseif arg.rff.snr>0.3
  tmp = strcat(...
    'Detected max unity-m0 signal >0.3: are you sure?\n',...
    'Overestimating may induce improper m0 sampling and krr errors...');
  warn(tmp);
elseif arg.rff.snr<0.01
  tmp = strcat(...
    'Detected max unity-m0 signal <0.01: are you sure?\n',...
    'Underestimating will cause inefficient m0 sampling and krr underperformance...');
  warn(tmp);
end

% set noise standard deviation
if isempty(arg.rff.std)
  warn('Detected empty rff.std; automatically setting to 3.8607e-4...');
  arg.rff.std = 3.8607e-4;
end

% set kernel lengthscales
if isempty(arg.rff.len.y)
  warn('Detected empty rff.len.y; automatically setting to ones(D,1)');
  arg.rff.len.y = ones(D,1);
end
if isempty(arg.rff.len.nu)
  arg.rff.len.nu = mean(transpose([arg.nu{:}]),2);
end
tmp = [arg.rff.len.y; arg.rff.len.nu];
arg.rff.len.q = max(tmp,eps);

% generate training points
samp = train(arg.rff, arg.dist, arg.nu, P, dim, arg.bool);

% bias and covariance
% conditioned on training samples, x, nu
stat = bias_cov_helper(x, arg.nu, P, w, samp, arg.rff, arg.inv.krr, dim, arg.bool);
      
% embed for output
stat.bias = embed(stat.bias, arg.mask.est);
stat.cov  = embed(stat.cov, arg.mask.est);

% use mask.disp for display
stat.bias = embed(masker(stat.bias, arg.mask.disp), arg.mask.disp);
stat.cov  = embed(masker(stat.cov, arg.mask.disp), arg.mask.disp);
  end

    
  function [samp] = train(rff, dist, nu, P, dim, bool)
%|function [samp] = train(rff, dist, nu, P, dim, bool)
%|
%|  kernel ridge regression training
%|
%|  inputs
%|    rff       [1x1 struct]    random fourier features object
%|     .snr     [1]               estimate of max sig for unity m0    
%|     .std     [1]               noise std dev in training data
%|     .len.y   [D]               kernel data length scales           
%|     .len.nu  [N]               kernel known param length scales     
%|     .c       [1]               global kernel length scale parameter 
%|     .H       [1]               embedding dimension       
%|     .K       [1]               number of training samples
%|    dist      {L cell}        sampling distribution object 
%|     .supp    [2]               [lb ub] distribution support        
%|     .prior   {1}               ('unif', 'logunif') distribution     
%|    nu        {N cell}        known parameters
%|    P         [1x1 struct]    scan parameters
%|                                1st field: data type (ir,se,sp,de)
%|                                2nd field: scan parameter, as appropriate
%|     .*.tr    [S.*]               repetition times                                          ms
%|     .*.ti    [S.ir]              inversion times                                           ms
%|     .*.te    [S.* E.*]           echo times                                                ms
%|     .*.ainv  [S.ir]              nominal effective flip angle of inversion                 rad
%|     .*.aex   [S.*]               nominal flip angle of excitation                          rad
%|     .*.aref  [S.*]               nominal flip angle of refocusing                          rad
%|    dim       [1x1 struct]    object containing dimension info
%|    bool      [1x1 struct]    boolean variables
%|     .reset   false|true        reset random number generator during training
%|     .exchg   false|true        estimate exchange map                
%|     .mag.*   false|true        using magnitude (spgr,dess) data
%|     .chat    false|true        verbosity  
%|     .rfftst  false|true        show kernel approximation
%|
%|  outputs
%|    samp      [1x1 struct]    krr training samples
%|     .X       [L K]             regressands
%|     .Q       [D+N K]           regressors
%|
%|  version control
%|    1.1       2017-06-06      adapted from mri_multicomp_map(...)
%|    1.2       2017-06-12      rff.snr now controls m0 distribution sampling

% optional: reset random number generator
if bool.reset
  % reset random number generator for reproducible results
  rng('default');
end

% sample x and set population means
xsamp = cell(dim.L,1);
for l = 1:dim.L
  switch dist{l}.prior
    case 'unif'
      xsamp{l} = random('unif', dist{l}.supp(1), dist{l}.supp(2), [rff.K, 1]);

      % population mean is (b-a)/2
      dist{l}.mean = mean(dist{l}.supp);
    case 'logunif'
      tmp = log(dist{l}.supp);
      tmp2 = random('unif', tmp(1), tmp(2), [rff.K, 1]);
      xsamp{l} = exp(tmp2);

      % population mean is (b-a)/(log(b)-log(a))
      dist{l}.mean = div0((dist{l}.supp(2)-dist{l}.supp(1)),tmp(2)-tmp(1));
    otherwise
      error('Unknown distribution prior option for latent parameter %u.', l);
  end
end
samp.X = transpose([xsamp{:}]);                                         % [L K]

% sample nu
nusamp = cell(dim.N,1);
for n = 1:dim.N
  if max(nu{n})-min(nu{n}) < eps
    % special case: delta function
    nusamp{n} = ones(rff.K,1) * mean(nu{n});
  else
    % kernel density estimate from known map
    kde = fitdist(nu{n}, 'kernel');

%     % code for including support later
%     tmp = masker(nu{n}, nu{n}>=0.5 & nu{n}<=2);
%     kde = fitdist(tmp, 'kernel', 'support', [0.5 2]);
    
    % sample from kde
    nusamp{n} = random(kde, [rff.K, 1]);
  end
end

% noiseless complex training data
tmp = bool;
tmp.mag.sp = 0;
tmp.mag.de = 0;
ysamp = f(xsamp, nusamp, P, dim, tmp);                                  % [DK]

% add noise and use magnitude data
tmp = randn(length(ysamp), 2); 
tmp = complex(tmp(:,1), tmp(:,2));
tmp = tmp * rff.std;
ysamp = ysamp + tmp;
ysamp = abs(ysamp);

% training inputs
ysamp = reshape(ysamp, [dim.D rff.K]);                                  % [D K]
tmp = transpose([nusamp{:}]);                                           % [N K]
samp.Q = [ysamp; tmp];                                                  % [D+N K]
  end
  
  
  function [stat] = bias_cov_helper(x, nu, P, w, samp, rff, reg, dim, bool)
%|function [stat] = bias_cov_helper(x, nu, P, w, samp, rff, reg, dim, bool)
%|
%|  bias_cov helper function
%|  
%|  inputs
%|    x         {L cell}        true latent parameters
%|    nu        {K cell}        known parameters
%|    P         [1x1 struct]    scan parameters
%|                                1st field: data type (ir,se,sp,de)
%|                                2nd field: scan parameter, as appropriate
%|     .*.tr    [S.*]               repetition times                                          ms
%|     .*.ti    [S.ir]              inversion times                                           ms
%|     .*.te    [S.* E.*]           echo times                                                ms
%|     .*.ainv  [S.ir]              nominal effective flip angle of inversion                 rad
%|     .*.aex   [S.*]               nominal flip angle of excitation                          rad
%|     .*.aref  [S.*]               nominal flip angle of refocusing
%|    w         [D]             dataset weights
%|    samp      [1x1 struct]    krr training samples
%|     .X       [L K]             regressands
%|     .Q       [D+N K]           regressors
%|    rff       [1x1 struct]    random fourier features object
%|     .snr     [1]               estimate of max sig for unity m0    
%|     .std     [1]               noise std dev in training data
%|     .len.y   [D]               kernel data length scales           
%|     .len.nu  [N]               kernel known param length scales     
%|     .c       [1]               global kernel length scale parameter 
%|     .H       [1]               embedding dimension       
%|     .K       [1]               number of training samples
%|    reg       [1]             krr regularization parameter
%|    dim       [1x1 struct]    object containing dimension info
%|    bool      [1x1 struct]    boolean variables
%|     .reset   false|true        reset random number generator during training
%|     .exchg   false|true        estimate exchange map                
%|     .mag.*   false|true        using magnitude (spgr,dess) data
%|     .chat    false|true        verbosity      
%|
%|  outputs
%|    stat      [1x1 struct]    object statistics
%|     .bias    [V L]             conditional bias
%|     .cov     [V L L]           conditional covariance
%|
%|  version control
%|    1.1       2017-09-04      implemented bias
%|    1.2       2017-09-05      implemented covariance
%|    1.3       2017-09-08      improved memory management

% trick: scale lengthscales to re-weight data
%   w>1 -> rff.len lower  -> data emphasized more
%   w<1 -> rff.len higher -> data emphasized less
rff.len.y = div0(rff.len.y, max(w,eps));

% lengthscale matrices
Lam.a   = diag(bsxfun(@max, rff.len.y, eps));                           % [D D]
Lam.nu  = diag(bsxfun(@max, rff.len.nu, eps));                          % [N N]
Lam.q   = blkdiag(Lam.a, Lam.nu);                                       % [D+N D+N]

% apply Gram matrix inversion separately to release memory
[krr, Xbar, Kbar] = krr_fun(samp, Lam, rff, reg);

% noiseless complex test signals
tmp = bool;
tmp.mag.sp = 0;
tmp.mag.de = 0;
ytest = f(x, nu, P, dim, tmp);                                          % [DV]
ytest = reshape(ytest, [dim.D dim.V]);                                  % [D V]

% conditional mean/covariance of noisy magnitude test signal (given x,nu)
tmp = abs(ytest).^2;
tmp = bsxfun(@plus, tmp, rff.std^2);
a.mu = sqrt(tmp);                                                       % [D V]
a.cov = diag(ones(dim.D,1) * rff.std^2);                                % [D D]

% conditional bias (given x,nu,X,Q)
tmp = Lam.a * (Lam.a^(-2)+a.cov^(-1))^(1/2) * a.cov^(1/2);
tmp = blkdiag(tmp, Lam.nu);
k = kernel(samp.Q, [a.mu; transpose([nu{:}])], tmp);                    % [K V]
k = k * sqrt(div0(1, det(a.cov/(Lam.a^2) + speye(dim.D))));
k = bsxfun(@minus, k, Kbar);                                            % [K V]
tmp2 = krr * k;                                                         % [L V]
tmp2 = bsxfun(@plus, tmp2, Xbar);                                       % [L V]
stat.bias = transpose(tmp2) - [x{:}];                                   % [V L]

% conditional covariance (given x,nu,X,Q)
knu = kernel(transpose([nu{:}]), samp.Q(dim.D+1:end,:), Lam.nu);        % [V K]
c1 = div0(1,sqrt(det(2*a.cov/(Lam.a^2) + speye(dim.D))));               % [1]
c2 = div0(1,det(a.cov/(Lam.a^2) + speye(dim.D)));                       % [1]
tmp2 = Lam.a*sqrt(2);                                                   % [D D]
tmp3 = tmp2 * (2*Lam.a^(-2)+a.cov^(-1))^(1/2) * a.cov^(1/2);            % [D D]
C1 = blkdiag(tmp2,tmp3);                                                % [2D 2D]
tmp3 = tmp2 * (Lam.a^(-2)+a.cov^(-1))^(1/2) * a.cov^(1/2);
C2 = blkdiag(tmp3,tmp3);                                                % [2D 2D]
stat.cov = NaN(dim.V, dim.L, dim.L);
for v = 1:dim.V
  tmp2 = bsxfun(@minus, samp.Q(1:dim.D,:), a.mu(:,v));                  % [D K]
  tmp3 = [tmp2; tmp2];                                                  % [2D K]
  tmp4 = [tmp2; -tmp2];                                                 % [2D K]
  tmp5 = c1*kernel(tmp3, tmp4, C1);                                     % [K K]
  tmp5 = tmp5 - c2*kernel(tmp3, tmp4, C2);                              % [K K]
  tmp6 = krr * spdiags(col(knu(v,:)), 0, rff.K, rff.K);                 % [L K]
  stat.cov(v,:,:) = tmp6 * (tmp5 * transpose(tmp6));                    % [L L]
end
  end
  
  
  function [krr, Xbar, Kbar] = krr_fun(samp, Lam, rff, reg)
%|function [krr, Xbar, Kbar] = krr_fun(samp, Lam, rff, reg)
%|
%|  applies kernel training data operations and then releases memory
%|
%|  inputs
%|    samp      [1x1 struct]    krr training samples
%|     .X       [L K]             regressands
%|     .Q       [D+N K]           regressors
%|    Lam       [1x1 struct]    smoothing lengthscale object
%|     .q       [D+N D+N]         for data and known parameters 
%|    rff       [1x1 struct]    random fourier features object
%|     .snr     [1]               estimate of max sig for unity m0    
%|     .std     [1]               noise std dev in training data
%|     .len.y   [D]               kernel data length scales           
%|     .len.nu  [N]               kernel known param length scales     
%|     .c       [1]               global kernel length scale parameter 
%|     .H       [1]               embedding dimension       
%|     .K       [1]               number of training samples
%|    reg       [1]             krr regularization parameter
%|
%|  outputs
%|    krr       [L K]           krr kernel training operations
%|    Xbar      [L]             sample mean of regressands
%|    Kbar      [K]             sample mean over columns of Gram matrix
%|
%|  version control
%|    1.1       2017-09-08      original

% demeaned regressands
Xbar = mean(samp.X,2);                                                  % [L]
XM = bsxfun(@minus, samp.X, Xbar);                                      % [L K]  
 
% apply kernel
K = kernel(samp.Q, samp.Q, Lam.q);                                      % [K K]
Kbar = mean(K,2);                                                       % [K]
K = bsxfun(@minus, K, Kbar);                                            % [K K]
krr = XM / (K + rff.K*reg*speye(rff.K));                                % [L K]
  end

  
  function k = kernel(q1, q2, Lam)
%|function k = kernel(q1, q2, Lam)
%|
%|  gaussian kernel evaluation
%|
%|  inputs
%|    q1        [Q N]           1st kernel argument
%|    q2        [Q M]           2nd kernel argument
%|    Lam       [Q Q]           smoothing lengthscale
%|
%|  outputs
%|    k         [N M]           kernel evaluations
%|
%|  version control
%|    1.1       2017-08-31      original

% bandwidth matrix can have near-zeros
warning('off', 'MATLAB:nearlySingularMatrix');
tmp = inv(Lam);
warning('on', 'MATLAB:nearlySingularMatrix');

k = NaN(size(q1,2), size(q2,2));
for i = 1:size(q1,2)
  dq = bsxfun(@minus, q1(:,i), q2);
  dq = tmp * dq;
  k(i,:) = exp(-(sum(abs(dq.^2),1))/2);
end
  end


  function [fval] = f(x, nu, P, dim, bool)
%|function [fval] = f(x, nu, P, dim, bool)
%|
%|  signal model function evaluation
%|
%|  inputs
%|    x         {L cell}        latent object parameters with cells size [V]
%|    nu        {N cell}        known parameters with cells size [Vloc]
%|    P         [1x1 struct]    scan parameters
%|                                1st field: data type (ir,se,sp,de)
%|                                2nd field: scan parameter, as appropriate
%|     .*.tr    [S.*]               repetition times                                          ms
%|     .*.ti    [S.ir]              inversion times                                           ms
%|     .*.te    [S.* E.*]           echo times                                                ms
%|     .*.ainv  [S.ir]              nominal effective flip angle of inversion                 rad
%|     .*.aex   [S.*]               nominal flip angle of excitation                          rad
%|     .*.aref  [S.*]               nominal flip angle of refocusing                          rad
%|    dim       [1x1 struct]    object containing dimension info  
%|    bool      [1x1 struct]    boolean variables               
%|     .mag.*   false|true        using magnitude (ir,se,sp,de) data     
%|
%|  outputs
%|    fval      [DV]            signal model function evaluations
%|
%|  version control
%|    1.1       2016-06-02      original
%|    1.2       2016-06-06      added magnitude signal options
%|    1.3       2016-06-08      compatibility with spgr/dess
%|    1.4       2016-06-27      new se-ir signal model
%|    1.5       2016-06-28      nu now expected as cell (previously struct)

% initialize with matrix format for convenient indexing
dim.Vloc = numel(x{1});                                                 % [D Vloc]
fval = NaN(dim.D, dim.Vloc);
row = 1;

% se-ir function calls
for s = 1:dim.S.ir
  fval(row,:) = IR_fun_v4(...
    x{1:3},...
    P.ir.tr(s), P.ir.ti(s), P.ir.te(s,1),...
    'inveff', x{4},...
    'kap', nu{1},...
    'wf', nu{2},...
    'flip_inv', P.ir.ainv(s),...
    'flip_ex', P.ir.aex(s),...
    'flip_ref', P.ir.aref(s),...
    'mag', bool.mag.ir);
  row = row+1;
end

% se function calls
for s = 1:dim.S.se
  fval(row,:) = SE_fun_v4(...
    x{1:3},...
    P.se.tr(s), P.se.te(s,1),...
    'kap', nu{1},...
    'wf', nu{2},...
    'flip_ex', P.se.aex(s),...
    'flip_ref', P.se.aref(s),...
    'mag', bool.mag.se);
  row = row+1;
end

% spgr function calls
for s = 1:dim.S.sp
  fval(row,:) = spgr_fun_v2(...
    x{1:3},...
    P.sp.aex(s), P.sp.tr(s), P.sp.te(s,1),...
    'kap', nu{1},...
    'dw', nu{2},...
    'R2p', nu{3},...
    'mag', bool.mag.sp);
  row = row+1;
end

% dess function calls
for s = 1:dim.S.de
  [fval(row,:), fval(row+dim.S.de,:)] = dess_fun_v2(...
    x{1:3},...
    P.de.aex(s), P.de.tr(s), P.de.te(s,1), P.de.te(s,2),...
    'kap', nu{1},...
    'dw', nu{2},...
    'R2p', nu{3},...
    'mag', bool.mag.de);
  row = row+1;
end

% reshape for output
fval = col(fval);                                                       % [DVloc]
end


  function [spX] = spmat(fullX)
%|function [spX] = spmat(fullX)
%|  
%|  efficiently construct sparse block-diagonal matrix
%|    useful when C >> A*B
%|    
%|  input
%|    fullX     [A B C]         full version
%|
%|  output
%|    spX       [AC BC sparse]  sparse block-diagonal version
%|
%|  version control
%|    1.1       2016-10-27      original

% array size
[A, B, C] = size(fullX);

% index/val instantiation
T = numel(fullX);
i = NaN(T,1);
j = NaN(T,1);
s = NaN(T,1);

% populate 
row = 1;
for a = 1:A
  for b = 1:B
    i(row : row+C-1) = (a : A : A*C)';                                  % [C]
    j(row : row+C-1) = (b : B : B*C)';                                  % [C]
    s(row : row+C-1) = fullX(a,b,:);                                    % [C]
    row = row + C;
  end
end

% sparse matrix
spX = sparse(i, j, s, A*C, B*C, T);                                     % [AC BC sparse]
  end
