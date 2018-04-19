  function [x, t] = mri_m0t1t2inveff_map(y, P, varargin)
%|function [x, t] = mri_m0t1t2inveff_map(y, P, varargin)
%|
%|  regularized least-squares m0,t1,t2(,inveff) estimation 
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
%|    y         [1x1 struct]    coil-combined image data
%|     .ir      [(odims) S.ir 1]  spin-echo inversion-recovery
%|     .se      [(odims) S.se 1]  spin-echo
%|     .sp      [(odims) S.sp 1]  spoiled gradient-recalled echo
%|     .de      [(odims) S.de 2]  dual-echo steady state 
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
%|     .est     [(odims)]         over which to perform estimation      def: imdilate(mask.disp)
%|     .noise   [(odims)]         over which to estimate noise std dev  def: ~mask.est
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
%|    x0        [1x1 struct]    latent parameter initial estimates      def: from ML
%|     .m0      [(odims)]         spin density 
%|     .t1      [(odims)]         spin-lattice relaxation time                                ms
%|     .t2      [(odims)]         spin-spin relaxation time                                   ms
%|     .inveff  [(odims)]         inversion efficiency (S.ir>0 only)
%|    meth      [1x1 struct]    estimation methods                  
%|     .init    {1}               initialization                        def: 'krr'
%|     .iter    {1}               iterative local optimization          def: 'pgpm'  
%|    dist.x    [1x1 struct]    latent parameter sampling distribution object (ignored if x0 set!)  
%|                                2nd field: latent parameter (m0,t1,t2(,inveff))
%|     .*.supp  [2]               [lb ub] distribution support          def: see below
%|     .*.nsamp [1]               (vpm) number of dict samples          def: see below
%|     .*.prior {1}               ('unif', 'logunif') distribution      def: see below 
%|    dist.nu   [1x1 struct]    known parameter sampling distribution object
%|                                2nd field: known parameter (kap,b0,r2p)
%|     .*.supp  [2]               [lb ub] distribution support          def: see below
%|    kmean     [1x1 struct]    (vpm,pgpm) kmeans object to pool x0/nu
%|     .C       [1]               number of kmeans clusters             def: 10
%|     .opt     {1 2*nopt}        optional arguments to kmeans          def: see below
%|    rff       [1x1 struct]    (krr) random fourier features object
%|     .snr     [1]               estimate of max sig for unity m0      def: 0.1
%|     .std     [1]               noise std dev in training data        def: est from noise 
%|     .len     [D+N]             kernel input length scales            def: from data
%|     .c       [1]               global kernel length scale parameter  def: 2^0
%|     .H       [1]               embedding dimension                   def: 10^4
%|     .K       [1]               number of training samples            def: 10^6
%|    train     [1x1 struct]    (krr) training parameter object         def: from training
%|     .mean.z  [H]               sample mean of feature maps
%|     .mean.x  [L]               sample mean of x
%|     .cov.zz  [H H]             sample auto-cov of feature maps
%|     .cov.xz  [L H]             sample cross-cov b/w x and feature maps
%|     .cov.xx  [L L]             sample auto-cov of latent parameters x
%|     .freq    [H D+N]           random 'frequency' vector
%|     .ph      [H]               random phase vector
%|    inv       [1x1 struct]    matrix inversion reg strength object
%|     .krr     [1]               kernel ridge regression               def: 10^-8
%|     .lm      [1]               levenberg-marquardt                   def: 10^-8
%|    precon    [1x1 struct]    preconditioning parameter object        
%|     .update  [1]               number of iter to update precon       def: 5
%|     .hessreg [1]               hessian regularization parameter      def: 10^-10
%|    line      [1x1 struct]    backtracking line search object
%|     .step0   [1]               initial step size                     def: 1
%|     .prop    [1]               step size reduction factor            def: 0.5
%|    boxcon    [1x1 struct]    box constraints for iter opt   
%|                                1st field: latent parameter (m0,t1,t2(,inveff))
%|    reg       [1x1 struct]    edge-preserving (hyper3) regularizers (see Reg1.m)   
%|                                1st field: latent parameter (m0,t1,t2(,inveff))
%|     .*.pot   {1 npotarg}       potential function arguments          def: {'hyper3', delta}      
%|     .*.beta  [1]               strength                              def: scales with D
%|    stop      [1x1 struct]    stop criteria
%|     .iter    [1]               maximum number of iterations          def: 100
%|     .wghtx   [L]               provides weighting in weighted norm   def: from boxcon
%|     .tolx    [1]               compare vs wnorm(x-xprev)/wnorm(x)    def: 10^-7
%|    bool      [1x1 struct]    boolean variables
%|     .mag.*   false|true        using magnitude (ir,se,sp,de) data    def: true
%|     .chat    false|true        verbosity                             def: false
%|     .norm    false|true        normalize data for scale-invar reg    def: true
%|     .m0dess  false|true        (mom) use dess data for m0 est        def: false
%|     .reset   false|true        (krr) reset rng while sampling        def: true
%|     .rfftst  false|true        (krr) show kernel approximation       def: false
%|     .nuclip  false|true        (krr) clip nu sampling distribtion    def: true
%|     .reg     false|true        use regularization                    def: false
%|     .precon  false|true        use preconditioner (hessian approx)   def: true
%|     .disp    false|true        show image updates                    def: false
%|    disp      [1x1 struct]    display ranges for image updates       
%|     .m0      [2]               abs(m0) display range                 def: [0 2*median(abs(m0))]
%|     .t1      [2]               t1 display range                      def: [0 2000]         ms
%|     .t2      [2]               t2 display range                      def: [0 200]          ms
%|     .inveff  [2]               inveff display range (S.ir>0 only)    def: [0.5 1]      
%|
%|  outputs
%|    x         [1x1 struct]    object estimates
%|     .init.*  [(odims)]         initial estimates
%|     .iter.*  [(odims)]         iterative estimates, if computed
%|                                2nd field: latent parameter (m0,t1,t2(,inveff))
%|    t.*       [1x1 struct]    run times (init,iter)
%|
%|  copyright 2016, gopal nataraj, university of michigan
%|
%|  version control
%|    1.1       2016-06-01      skeleton code (mri_irse_se_m0t1t2_map.m)
%|    1.2       2016-06-05      ml estimator working
%|    1.3       2016-06-06      rls estimator working
%|    1.4       2016-06-07      separate estimation/display masks, omit zero-weight data
%|    2.1       2016-06-08      compatibility with spgr/dess 
%|    2.2       2016-06-13      display range options
%|    2.3       2016-06-14      warnings when nu.* maps unset
%|    2.4       2016-06-22      normalize data *before* ml est (previously after)
%|    3.1       2016-06-27      jointly estimating inveff also when profiles contain ir scans
%|    4.1       2016-06-28      x0 may be partially preset; vpm now clusters over nu
%|    4.2       2016-07-02      switched from proj gauss-newton to proj levenberg-marquardt
%|    4.3       2016-07-17      if nu and preset x0 are all uniform, set kmean.C to 1
%|    5.1       2017-06-06      added krr init option
%|    5.2       2017-06-12      rff.snr now controls m0 distribution sampling
%|    5.3       2017-06-17      added krr regularization strengh option
%|    5.4       2017-07-10      switched output format from e.g. x.m0.init to x.init.m0
%|    5.5       2017-09-14      added krr nu distribution clipping
%|    6.1       2017-12-19      added method-of-moments init option for spgr/dess
%|    6.2       2018-01-12      pgpm now loops over tissue clusters for separate line searches
%|    6.3       2018-02-06      if rff.std set and bool.norm=1, now rescales rff.std

% object dimensions
tmp = isfield(y, {'ir','se','sp','de'});
if sum(tmp)==0
  error('Detected no data?!');
end
for i = 1:length(tmp)
  if tmp(i)
    switch i
      case 1
        tmp2 = size(y.ir);
        dim.odims = tmp2(1:end-1);
      case 2
        tmp2 = size(y.se);
        dim.odims = tmp2(1:end-1);
      case 3
        tmp2 = size(y.sp);
        dim.odims = tmp2(1:end-1);
      case 4
        tmp2 = size(y.de);
        dim.odims = tmp2(1:end-2);
    end
  end
end
  
% create empty arrays for missing fields
for i = 1:length(tmp)
  if ~tmp(i)
    switch i
      case 1, y.ir = zeros([dim.odims 0 1]);
      case 2, y.se = zeros([dim.odims 0 1]);
      case 3, y.sp = zeros([dim.odims 0 1]);
      case 4, y.de = zeros([dim.odims 0 2]);
    end
  end
end

% number of scans
dim.S.ir = size(y.ir, length(dim.odims)+1);
dim.S.se = size(y.se, length(dim.odims)+1);                                            
dim.S.sp = size(y.sp, length(dim.odims)+1);
dim.S.de = size(y.de, length(dim.odims)+1);

% number of echoes
dim.E.ir = size(y.ir, length(dim.odims)+2);
dim.E.se = size(y.se, length(dim.odims)+2);                                            
dim.E.sp = size(y.sp, length(dim.odims)+2);
dim.E.de = size(y.de, length(dim.odims)+2);

% number of latent parameters
if dim.S.ir>0
  dim.L = 4;
else
  dim.L = 3;
end

% default values
arg.mask.est = [];
arg.mask.disp = [];
arg.mask.noise = [];

arg.nu.kap = [];
arg.nu.b0 = [];
arg.nu.r2p = [];

arg.wght.ir = ones(dim.S.ir, dim.E.ir);
arg.wght.se = ones(dim.S.se, dim.E.se);
arg.wght.sp = ones(dim.S.sp, dim.E.sp);
arg.wght.de = ones(dim.S.de, dim.E.de);

arg.thresh = 0.05;

arg.x0.m0 = [];
arg.x0.t1 = [];
arg.x0.t2 = [];
if dim.L>3
  arg.x0.inveff = [];
end

arg.meth.init = 'krr';
arg.meth.iter = 'pgpm';

arg.dist.x.m0.supp = [];
arg.dist.x.m0.nsamp = 1;
arg.dist.x.m0.prior = 'unif';
arg.dist.x.t1.supp = [100 3000].';
arg.dist.x.t1.nsamp = 300;
arg.dist.x.t1.prior = 'logunif';
arg.dist.x.t2.supp = [1 700].';
arg.dist.x.t2.nsamp = 300;
arg.dist.x.t2.prior = 'logunif';
if dim.L>3
  arg.dist.x.inveff.supp = [0.51 1].';
  arg.dist.x.inveff.nsamp = 50;
  arg.dist.x.inveff.prior = 'unif';
end
arg.dist.nu.kap.supp = [0.5 2].';
arg.dist.nu.b0.supp = [0 0].';
arg.dist.nu.r2p.supp = [0 0].';

arg.kmean.C = 10;
arg.kmean.opt = {...
  'EmptyAction', 'singleton',...
  'MaxIter', 1000};

arg.rff.snr = 0.1;
arg.rff.std = [];
arg.rff.len = [];
arg.rff.c = 2^0;
arg.rff.H = 10^4;
arg.rff.K = 10^6;

arg.train = [];

arg.inv.krr = 10^-8;
arg.inv.lm = 10^-8;

arg.precon.update = 5;
arg.precon.hessreg = 10^-10;

arg.line.step0 = 1;
arg.line.prop = 0.5;

arg.boxcon.m0 = [-Inf Inf].';
arg.boxcon.t1 = [100 3000].';
arg.boxcon.t2 = [1 700].';
if dim.L>3
  arg.boxcon.inveff = [0.5 1].';
end

arg.reg.m0.pot  = {'hyper3', 2^-2};
arg.reg.m0.beta = [];
arg.reg.t1.pot  = {'hyper3', 2^5};
arg.reg.t1.beta = [];
arg.reg.t2.pot  = {'hyper3', 2^2};
arg.reg.t2.beta = [];
if dim.L>3
  arg.reg.inveff.pot = {'hyper3', 2^-3};
  arg.reg.inveff.beta = [];
end

arg.stop.iter = 100;
arg.stop.wghtx = [];
arg.stop.tolx = 10^-7;

arg.bool.mag.ir = 1;
arg.bool.mag.se = 1;
arg.bool.mag.sp = 1;
arg.bool.mag.de = 1;
arg.bool.chat   = 0;
arg.bool.norm   = 1;
arg.bool.m0dess = 0;
arg.bool.reset  = 1;
arg.bool.rfftst = 0;
arg.bool.nuclip = 1;
arg.bool.reg    = 0;
arg.bool.precon = 1;
arg.bool.disp   = 0;

arg.disp.m0 = [];
arg.disp.t1 = [0 2000];
arg.disp.t2 = [0 200];
if dim.L>3
  arg.disp.inveff = [0.5 1];
end

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

% noise mask
if isempty(arg.mask.noise)
  arg.mask.noise = ~arg.mask.est;
end

% set V to number of voxels in mask.est
dim.V = numel(arg.mask.est(arg.mask.est));

% check for ir complex data
if ~isreal(y.ir) && arg.bool.mag.ir
  warn('\nDetected complex IR data but bool.mag.ir set true!? Reverting to default false.');
  arg.bool.mag.ir = false;
elseif ~isempty(y.ir) && isreal(y.ir) && ~arg.bool.mag.ir
  warn('\nDetected pure-real IR data; recommend setting bool.mag.ir true for speed.');
end

% check for se complex data
if ~isreal(y.se) && arg.bool.mag.se
  warn('\nDetected complex SE data but bool.mag.se set true!? Reverting to default false.');
  arg.bool.mag.se = false;
elseif ~isempty(y.se) && isreal(y.se) && ~arg.bool.mag.se
  warn('\nDetected pure-real SE data; recommend setting bool.mag.se true for speed.');
end

% check for spgr complex data
if ~isreal(y.sp) && arg.bool.mag.sp
  warn('\nDetected complex SPGR data but bool.mag.sp set true!? Reverting to default false.');
  arg.bool.mag.sp = false;
elseif ~isempty(y.sp) && isreal(y.sp) && ~arg.bool.mag.sp
  warn('\nDetected pure-real SPGR data; recommend setting bool.mag.sp true for speed.');
end

% check for dess complex data
if ~isreal(y.de) && arg.bool.mag.de
  warn('\nDetected complex DESS data but bool.mag.de set true!? Reverting to default false.');
  arg.bool.mag.de = false;
elseif ~isempty(y.de) && isreal(y.de) && ~arg.bool.mag.de
  warn('\nDetected pure-real DESS data; recommend setting bool.mag.de true for speed.');
end

% omit ir scans with zero weight
keep.ir = arg.wght.ir(:,1)>0;
if ~all(keep.ir)
  y.ir = y.ir(:,:,keep.ir,:);
  P.ir.tr = P.ir.tr(keep.ir);
  P.ir.ti = P.ir.ti(keep.ir);
  P.ir.te = P.ir.te(keep.ir, dim.E.ir);
  P.ir.ainv = P.ir.ainv(keep.ir);
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
    arg.dist.x = rmfield(arg.dist.x, 'inveff');
    arg.boxcon = rmfield(arg.boxcon, 'inveff');
    arg.reg = rmfield(arg.reg, 'inveff');
    arg.disp = rmfield(arg.disp, 'inveff');
    if arg.bool.chat
      warn('\nNo nonzero-weight IR scans remain; omitting inveff estimation.');
    end
  end
end

% omit se scans with zero weight
keep.se = arg.wght.se(:,1)>0;
if ~all(keep.se)
  y.se = y.se(:,:,keep.se,:);
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
  y.sp = y.sp(:,:,keep.sp,:);
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
  y.de = y.de(:,:,keep.de,:);
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

% vectorize data
y = cat(length(dim.odims)+1,...
  reshape(y.ir, [dim.odims dim.Dir]),...
  reshape(y.se, [dim.odims dim.Dse]),...
  reshape(y.sp, [dim.odims dim.Dsp]),...
  reshape(y.de, [dim.odims dim.Dde]));                                  % [(odims) D]
if strcmp(arg.meth.init, 'krr')
  n = masker(y, arg.mask.noise);                                        % [V_bg D]
  n = col(transpose(n));                                                % [DV_bg]
end
y = masker(y, arg.mask.est);                                            % [V D]
y = col(transpose(y));                                                  % [DV]

% trick:  normalize data by median of sum of magnitude non-background values (over datasets)
%         effective regularization strength is then scale-invariant
if arg.bool.norm
  tmp = sum(reshape(abs(y), [dim.D dim.V]), 1);
  scale = median(tmp(tmp > arg.thresh * max(col(tmp))));
  y = div0(y,scale);
  if ~isempty(arg.x0.m0)
    arg.x0.m0 = div0(arg.x0.m0,scale);
  elseif strcmp(arg.meth.init, 'krr')
    n = div0(n,scale);
  end
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
arg.nu.kap  = arg.nu.kap(arg.mask.est);                                 % [V]
arg.nu.b0   = arg.nu.b0(arg.mask.est);                                  % [V]  
arg.nu.r2p  = arg.nu.r2p(arg.mask.est);                                 % [V]

% vectorize data weights
w = [...
  col(arg.wght.ir);...
  col(arg.wght.se);...
  col(arg.wght.sp);...
  col(arg.wght.de)];                                                    % [D]

% warn user if attempt is made to constrain m0 estimation
if ~all(isinf(arg.boxcon.m0))
  warn('Unconstrained m0 estimation: ignoring arg.boxcon.m0.* settings.');
  arg.boxcon.m0 = [-Inf Inf].';
end

% to enumerate conveniently, switch from struct to cell
arg.nu = struct2cell(arg.nu);                                           % {N} cell of [V] arrays
arg.dist.x = struct2cell(arg.dist.x);                                   % {L} cell of [1] structs
arg.dist.nu = struct2cell(arg.dist.nu);                                 % {N} cell of [1] structs
arg.x0 = struct2cell(arg.x0);                                           % {L} cell of [(odims)] arrays
arg.boxcon = struct2cell(arg.boxcon);                                   % {L} cell of [2] arrays
arg.disp = struct2cell(arg.disp);                                       % {L} cell of [1 2] arrays
dim.N = length(arg.nu);

% if any x0 fields left unspecified, initialization necessary
tmp = 0;
for l = 1:dim.L
  tmp = tmp + ~isempty(arg.x0{l});
end
if tmp<dim.L
  % choose initialization method
  switch arg.meth.init
    case 'vpm'
      % warn user if attempt is made to constrain m0 ml estimation
      if ~isempty(arg.dist.x{1}.supp) || arg.dist.x{1}.nsamp~=1
        warn('\nVarPro ML estimation fixes unity m0: ignoring arg.dist.x.m0.* settings.');
        arg.dist.x{1}.supp  = [];
        arg.dist.x{1}.nsamp   = 1;
        arg.dist.x{1}.prior   = 'unif';
      end
      arg.dist.x{1}.supp = [1 1]';
      
      % use only one k-means cluster if all nu and preset x0 fields uniform
      tmp = true;
      for n = 1:dim.N
        if min(arg.nu{n})~=max(arg.nu{n})
          tmp = false;
          break;
        end
      end
      if tmp
        for l = 1:dim.L
          if ~isempty(arg.x0{l})
            tmp2 = minmax(masker(arg.x0{l}, arg.mask.est));
            if min(tmp2)~=max(tmp2)
              tmp = false;
              break;
            end
          end
        end
      end
      if tmp
        fprintf('Detected uniform known maps: setting kmean.C to 1.\n');
        arg.kmean.C = 1;
      end

      % max-likelihood estimation via variable projection method
      if arg.bool.chat
        fprintf('\n=============================================================');
        fprintf('\nInitialization via VarPro...');
        fprintf('\n=============================================================\n');
      end
      [arg.x0, t.init] = varpro(...
        y, arg.x0, arg.nu, P, w,...
        arg.kmean, arg.dist.x, dim, arg.bool, arg.mask.est);            % {L} cell w/ [V] arrays
    case 'mom'
      % check for at least 2 spgr and 1 dess scan
      if dim.Dsp<2 && isempty(arg.x0{2})
        error('need at least 2 spgr scans for method-of-moments t1 estimation!');
      end
      if dim.Dde<1 && isempty(arg.x0{3})
        error('need at least 1 dess scan for method-of-moments t2 estimation!');
      end
      
      % omit ir/se scans and warn accordingly 
      tmp = reshape(y, [dim.D dim.V]);
      if dim.Dir>0
        tmp2 = strcat(...
          'Method-of-moments t1/t2 estimation from ir scans not yet implemented.\n',...
          'Omitting ir scans from method-of-moments initialization...');
        warn(tmp2);
        tmp = tmp(dim.Dir+1:end,:);
      end
      if dim.Dse>0
        tmp2 = strcat(...
          'Method-of-moments t1/t2 estimation from se scans not yet implemented.\n',...
          'Omitting se scans from method-of-moments initialization...');
        warn(tmp2);
        tmp = tmp(dim.Dse+1:end,:);
      end
      
      % method-of-moments estimation via algebraic manipulations
      if arg.bool.chat
        fprintf('\n=============================================================');
        fprintf('\nInitialization via Method-of-Moments Estimation...');
        fprintf('\n=============================================================\n');
      end
      [arg.x0, t.init] = mom(...
        tmp, arg.x0, arg.nu, P, w, dim, arg.bool, arg.mask.est);        % {L} cell w/ [V] arrays
    case 'krr'
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
      
      % set training data noise std dev
      if isempty(arg.rff.std)
        if isempty(n)
          if arg.bool.norm
            arg.rff.std = div0(3.8607e-4,scale);
          else
            arg.rff.std = 3.8607e-4;
          end
          warn('No noise mask given: rff.std set automatically to %0.8f.', arg.rff.std);
        else
          if isreal(n)
            % n is assumed ~rayleigh(\sigma)
            % rff.std is a (slightly biased) estimate of \sigma
            tmp = sum(n.^2);
            tmp = div0(tmp,2*length(n));
            arg.rff.std = sqrt(tmp);
          else
            error('todo: estimate rff.std for complex data');
          end
        end  
      elseif arg.bool.norm
        arg.rff.std = div0(arg.rff.std,scale);
        warn('Detected rff.std and bool.norm=1: normalized rff.std by %0.8f.', scale);
      end
      
      % set kernel length scales 
      if isempty(arg.rff.len)
        tmp = reshape(y, [dim.D dim.V]);
        tmp = [tmp; transpose([arg.nu{:}])];
        arg.rff.len = mean(tmp,2);
        arg.rff.len = max(arg.rff.len,eps);
      end
      
      % kernel regression
      if arg.bool.chat
        fprintf('=============================================================');
        fprintf('\nInitialization via Kernel Regression...');
        fprintf('\n=============================================================\n');
      end
      [arg.x0, t.init] = krr(...
        y, arg.x0, arg.nu, P, w,...
        arg.rff, arg.train, arg.dist, arg.inv.krr,...
        dim, arg.bool, arg.mask.est);                                 % {L} cell w/ [V] blocks
    otherwise
      error('Unknown initialization method requested!');
  end
  if arg.bool.chat  
    fprintf('========================================================================');
    fprintf('\n                                    ...done in %0.3f seconds.', t.init);
    fprintf('\n========================================================================\n');
  end
else
  for l = 1:dim.L
    arg.x0{l} = masker(arg.x0{l}, arg.mask.est);
  end
end

if arg.bool.reg
  % if not preset, set regularizer strengths proportional to D
  if isempty(arg.reg.m0.beta)
    arg.reg.m0.beta = dim.D * 2^-26;
  end
  if isempty(arg.reg.t1.beta)
    arg.reg.t1.beta = dim.D * 2^-21;
  end
  if isempty(arg.reg.t2.beta)
    arg.reg.t2.beta = dim.D * 2^-23;
  end
  if dim.L>3 && isempty(arg.reg.inveff.beta)
    arg.reg.inveff.beta = dim.D * 2^-26;
  end
  
  % to enumerate conveniently, switch from struct to cell
  arg.reg = struct2cell(arg.reg);                                         % {L}

  % regularizer objects
  for l = 1:dim.L
    arg.reg{l}.R = Reg1(arg.mask.est,...
      'pot_arg', arg.reg{l}.pot,...
      'beta', arg.reg{l}.beta,...
      'type_penal', 'mat');
  end
end
  
% set iterative stopping criterion weights
arg.stop.wghtx = NaN(dim.L,1);
for l=1:dim.L
  if sum(isinf(arg.boxcon{l}))>0
    tmp = arg.x0{l};
    tmp = median(tmp(tmp > arg.thresh * max(col(tmp))));
  else
    tmp = mean(arg.boxcon{l});
  end
  arg.stop.wghtx(l) = div0(1,tmp);
end

% choose iterative method
switch arg.meth.iter
  case 'plm'
    if arg.bool.chat
      fprintf('\n=============================================================');
      fprintf('\nIterative estimation via projected levenberg-marquardt...');
      fprintf('\n=============================================================\n');
    end
    [tmp, t.iter] = plm(...
      y, arg.x0, arg.nu, P, w, arg.inv.lm, arg.boxcon,...
      arg.reg, arg.stop, dim, arg.bool, arg.mask.est, arg.disp);        % {L} cell w/ [V] blocks
  case 'pgpm'
    if arg.bool.chat
      fprintf('\n=============================================================');
      fprintf('\nIterative estimation via gradient projection method...');
      fprintf('\n=============================================================\n');
    end
    [tmp, t.iter] = pgpm(...
      y, arg.x0, arg.nu, P, w, arg.kmean, arg.precon,...
      arg.line, arg.boxcon, arg.reg, arg.stop, dim, arg.bool);          % {L} cell w/ [V] blocks
  otherwise
    error('Unknown iterative method!');
end
if arg.bool.chat
  fprintf('==================================================');
  fprintf('\n                        ...done in %0.3f seconds.', t.iter);
  fprintf('\n==================================================\n\n');
end

% trick: rescale m0 for output
if arg.bool.norm
  arg.x0{1} = arg.x0{1} * scale;
  tmp{1} = tmp{1} * scale;
end

% embed for output
x.init.m0 = embed(arg.x0{1}, arg.mask.est);                             % [(odims)]
x.init.t1 = embed(arg.x0{2}, arg.mask.est);                             % [(odims)]
x.init.t2 = embed(arg.x0{3}, arg.mask.est);                             % [(odims)]
if dim.L>3  
  x.init.inveff = embed(arg.x0{4}, arg.mask.est);                       % [(odims)]
end

x.iter.m0 = embed(tmp{1}, arg.mask.est);                                % [(odims)]
x.iter.t1 = embed(tmp{2}, arg.mask.est);                                % [(odims)]
x.iter.t2 = embed(tmp{3}, arg.mask.est);                                % [(odims)]
if dim.L>3  
  x.iter.inveff = embed(tmp{4}, arg.mask.est);                          % [(odims)]
end

% use mask.disp for display
x.init.m0(~arg.mask.disp) = 0;                                          % [(odims)]
x.init.t1(~arg.mask.disp) = 0;                                          % [(odims)]
x.init.t2(~arg.mask.disp) = 0;                                          % [(odims)]
if dim.L>3
  x.init.inveff(~arg.mask.disp) = 0;                                    % [(odims)]
end

x.iter.m0(~arg.mask.disp) = 0;                                          % [(odims)]
x.iter.t1(~arg.mask.disp) = 0;                                          % [(odims)]
x.iter.t2(~arg.mask.disp) = 0;                                          % [(odims)]
if dim.L>3
  x.iter.inveff(~arg.mask.disp) = 0;                                    % [(odims)]
end
  end
