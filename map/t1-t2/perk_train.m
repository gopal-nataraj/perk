  function train = perk_train(rff, dist, y, w, nu, P, dim, bool)
%|function train = perk_train(rff, dist, y, w, nu, P, dim, bool)
%|
%|  PERK training
%|
%|  inputs
%|    rff       [1x1 struct]    random fourier features object
%|     .snr     [1]               estimate of max sig for unity m0    
%|     .std     [1]               noise std dev in training data
%|     .len     [D+N]             kernel input length scales
%|     .c       [1]               global kernel length scale parameter 
%|     .H       [1]               embedding dimension       
%|     .K       [1]               number of training samples
%|    dist    	[1x1 struct]    sampling distribution object (ignored if x0 set!) 
%|     .x       {L cell}          latent parameter object
%|      .supp   [2]                 [lb ub] distribution support        
%|      .prior  {1}                 ('unif', 'logunif') distribution     
%|     .nu      {K cell}          known parameter object
%|      .supp   [2]                 [lb ub] distribution support     
%|    y         [DV]            coil-combined image data
%|    w         [D]             dataset weights
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
%|     .mag.*   false|true        using magnitude (spgr,dess) data
%|     .chat    false|true        verbosity  
%|     .reset   false|true        reset random number generator during training   
%|     .rfftst  false|true        show kernel approximation
%|     .nuclip  false|true        (perk) clip nu sampling distribtion   
%|
%|  outputs
%|    train     [1x1 struct]    perk training parameter object (if empty, will train)
%|     .mean.z  [H]               sample mean of feature maps
%|     .mean.x  [L]               sample mean of x
%|     .cov.zz  [H H]             sample auto-cov of feature maps
%|     .cov.xz  [L H]             sample cross-cov b/w x and feature maps
%|     .cov.xx  [L L]             sample auto-cov of latent parameters x
%|     .freq    [H D+N]           random 'frequency' vector
%|     .ph      [H]               random phase vector
%|  
%|  version control
%|    1.1       2017-06-06      adapted from mri_multicomp_map(...)
%|    1.2       2017-06-12      rff.snr now controls m0 distribution sampling
%|    1.3       2017-09-14      added perk nu distribution clipping

% optional: reset random number generator
if bool.reset
  % reset random number generator for reproducible results
  rng('default');
end

% sample x and set population means
xsamp = cell(dim.L,1);
for l = 1:dim.L
  % special case: set x{1}=m0 based on data scale
  if l==1 && isempty(dist.x{l}.supp)
    % assumption: prior on m0 ~unif(0,max(y)/rff.snr)
    tmp = div0(max(y),rff.snr);
    dist.x{l}.supp = [eps tmp].';
    xsamp{l} = random('unif', dist.x{l}.supp(1), dist.x{l}.supp(2), [rff.K, 1]);
    
    % population mean is (b-a)/2
    dist.x{l}.mean = mean(dist.x{l}.supp);
  else
    if l==1
      warn('Training w/ preset m0 range could cause KRR errors due to scale mismatch!');
    end
    switch dist.x{l}.prior
      case 'unif'
        xsamp{l} = random('unif', dist.x{l}.supp(1), dist.x{l}.supp(2), [rff.K, 1]);
        
        % population mean is (b-a)/2
        dist.x{l}.mean = mean(dist.x{l}.supp);
      case 'logunif'
        tmp = log(dist.x{l}.supp);
        tmp2 = random('unif', tmp(1), tmp(2), [rff.K, 1]);
        xsamp{l} = exp(tmp2);
        
        % population mean is (b-a)/(log(b)-log(a))
        dist.x{l}.mean = div0((dist.x{l}.supp(2)-dist.x{l}.supp(1)),tmp(2)-tmp(1));
      otherwise
        error('Unknown distribution prior option for latent parameter %u.', l);
    end
  end
end

% sample nu
nusamp = cell(dim.N,1);
for n = 1:dim.N
  if max(nu{n})-min(nu{n}) < eps
    % special case: delta function
    nusamp{n} = ones(rff.K,1) * mean(nu{n});
  else
    if bool.nuclip
      % fit using only samples within support
      tmp = masker(nu{n}, nu{n}>=dist.nu{n}.supp(1) & nu{n}<=dist.nu{n}.supp(2));

      % kernel density estimate, with finite support
      kde = fitdist(tmp, 'kernel',...
        'support', [dist.nu{n}.supp(1), dist.nu{n}.supp(2)]);
    else
      % kernel density estimate, with possibly infinite support
      kde = fitdist(nu{n}, 'kernel');
    end
    
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
q = [ysamp; tmp];                                                       % [D+N K]
rff.Q = size(q,1);

% random fourier features
% to approximate gaussian kernel N(0, Sigma),
%   1. construct rff.cov = inv((2*pi)^2*Sigma)
%   1. draw rff.freq from N(0, rff.cov)
%   2. draw rff.ph from unif(0,1)
tmp = rff.len * (2*pi*rff.c);
tmp = div0(1,tmp.^2);                                                   % [D+N D+N]
rff.cov = spdiags(tmp, 0, rff.Q, rff.Q);                                % [D+N D+N sparse]
train.freq = randn(rff.H, rff.Q) * chol(rff.cov);                       % [H D+N]
train.ph = rand(rff.H,1);                                               % [H]

% feature maps
z = rff_map(q, rff, w, train.freq, train.ph, dim);                      % [H K]
if bool.rfftst
  rff_map_test(q(:,1:100), rff, w, train.freq, train.ph, dim);
end

% sample means
train.mean.z = mean(z,2);                                               % [H]
tmp = transpose([xsamp{:}]);                                            % [L K]
train.mean.x = mean(tmp,2);                                             % [L]

% sample covariances
tmp = bsxfun(@minus, tmp, train.mean.x);
z = bsxfun(@minus, z, train.mean.z);
train.cov.zz = div0(z * z', rff.K);                                     % [H H]
train.cov.xz = div0(tmp * z', rff.K);                                   % [L H]
train.cov.xx = div0(tmp * tmp', rff.K);                                 % [L L]

% free memory
clear('z');
end