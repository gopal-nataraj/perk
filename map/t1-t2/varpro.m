  function [x_ml, t_ml] = varpro(y, x0, nu, P, w, kmean, dist, dim, bool, mask)
%|function [x_ml, t_ml] = varpro(y, x0, nu, P, w, kmean, dist, dim, bool, mask)
%|
%|  maximum-likelihood latent parameter estimation via variable projection method
%|    preclusters known parameter maps and computes separate dictionaries for each cluster
%|
%|  inputs
%|    y         [DV]            coil-combined image data
%|    x0        {L cell}        latent parameter initial estimates
%|                                any subset of {m0,t1,t2,inveff}  
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
%|    w         [D]             dataset weights 
%|    kmean     [1x1 struct]    kmeans object to pool x0/nu 
%|    dist      [1x1 struct]    sampling distribution object (ignored if x0 set!)
%|                                1st field: latent parameter (m0,t1,t2(,inveff))
%|     .*.supp  [2]               [lb ub] distribution endpoints        
%|     .*.nsamp [1]               number of dict samples      
%|     .*.prior {1}               ('unif', 'logunif') distribution       
%|    dim       [1x1 struct]    object containing dimension info
%|    bool      [1x1 struct]    boolean variables        
%|     .mag.*   false|true        using magnitude (spgr,dess) data
%|     .chat    false|true        verbosity  
%|    mask      [(odims)]       defines region over which to estimate
%|
%|  outputs
%|    x_ml      {L cell}        latent object parameters
%|    t_ml      [1]             run time
%|  
%|  copyright 2016, gopal nataraj, university of michigan
%|
%|  version control
%|    1.1       2016-06-01      original
%|    1.2       2016-06-08      modified how dataset weights are passed
%|    1.3       2016-06-27      now handles inveff estimation when ir data passed
%|    2.1       2016-06-28      cluster known parameters, compute separate dictionary per cluster
%|    2.2       2017-06-06      for large K, now computes dictionary elements on-the-fly

% reshape y as required here
y = reshape(y, [dim.D dim.V]);

% transfer prespecified x0 maps as known parameters
x_ml = cell(dim.L,1);
x.known = false(dim.L,1);
for l=1:dim.L
  if ~isempty(x0{l})
    if bool.chat
      fprintf('Detected prespecified l=%u latent parameter; omitting from ML estimation...\n', l);
    end
    x_ml{l} = masker(x0{l}, mask);
    x.known(l) = true;
  else
    x_ml{l} = NaN(dim.V,1);
  end
end

% cluster nu and known x maps
[cl.idx, tmp] = ...
  kmeans([nu{:}, x_ml{x.known}], kmean.C, kmean.opt{:});                % [V], [C N+L_known]
cl.cent.nu = tmp(:,1:dim.N);                                            % [C N]
cl.cent.x0 = tmp(:,dim.N+1:end);                                        % [C L_known]

% weighting matrices
Wd = spdiags(w, 0, dim.D, dim.D);                                       % [D D sparse]
W = spdiags(col(repmat(w, [1 dim.V])), 0, dim.D*dim.V, dim.D*dim.V);    % [DV DV sparse]

% run varpro cluster-by-cluster
tic;
for c = 1:kmean.C
  if bool.chat
    fprintf('VarPro: cluster %2u of %2u...', c, kmean.C);
  end
  
  % create dictionary sample locations
  % compute number of dictionary elements, K
  x.samp = cell(1,dim.L);
  x.nsamp = NaN(1,dim.L);
  dim.K = 1;
  for l = 1:dim.L
    % if x{l} known, use cth cluster value
    if x.known(l)
      x.samp{l} = cl.cent.x0(c, sum(x.known(1:l)));                     % [1]
      x.nsamp(l) = 1;
    % otherwise, sample a spread of values
    else
      switch dist{l}.prior
        case 'unif'
          x.samp{l} = col(linspace(dist{l}.supp(1), dist{l}.supp(2), dist{l}.nsamp));
        case 'logunif'
          x.samp{l} = col(logspace(...
            log10(dist{l}.supp(1)), log10(dist{l}.supp(2)), dist{l}.nsamp));
        otherwise
          error('Unknown distribution prior option for latent parameter %u.', l);
      end
      dim.K = dim.K * dist{l}.nsamp;
      x.nsamp(l) = dist{l}.nsamp;
    end
  end

  % for small-medium K, compute and store dictionary signals via vector ops
  if dim.K <= 10^7
    % store sample locations
    x.grid = cell(1,dim.L);                                             % {1 L}, w/ lth element:
    [x.grid{:}] = ndgrid(x.samp{:});                                    % [x.nsamp(1) ... x.nsamp(L)]
    for l = 1:dim.L
      x.grid{l} = col(x.grid{l});                                       % [dim.K]
    end
    x.grid = col(x.grid);                                               % {L}

    % construct dictionary over all sample locations
    tmp = cell(dim.N,1);
    for n = 1:dim.N
      tmp{n} = ones(dim.K,1) * cl.cent.nu(c,n);
    end
    D = reshape(f(x.grid, tmp, P, dim, bool), [dim.D dim.K]);           % [D K]
    if bool.chat
      fprintf('computed dict of K=%u elements...', dim.K);
    end

    % variable-projection method
    %   effectively skips parameters (m0, l=1) linearly related to data
    maxProd = zeros(sum(cl.idx==c),1);
    bestIdx = ones(sum(cl.idx==c),1);
    for k = 1:dim.K
      % compute kth inner product
      tmp = abs(D(:,k)' * Wd * D(:,k));                                 % [1 1]
      tmp = D(:,k)' * Wd * y(:,cl.idx == c) / sqrt(tmp);                % [1 Vc]
      kthProd = transpose(abs(tmp).^2);                                 % [Vc]

      % update indices whose kth inner product is larger than max so far
      update = kthProd > maxProd;                                       % [Vc]
      maxProd(update) = kthProd(update);                                % [Vc]
      bestIdx(update) = k;
    end
    
  % for large K, compute dictionary elements on the fly
  else
    if bool.chat
      tmp = '\ndict of K=%u atoms is too large to store; will compute on the fly...\n';
      fprintf(tmp, dim.K);
    end
    
    % variable projection method for cth cluster
    %   effectively skips parameter (m0, l=1) linearly related to data
    maxProd = zeros(sum(cl.idx==c),1);
    bestIdx = ones(sum(cl.idx==c),1);
    idxk = cell(1,dim.L);
    x.k = cell(dim.L,1);
    for k = 1:dim.K
      if bool.chat && rem(k,100)==0
        fprintf('done with %u of %u dict elements...\n', k, dim.K);
      end
      
      % extract kth x combination
      [idxk{:}] = ind2sub(x.nsamp, k);                                  % [1 dim.L]
      for l = 1:dim.L
        x.k{l} = x.samp{l}(idxk{l});
      end

      % evaluate kth dictionary element
      dk = f(x.k, num2cell(col(cl.cent.nu(c,:))), P, dim, bool);        % [D]

      % compute kth inner product
      tmp = abs(dk' * Wd * dk);                                         % [1]
      tmp = dk' * Wd * y(:,cl.idx == c) / sqrt(tmp);                    % [1 Vc]
      kProd = transpose(abs(tmp).^2);                                   % [Vc]

      % update indices whose kth inner product is larger than max so far
      update = kProd > maxProd;                                         % [Vc]
      maxProd(update) = kProd(update);                                  % [Vc]
      bestIdx(update) = k;                                          
    end
  end
  
  % extract indices for max-likelihood maps
  % 	x{1} invalid here (all ones)
  %   spurious value required to call rgradx_f(...) 
  x.idx = cell(1,dim.L);
  [x.idx{:}] = ind2sub(x.nsamp, bestIdx);
  for l = 1:dim.L
    if ~x.known(l)
      x_ml{l}(cl.idx==c) = x.samp{l}(x.idx{l});                         % [Vc]
    end
  end
  if bool.chat
    fprintf('mapped %6u of %6u voxels.\n', sum(cl.idx<=c), dim.V);
  end
end
t_ml = toc;

% least-squares estimates (m0, l=1) for remaining parameters
% note we can now account for spatial variation of known parameters 
rgradx_f_ml = rgradx_f(x_ml, nu, P, dim, bool);                         % [DV LV sparse]
Am = rgradx_f_ml(:,1:dim.L:end);                                        % [DV V sparse]
num = Am'*W*col(y);                                                     % [V]
den = (Am'*W*Am)*ones(dim.V,1);                                         % [V]
x_ml{1} = div_safe(num, den);                                           % [V]
  end