  function [x, t] = pgpm(y, x0, nu, P, w, kmean, precon, line, boxcon, reg, stop, dim, bool)
%|function [x, t] = pgpm(y, x0, nu, P, w, kmean, precon, line, boxcon, reg, stop, dim, bool)
%|
%|  iterative latent parameter estimation via gradient projection method
%|    uses diagonal preconditioner to improve convergence
%|    uses backtracking line search to ensure monotone descent in cost
%|
%|  inputs
%|    y         [DV]            coil-combined image data
%|    x0        {L cell}        init est of latent parameters (cells size [V])
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
%|    precon    [1x1 struct]    preconditioning parameter object
%|     .update  [1]               number of iter to update precon
%|     .hessreg [1]               hessian regularization parameter  
%|    line      [1x1 struct]    backtracking line search object
%|     .step0   [1]               initial step size
%|     .prop    [1]               step size reduction factor
%|    boxcon    [1x1 struct]    [lb ub] box constraints
%|                                1st field: latent parameter (m0,t1,t2,(inveff))
%|    reg       [1x1 struct]    edge-preserving (hyper3) regularizers (see Reg1.m)   
%|                                1st field: latent parameter (m0,t1,t2,(inveff))
%|     .*.pot   {1 npotarg}       potential function arguments                  
%|     .*.beta  [1]               strength                              
%|    stop      [1x1 struct]    stop criteria
%|     .iter    [1]               maximum number of iterations          
%|     .wghtx   [L]               provides weighting in weighted norm   
%|     .tolx    [1]               compare vs wnorm(x-xprev)/wnorm(x)    
%|    dim       [1x1 struct]    object containing dimension info
%|    bool      [1x1 struct]    boolean variables
%|     .mag.*   false|true        using magnitude (ir,se,spgr,dess) data
%|     .chat    false|true        verbosity                             
%|     .norm    false|true        normalize data for scale-invar reg
%|     .reg     false|true        use regularization                   
%|     .precon  false|true        use preconditioner            
%|
%|  outputs  
%|    x         {L cell}        latent object parameter estimates
%|    t         [1]             run time
%|
%|  version control
%|    1.1       2017-06-06      adapted from mri_multicomp_map(...)
%|    1.2       2017-12-21      now initialization is also displayed
%|    2.1       2018-01-12      now loops over tissue clusters for separate line searches
%|    2.2       2018-01-14      now allow multiple levels of chat

% improve initialization by projecting onto box constraint
for l=1:dim.L
  x0{l} = max(x0{l}, boxcon{l}(1));
  x0{l} = min(x0{l}, boxcon{l}(2));
end

% reshape y as required here
y = reshape(y, [dim.D dim.V]);

% if regularizing, omit clustering
if bool.reg
  warn('multi-cluster regularization not implemented; reverting to single cluster.\n');
  kmean.C = 1;
end

% cluster x0 and nu maps
x = x0;
if stop.iter==0
  t = 0;
  return;
elseif kmean.C>=dim.V
  kmean.C = dim.V;
  cl.idx = ones(dim.V,1);
elseif kmean.C==1
  cl.idx = ones(dim.V,1);
else
  cl.idx = kmeans([x0{:} nu{:}], kmean.C, kmean.opt{:});                % [V], [C L+N]
end

% run pgpm cluster-by-cluster
tic;
for c = 1:kmean.C
  if bool.chat>=1
    fprintf('\npgpm: starting iterations on cluster %2u of %2u...', c, kmean.C);
  end
  
  % cluster-specific variables
  yc = col(y(:,cl.idx==c));                                             % [DVc]
  xc = cell(dim.L,1);
  for l = 1:dim.L
    xc{l} = x{l}(cl.idx==c);                                            % [Vc]
  end
  nuc = cell(dim.N,1);
  for n = 1:dim.N
    nuc{n} = nu{n}(cl.idx==c);                                          % [Vc]
  end
  dimc = dim;
  dimc.V = length(nuc{1});
  tmp = col(repmat(w, [1 dimc.V]));
  Wc = spdiags(tmp, 0, dim.D*dimc.V, dim.D*dimc.V);                     % [DVc DVc sparse] 
  flag = 0;
  
  for i=1:stop.iter
    % store previous cost
    psiPrev = Psi(yc, xc, nuc, P, Wc, dimc, reg, bool);
    if bool.chat>=2
      fprintf('\n  cost is %0.6e at iter %u of %u',...
        psiPrev, i, stop.iter);
    end
    
    % descent direction
    % for initial iterations, update preconditioner
    if i<precon.update
      [d, D] = pgpm_dir(...
        yc, xc, nuc, P, Wc, dimc, reg, precon.hessreg, bool);
    else
      d = pgpm_dir(...
        yc, xc, nuc, P, Wc, dimc, reg, precon.hessreg, bool, 'D', D);
    end
    
    % backtracking line search
    step = line.step0;
    while true
      % take step along descent direction
      xcv = col(transpose([xc{:}]));                                    % [LVc]
      xcvCand = xcv + (d*step);                                         % [LVc]
      xcCand = num2cell(reshape(xcvCand, [dimc.L, dimc.V]), 2);         % {L} cell w/ [1 Vc] blocks
      
      % project onto box constraints (omit m0)
      for l = 1:dim.L
        if l==1
          xcCand{l} = col(xcCand{l});                                   % m0 now w/ [Vc] block
        else
          xcCand{l} = bsxfun(@max, col(xcCand{l}), boxcon{l}(1));       % {L} cell w/ [Vc] blocks
          xcCand{l} = bsxfun(@min, col(xcCand{l}), boxcon{l}(2));       % {L} cell w/ [Vc] blocks
        end
        xcvCand(l:dim.L:end) = xcCand{l};
      end
      
      % if candidate iterate descends cost, exit to update iterate
      psiCand = Psi(yc, xcCand, nuc, P, Wc, dimc, reg, bool);
      if psiCand<psiPrev
        break;
      % else reduce step size and try to reduce cost again
      else
        if bool.chat>=3
          fprintf('\n    candidate cost %0.6e; step reduced from %0.2e to %0.2e',...
            psiCand, step, step*line.prop);
        end
        step = step*line.prop;
      end
      
      % if step gets too small, break twice to next cluster
      num = wnorm2(xcvCand, stop.wghtx, dimc);
      den = wnorm2(d, stop.wghtx, dimc);
      stepTol = stop.tolx * div0(num,den);
      stepTol = max(stepTol, eps);
      if step<stepTol
        if bool.chat>=3
          fprintf('\n    return early: step < stepTol = %0.6e', stepTol);
        end
        flag = 1;
        break;
      end
    end
    if flag
      break;
    end
    
    % update iterate
    xc = xcCand;
    
    % if minimal change in iterate, break to next cluster
    num = wnorm2(xcvCand-xcv, stop.wghtx, dimc);
    den = wnorm2(xcv, stop.wghtx, dimc);
    tmp = div0(num,den);
    if tmp<stop.tolx
      if bool.chat>=2
        fprintf('\n  return early: wnorm(x-xprev)/wnorm(x) < xTol=%0.3e', stop.tolx); 
      end
      break;
    end
  end
  if bool.chat>=2
    fprintf('\npgpm: cluster %2u of %2u done: mapped %6u of %6u voxels.\n',...
      c, kmean.C, sum(cl.idx<=c), dim.V);
  elseif bool.chat>=1
    fprintf('done with %6u of %6u voxels.', sum(cl.idx<=c), dim.V);
  end
  
  % store cluster estimate for output
  for l = 1:dim.L
    x{l}(cl.idx==c) = xc{l};
  end
end
t = toc;
end