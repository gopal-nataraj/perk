  function [x, t] = plm(y, x0, nu, P, w, mu0, boxcon, reg, stop, dim, bool, mask, disp)
%|function [x, t] = plm(y, x0, nu, P, w, mu0, boxcon, reg, stop, dim, bool, mask, disp)
%|
%|  iterative latent parameter estimation via projected levenberg-marquardt method
%|    uses step-halving line search to ensure monotone descent in cost
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
%|    mu0       [1]             levenberg-marquardt parameter         
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
%|     .mag.*   false|true        using magnitude (spgr,dess) data
%|     .chat    false|true        verbosity                             
%|     .norm    false|true        normalize data for scale-invar reg
%|     .reg     false|true        use regularization                   
%|     .precon  false|true        use preconditioner (hessian approx)   
%|     .disp    false|true        show image updates
%|    mask      [(odims)]       binary object mask
%|    disp      {L cell}        display ranges for image updates             
%|
%|  outputs  
%|    x         {L cell}        latent object parameters
%|    t         [1]             run time
%|  
%|  version control
%|    1.1       2016-06-06      original
%|    1.2       2016-06-13      display range options
%|    1.3       2016-06-27      now handles inveff estimation when ir data passed
%|    1.4       2016-07-02      using new plm_dir(...) for proj lev-marq method
%|    2.1       2017-06-07      project initialization onto box constraint
%|    2.2       2018-01-07      now passes weight vector instead of sparse weight matrix

% improve initialization by projecting onto box constraint
for l=1:dim.L
  x0{l} = max(x0{l}, boxcon{l}(1));
  x0{l} = min(x0{l}, boxcon{l}(2));
end
x = x0;                                                                 % {L} cell w/ [V] blocks

% weighting matrix
W = spdiags(col(repmat(w, [1 dim.V])), 0, dim.D*dim.V, dim.D*dim.V);    % [DV DV sparse]

tic;
for n = 1:stop.iter 
  d = plm_dir(y, x, nu, P, W, dim, reg, mu0, bool);
  psiPrev = Psi(y, x, nu, P, W, dim, reg, bool);
  step = 1;
  if bool.chat && rem(n,1)==0
    fprintf('RLS via proj LM w/ backtrack: cost %0.6e at iter %u of %u:\n',...
      psiPrev, n, stop.iter);
  end
  if isnan(psiPrev)
    error('Cost is NaN!?');
  end
  
  % backtracking line search
  while true
    xv = col(transpose([x{:}]));                                        % [LV]
    xvCand = xv + d * step;                                             % [LV]
    xCand = num2cell(reshape(xvCand, [dim.L dim.V]), 2);                % {L} cell w/ [1 V] blocks 
    
    % project onto box contraints (omit m0)
    for l = 1:dim.L
      if l==1
        xCand{l} = col(xCand{l});                                       % m0 now w/ [V] block
      else
        xvCand(l:dim.L:end) = max(xvCand(l:dim.L:end), boxcon{l}(1));
        xvCand(l:dim.L:end) = min(xvCand(l:dim.L:end), boxcon{l}(2));
        xCand{l} = max(col(xCand{l}), boxcon{l}(1));                    % {L} cell w/ [V} blocks
        xCand{l} = min(col(xCand{l}), boxcon{l}(2));                    % {L} cell w/ [V] blocks
      end
    end
    psiCand = Psi(y, xCand, nu, P, W, dim, reg, bool);
    
    % if cost descends, ready to proceed to next iteration
    if psiCand < psiPrev
      break;
      
    % else halve the step size and try to reduce cost again
    % if step gets too small, return current x
    else
      if bool.chat && rem(n,1)==0
        fprintf('\tCost %0.6e; step halved from %0.2e to %0.2e.\n', psiCand, step, step/2);
      end
      if isnan(psiCand)
        error('Cost is NaN!?');
      end
      step = step/2;
      
      % if step gets too small, return current x
      num = wnorm2(xvCand, stop.wghtx, dim);
      den = wnorm2(d, stop.wghtx, dim);
      stepTol = stop.tolx * div0(num, den);
      if step < stepTol
        t = toc;
        if bool.chat
          if rem(n,1)~=0
            fprintf('RLS via proj LM w/ backtrack: cost %0.6e at iter %u of %u:\n',...
              psiPrev, n, stop.iter);
          end
          fprintf('\tReturn early: step < stepTol = %0.6e.\n',  stepTol);
        end
        return;
      end
    end
  end

  % update x
  x = xCand;
  
  % return early if minimal change in solution
  num = wnorm2(xvCand-xv, stop.wghtx, dim);
  den = wnorm2(xv, stop.wghtx, dim);
  xvChange = div0(num,den);
  if xvChange < stop.tolx
    t = toc;
    if bool.chat
      fprintf('Return early: norm(xnew-xold)/norm(xnew) < xTol = %0.3e.\n', stop.tolx);
    end
    return;
  end
  
  % display images
  if bool.disp && rem(n,1)==0
    % if disp{1} (m0) unset, fix once
    clf;
    if isempty(disp{1})
      disp{1} = [0 2*median(abs(x{1}))];
    end
    figure(1); 
    im(embed(abs(x{1}), mask), disp{1}, 'cbar'); 
    drawnow;
    for l=2:dim.L
      figure(l); 
      im(embed(x{l}, mask), disp{l}, 'cbar'); 
      drawnow;
    end
  end
end
t = toc;
end