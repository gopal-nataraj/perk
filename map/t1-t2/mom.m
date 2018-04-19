  function [x, t] = mom(y, x0, nu, P, w, dim, bool, mask)
%|function [x, t] = mom(y, x0, nu, P, w, dim, bool, mask)
%|
%|  method-of-moments t1/t2 estimation from spgr/dess scans
%|    unless flip angles are uniformly pi/2, t2 estimates will be biased for finite t1
%|    requires fixed tr,te for all spgr/dess scans
%|    neglects b0, r2p effects
%|  
%|  inputs
%|    y         [D V]           coil-combined image data
%|    x0        {1-L cell}      latent parameter initial estimates
%|                                any subset of {m0,t1,t2}  
%|    nu        {N cell}        known parameters
%|    P         [1x1 struct]    scan parameters
%|                                1st field: data type (sp,de)
%|                                2nd field: scan parameter, as appropriate
%|     .*.tr    [S.*]               repetition times                                          ms
%|     .*.te    [S.* E.*]           echo times                                                ms
%|     .*.aex   [S.*]               nominal flip angle of excitation                          rad
%|    w         [D]             dataset weights
%|    dim       [1x1 struct]    object containing dimension info
%|    bool      [1x1 struct]    boolean variables          
%|     .mag.*   false|true        using magnitude (spgr,dess) data
%|     .chat    false|true        verbosity  
%|     .m0dess  false|true        use dess data for m0 estimation (will cause m0 bias)
%|  
%|  outputs
%|    x         {L cell}        latent object parameter estimates
%|    t         [1]             run time
%|
%|  copyright 2017, gopal nataraj, university of michigan
%|
%|  version control
%|    1.1       2017-12-19      original
%|    1.2       2017-12-21      now ensures t1,t2 are pure-real

% check image data is pure real
if ~isreal(y)
  warn('Detected complex images; for method-of-moments estimation; using magnitude...');
  y = abs(y);
end

% start timer
tic;
x = cell(dim.L,1);

% t1 estimation
if ~isempty(x0{2})
  x{2} = masker(x0{2}, mask);
else
  % check spgr repetition times all equal
  tmp = bsxfun(@minus, P.sp.tr, P.sp.tr(1));
  if sum(tmp)>0
    error('All spgr repetition times must be equal!');
  end
  
  % check spgr echo times all equal
  tmp = bsxfun(@minus, P.sp.te, P.sp.te(1));
  if sum(tmp)>0
    error('All spgr echo times must be equal!');
  end
  
  % preprocess spgr data
  tmp = P.sp.aex * transpose(nu{1});                                    % [Dsp V]
  tmp2 = col(y(1:dim.Dsp,:) ./ sin(tmp));                               % [Dsp*V]
  tmp3 = reshape(y(1:dim.Dsp,:) ./ tan(tmp), [dim.Dsp 1 dim.V]);                       
  tmp3 = cat(2, tmp3, ones(size(tmp3)));                                % [Dsp 2 V]
  tmp3 = spmat(tmp3);                                                   % [Dsp*V 2V sparse]
  
  % spgr scan weights
  tmp4 = reshape(w(1:dim.Dsp), [dim.S.sp dim.E.sp]);
  tmp4 = mean(tmp4,2);
  tmp4 = spdiags(col(repmat(tmp4, [1 dim.V])),...
    0, dim.S.sp*dim.V, dim.S.sp*dim.V);
  
  % affine regression
  tmp5 = tmp3'*(tmp4*tmp2);                                             % [2V]
  tmp6 = tmp3'*(tmp4*tmp3);                                             % [2V 2V sparse]
  tmp7 = tmp6 \ tmp5;                                                   % [2V]
  tmp7 = reshape(tmp7, [2 dim.V]);                                      % [2 V]
  tmp8 = bsxfun(@max, col(tmp7(1,:)), eps);                             % [V]
  x{2} = div0(-ones(dim.V,1)*P.sp.tr(1),log(tmp8));                     % [V]
end

% t2 estimation
if ~isempty(x0{3})
  x{3} = masker(x0{3}, mask);
else
  % check dess repetition times all equal
  tmp = bsxfun(@minus, P.de.tr, P.de.tr(1));
  if sum(tmp)>0
    error('All dess repetition times must be equal!');
  end
  
  % check dess echo times all equal
  tmp = bsxfun(@minus, P.de.te, P.de.te(1,1));
  if sum(col(tmp))>0
    error('All dess echo times must be symmetric and equal!');
  end
  
  % preprocess dess data
  tmp = reshape(y(dim.Dsp+1:end,:), [dim.S.de dim.E.de dim.V]);
  tmp2 = col(tmp(:,2,:));                                               % [S.de*V]
  tmp3 = reshape(tmp(:,1,:), [dim.S.de 1 dim.V]);
  tmp3 = spmat(tmp3);                                                   % [S.de*V V sparse]
  
  % dess scan pair weights set as average
  tmp4 = reshape(w(dim.Dsp+1:end), [dim.S.de dim.E.de]);
  tmp4 = mean(tmp4,2);
  tmp4 = spdiags(col(repmat(tmp4, [1 dim.V])),...
    0, dim.S.de*dim.V, dim.S.de*dim.V);                                 % [S.de*V S.de*V sparse]
  
  % linear regression
  tmp5 = tmp3'*(tmp4*tmp2);                                             % [V]
  tmp6 = tmp3'*(tmp4*(tmp3*ones(dim.V,1)));                             % [V]
  tmp7 = div_safe(tmp5, tmp6);                                          % [V]
  tmp7 = bsxfun(@max, tmp7, eps);                                       % [V]
  
  % t2 estimation
  tmp8 = -2*(P.de.tr(1)-P.de.te(1,1));
  x{3} = div0(ones(dim.V,1)*tmp8,log(tmp7));                            % [V]
end

% m0 estimation via least-squares
if ~isempty(x0{1})
  x{1} = masker(x0{1}, mask);
else
  % omitting dess data will cause less bias but higher variance
  if ~bool.m0dess
    w(dim.Dsp+1:end) = 0;
  end
  
  % weight matrix
  W = spdiags(col(repmat(w, [1 dim.V])),...
    0, dim.D*dim.V, dim.D*dim.V);                                       % [DV DV sparse]
  
  % get system matrix only for spgr/dess scans
  dim.S.ir = 0;
  dim.S.se = 0;
  dim.D = dim.S.sp*dim.E.sp + dim.S.de*dim.E.de;
  x{1} = nan(dim.V,1);
  tmp = rgradx_f(x, nu, P, dim, bool);                                 
  Am = tmp(:,1:dim.L:end);
  
  % linear regression
  num = Am'*W*col(y);                                                   % [V]
  den = (Am'*W*Am)*ones(dim.V,1);                                       % [V]
  x{1} = div_safe(num, den);                                            % [V]
end

% record timer
t = toc;
  end