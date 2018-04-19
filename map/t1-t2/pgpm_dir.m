  function [d, D] = pgpm_dir(y, x, nu, P, W, dim, reg, lam, bool, varargin)
%|function [d, D] = pgpm_dir(y, x, nu, P, W, dim, reg, lam, bool, varargin)
%|
%|  computes a descent direction for gradient projection method
%|    includes optional diagonal preconditioner
%|  
%|  version control
%|    1.1       2016-06-06      original
%|    1.2       2016-06-28      corrected magnitude check 
%|    2.1       2017-06-06      better preconditioner from mri_multicomp_map(...)

% default values
arg.D = [];
arg = vararg_pair(arg, varargin);

% gradient
cgradx = Psi_cgradx(y, x, nu, P, W, dim, reg, bool);                    % [LV]
tmp = isinf(cgradx);
if sum(tmp)>0
  error('gradient has %u infinite entries!?', sum(tmp));
end

if bool.precon
  % if D not preset, set to diagonal of hessian
  if isempty(arg.D)
    H = Psi_hessx(y, x, nu, P, W, dim, reg, bool);                      % [LV LV sparse]
    tmp = full(diag(H));                                                % [LV]
    tmp = bsxfun(@plus, tmp, lam);                                   
    tmp = div0(1,tmp);
    arg.D = spdiags(tmp, 0, dim.L*dim.V, dim.L*dim.V);                  % [LV LV sparse]
  end
  
  % if data/m0 all pure-real, can use backslash directly
  tmp = ...
    (dim.S.ir>0 && ~bool.mag.ir) ||...
    (dim.S.se>0 && ~bool.mag.se) ||...
    (dim.S.sp>0 && ~bool.mag.sp) ||...
    (dim.S.de>0 && ~bool.mag.de);
  if ~tmp
    d = -(arg.D * cgradx);                                              % [LV]
    
  % else, descent direction is complex only for m0
  % real and imag components computed separately
  else
    dir.r = -(arg.D * real(cgradx));                                    % [LV]
    
    % extract only m0 portion of D and cgradx for imag component
    Dm = arg.D(1:dim.L:end, 1:dim.L:end);                               % [V V sparse]
    dm = -(Dm * imag(cgradx(1:dim.L:end)));                             % [V]
    dir.i = col(transpose([dm zeros(dim.V, dim.L-1)]));                 % [LV]
    
    % complex output
    d = complex(dir.r, dir.i);                                          % [LV]
  end
else
  d = -cgradx;
end

% set D, if needed
if nargout>1
  D = arg.D;
end
end