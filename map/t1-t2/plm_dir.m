  function [d] = plm_dir(y, x, nu, P, W, dim, reg, mu0, bool)
%|function [d] = plm_dir(y, x, nu, P, W, dim, reg, mu0, bool)
%|
%|  computes a descent direction for projected levenberg-marquardt
%|    should only be used with a line search
%|
%|  version control
%|    1.1       2016-06-06      original
%|    1.2       2016-06-28      corrected magnitude check 
%|    1.3       2016-07-02      added lev-marq mu param
%|    1.4       2016-09-17      added gradient check

% gradient
cgradx = Psi_cgradx(y, x, nu, P, W, dim, reg, bool);                    % [LV]
tmp = isinf(cgradx);
if sum(tmp)>0
  error('gradient has %u infinite entries!?', sum(tmp));
end

if bool.precon
  H = Psi_hessx(y, x, nu, P, W, dim, reg, bool);                        % [LV LV sparse]
  mu = mu0 * norm(cgradx,2);
  H = H + mu*speye(dim.L*dim.V);
  
  % if data/m0 all pure-real, can use backslash directly
  tmp = ...
    (dim.S.ir>0 && ~bool.mag.ir) ||...
    (dim.S.se>0 && ~bool.mag.se) ||...
    (dim.S.sp>0 && ~bool.mag.sp) ||...
    (dim.S.de>0 && ~bool.mag.de);
  if ~tmp
    d = -(H \ cgradx);                                                  % [LV]
  
  % else, descent direction is complex only for m0
  % real and imag components computed separately
  else
    % analogous calculation for real component
    dir.r = -(H \ real(cgradx));                                        % [LV]
    
    % extract only m0 portion of H and cgradx for imag component
    Hm = H(1:dim.L:end, 1:dim.L:end);                                   % [V V sparse]
    dm = -(Hm \ imag(cgradx(1:dim.L:end)));                             % [V]
    dir.i = col(transpose([dm zeros(dim.V, dim.L-1)]));                 % [LV]

    % complex output
    d = complex(dir.r, dir.i);                                          % [LV]
  end
else
  d = -cgradx;                                                          % [LV]
end
end