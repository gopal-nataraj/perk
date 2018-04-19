  function [rgrad] = rgradx_f(x, nu, P, dim, bool)
%|function [rgrad] = rgradx_f(x, nu, P, dim, bool)
%|
%|  signal model row gradient evaluation w.r.t. latent object parameters x
%|  
%|  inputs
%|    x         {L cell}        latent object parameters with cells size [V]
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
%|     .mag.*   false|true        using magnitude (ir,se,sp,de) data     
%|
%|  outputs
%|    rgrad     [DV LV sparse]  row gradient evaluation
%|
%|  version control
%|    1.1       2016-06-02      original
%|    1.2       2016-06-06      added magnitude signal options
%|    1.3       2016-06-08      compatibility with spgr/dess
%|    1.4       2016-06-27      new se-ir signal gradient
%|    1.5       2016-06-28      nu now expected as cell (previously struct)

% initialize with tensor format for convenient indexing
dim.Vloc = numel(x{1});
rgrad = zeros(dim.D, dim.Vloc, dim.L);                                  % [D Vloc L]
row = 1;

% se-ir function calls
for s = 1:dim.S.ir
  rgrad(row,:,:) = IR_gradx_v2(...
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
% note that if dim.L>3, derivative w.r.t. inveff left as zero
for s = 1:dim.S.se
  rgrad(row,:,1:3) = SE_gradx(...
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
% note that if dim.L>3, derivative w.r.t. inveff left as zero
for s = 1:dim.S.sp
  rgrad(row,:,1:3) = spgr_gradx(...
    x{1:3},...
    P.sp.aex(s), P.sp.tr(s), P.sp.te(s,1),...
    'kap', nu{1},...
    'dw', nu{2},...
    'R2p', nu{3},...
    'mag', bool.mag.sp);
  row = row+1;
end

% dess function calls
% note that if dim.L>3, derivative w.r.t. inveff left as zero
for s = 1:dim.S.de
  [rgrad(row,:,1:3), rgrad(row+dim.S.de,:,1:3)] = dess_gradx(...
    x{1:3},...
    P.de.aex(s), P.de.tr(s), P.de.te(s,1), P.de.te(s,2),...
    'kap', nu{1},...
    'dw', nu{2},...
    'R2p', nu{3},...
    'mag', bool.mag.de);
  row = row+1;
end

% construct sparse block-diagonal gradient matrix
rgrad = permute(rgrad, [1 3 2]);                                      % [D L Vloc]
rgrad = spmat(rgrad);                                                 % [DVloc LVloc sparse]
end