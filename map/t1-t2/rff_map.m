  function z = rff_map(q, rff, w, freq, ph, dim)
%|function z = rff_map(q, rff, w, freq, ph, dim)
%| 
%|  feature mapping via random fourier features
%|
%|  inputs
%|    q         [D+N K]         (lower-dimensional) input [data; known parameters]
%|    rff       [1x1 struct]    random fourier features object
%|     .std     [1]               noise std dev in training data
%|     .len     [D+N]             kernel input length scales
%|     .H       [1]               embedding dimension       
%|     .K       [1]               number of training samples
%|    w         [D]             dataset weights
%|    freq      [H D+N]         random 'frequency' vector
%|    ph        [H]             random phase vector
%|    dim       [1x1 struct]    object containing dimension info
%|
%|  outputs
%|    z         [H K]           (higher-dimensional) features
%|
%|  version control
%|    1.1       2016-09-02      original

% check freq and ph dimensions
if size(freq,1)~=rff.H || length(ph)~=rff.H
  error('Length mismatch: freq and/or ph not of length H!?');
end

% trick: scale lengthscales to re-weight data
%   w>1 -> rff.len lower  -> data emphasized more
%   w<1 -> rff.len higher -> data emphasized less
rff.len(1:dim.D) = div0(rff.len(1:dim.D), max(w,eps));

% random fourier features
tmp = freq * q;
tmp = bsxfun(@plus, tmp, ph);
tmp = cos(2*pi*tmp);
z = tmp * sqrt(div0(2,rff.H));  
end