  function rff_map_test(q, rff, w, freq, ph, dim)
%|function rff_map_test(q, rff, w, freq, ph, dim)
%|
%|  test to make sure that rff_map is correct
%|
%|  version control
%|    1.1       2016-09-05      original
%|    1.2       2016-09-16      suppress warning

% rff kernel approximation
z = rff_map(q, rff, w, freq, ph, dim);                                  % [H K]
k.rff = z' * z;                                                         % [K K]

% suppress singularity warning
warning('off', 'MATLAB:nearlySingularMatrix');

% exact gaussian kernel
tmp = rff.len * rff.c;
sqrtSig = spdiag(tmp, 0, size(q,1), size(q,1));
k.ex = NaN(size(k.rff));
for i = 1:size(k.rff,1)
  dq = bsxfun(@minus, q(:,i), q);
  dq = sqrtSig \ dq;
  k.ex(i,:) = exp(-(sum(abs(dq.^2),1))/2);
end
warning('on', 'MATLAB:nearlySingularMatrix');

% compare
str = sprintf('rff.std: %0.6f', rff.std);
figure;
  subplot(1,3,[1,2]); im(cat(3, k.rff, k.ex), [0 1], 'cbar', str); 
  subplot(1,3,3); im(abs(k.rff-k.ex), [0 0.2], 'cbar');
drawnow;  
end