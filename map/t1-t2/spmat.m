  function [spX] = spmat(fullX)
%|function [spX] = spmat(fullX)
%|  
%|  efficiently construct sparse block-diagonal matrix
%|    useful when C >> A*B
%|    
%|  input
%|    fullX     [A B C]         full version
%|
%|  output
%|    spX       [AC BC sparse]  sparse block-diagonal version
%|
%|  version control
%|    1.1       2016-10-27      original

% array size
[A, B, C] = size(fullX);

% index/val instantiation
T = numel(fullX);
i = NaN(T,1);
j = NaN(T,1);
s = NaN(T,1);

% populate 
row = 1;
for a = 1:A
  for b = 1:B
    i(row : row+C-1) = (a : A : A*C)';                                  % [C]
    j(row : row+C-1) = (b : B : B*C)';                                  % [C]
    s(row : row+C-1) = fullX(a,b,:);                                    % [C]
    row = row + C;
  end
end

% sparse matrix
spX = sparse(i, j, s, A*C, B*C, T);                                     % [AC BC sparse]
end