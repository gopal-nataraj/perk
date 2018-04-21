  function [Xkk_inv] = diag_pinv(X, k, max_cond_num)
% function [Xkk_inv] = diag_pinv(X, k, max_cond_num)
% diag_pinv.m   Efficiently computes diagonal elements of an inverse
%
% Inputs:
%   X           [n n]       Symmetric matrix
%   k           [1]         Diagonal element of inverse desired
%   max_val     [1]         Threshold at which pseudoinverse sets to zero
%
% Outputs:
%   Xkk_inv     [1]         pinv(X)_{k,k}
%
% Written by Gopal Nataraj; Copyright 2014

% Extract sub-block, X_k
X_k = [X(1:k-1,1:k-1) X(1:k-1,k+1:end); X(k+1:end,1:k-1) X(k+1:end,k+1:end)];
a_k = [X(1:k-1,k); X(k+1:end,k)];

% Submatrix inversion
Xkk_inv = 1 / (X(k,k) - a_k' * pinv(X_k) * a_k);
if abs(Xkk_inv) > max_cond_num
    Xkk_inv = Inf;          % True inverse
    % Xkk_inv = 0;            % Pseudoinverse 
end