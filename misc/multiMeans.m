function [means, SDs] = multiMeans(obj, centers, rad)
%MULTIMEANS.M Computes mean and SD over several ROIs over an obj
% Inputs:
%   obj         [nx ny]     object of interest
% 	centers     [M 2]       list of center coordinates, where M is #ROIs
%   rad         [1]         determines width of ROIs: (2*rad+1) x (2*rad+1)
% Output:
%   means       [M]         means over M ROIs
%   SDs         [M]         standard deviations over M ROIs
% 
% Written by: Gopal Nataraj
% Copyright 2013

% Variables
diam = 2*rad+1;
M = size(centers, 1);
ROI = zeros(diam, diam, M);

% Extract the regions of interest
for ii = 1:size(centers, 1)
    ROI(:,:,ii) = obj(centers(ii,1)-rad:centers(ii,1)+rad, ...
        centers(ii,2)-rad:centers(ii,2)+rad);
end

% Compute means and SDs
means = mean(reshape(ROI, [diam^2 M]), 1);
SDs = std(reshape(ROI, [diam^2 M]), 0, 1);
%means = squeeze(mean(mean(ROI, 1), 2));