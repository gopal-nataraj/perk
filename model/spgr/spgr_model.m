function [Smodel, Sdata, adata] = spgr_model(M, E1, y, cent, rad, amodel, flip)
% function [Smodel, Sdata] = signal_model(M, E1, y, cent, rad)
% Inputs:
%   M:      [nx ny]         Spin density
%   E1:     [nx ny]         exp(-TR/T1)
%   y:      [nx ny na]      Raw signal
%   cent:   [1 2]           Coordinates of interest
%   rad:    [1]             Radius of ROI
%   amodel  [nmod 1]        Desired fitting to model
%   flip    [nx ny na]      Flip angle map
% Outputs: 
%   Smodel: [nmod 1]        Model fit to parameters
%   Sdata:  [(2*rad+1)^2 na]Raw data
%   adata:  [na 1]          Mean flip angle over ROI

% Gather ROI data
M_roi = col(M(cent(1)-rad:cent(1)+rad, cent(2)-rad:cent(2)+rad));
E1_roi = col(E1(cent(1)-rad:cent(1)+rad, cent(2)-rad:cent(2)+rad));
y_roi = col(y(cent(1)-rad:cent(1)+rad, cent(2)-rad:cent(2)+rad, :));

% Take means of ROI
Mbar = mean(M_roi);
Ebar = mean(E1_roi);

% Get expected signal
Smodel = Mbar * (1-Ebar) * sin(amodel) ./ (1-Ebar*cos(amodel));

% Extract test data
adata = squeeze(mean(mean(flip(cent(1)-rad:cent(1)+rad, ...
    cent(2)-rad:cent(2)+rad, :), 1), 2));
    
Sdata = reshape(y_roi, [(2*rad+1)^2 size(y,3)]);