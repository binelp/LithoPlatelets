%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ETH Zurich, Switzerland
% Separation Processes Laboratory
%
% Project:  Lithoplatelets
% Year:     2021
% MATLAB:   R2019b, Windows 64bit
% Authors:  Anna Jaeggi (AJ)
%
% Purpose:
% Convert a list of angles or reduced quaternions to a distribution of
% angles/quaternions, or a probability density.
% 
% Input arguments: 
% - angles:         [n x 3] Array containing either euler angles or reduced
%                   quaternions.
% - lims:           [3 x 2] Array containing lower and upper bounds of
%                   distribution.
% - bins:           [1 x 3] Array containing number of bins along each
%                   dimension.
% - normalizeFlag:  Bool indicating wether to normalize at all, to the
%                   maximum value, or to the integral.
%  
% Output arguments:
% - AD:             Struct containing the density (F) and the bin
%                   boundaries and midpoints (grid).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [AD] = angleListToDist(angles, lims, nbins, normalizeFlag)
F = zeros(nbins);
% Create regular grid according lims and bins
for ii=1:3
    setGrid(ii).boundaries = linspace(lims(ii,1),lims(ii,2),nbins(ii)+1);
    setGrid(ii).y = (setGrid(ii).boundaries(1:end-1)+setGrid(ii).boundaries(2:end))/2;
end
AD.grid = setGrid;
% Count all angles in list in each of the bins
for ii = 1:nbins(1)
    for jj = 1:nbins(2)
        for kk = 1:nbins(3)
            
            ix1 = angles(:,1) < setGrid(1).boundaries(ii+1) & angles(:,1) > setGrid(1).boundaries(ii);
            ix2 = angles(:,2) < setGrid(2).boundaries(jj+1) & angles(:,2) > setGrid(2).boundaries(jj);
            ix3 = angles(:,3) < setGrid(3).boundaries(kk+1) & angles(:,3) > setGrid(3).boundaries(kk);
            ix = ix1 & ix2 & ix3;
            
            F(ii,jj,kk) = sum(ix);
        end
    end
end
% Normalization
if normalizeFlag == 1
    % Normalize to have integral 1 (probability)
    AD.F = F/sum(F,"all");
elseif normalizeFlag == 2
    % Normalize to have maximum value of 1
    AD.F = F/max(F,[],"all");
else
    % Don't normalize
    AD.F = F;
end
end