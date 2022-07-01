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
% Samples uniformly distributed quaternions from the space of all possible
% rotations by using rejection method sampling of the surface of a 4D sphere.
%
% Input arguments:
% - N:      Desired number of samples
%
% Output arguments:
% - pts:    [N x 4] Array containing sampled quaternions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pts] = uniformSampledQuats(N)
% Factor by which to oversample
fx = 1/N*100 + 1;
f = fx*2*2*2*2/(4/3*pi);
% Random points in 4D cube between -1,1
pts = (rand([ceil(N*f),4])-0.5)*2;
% Remove any outside 4D sphere (to avoid corner over density)
ix = (vecnorm(pts,2,2)<=1);
pts = pts(ix,:);
% Project all points inside 4D sphere onto 4D sphere surface
pts = pts./vecnorm(pts,2,2);
% If too many, delete a few
if length(pts) > N
    pts = pts(1:N,:);
end
end