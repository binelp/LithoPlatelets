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
% Maps quaternions (4D) to unique 3D number for easier visualization. It is
% based on the fact that for any quaternion q, |q|=1 and q = -q.
%
% Input arguments:
% - quat:       [N x 4] Original 4D quaternion
%
% Output arguments:
% - quatRed:    [N x 3] Reduced 3D "quaternion"
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [quatRed] = quat4DTo3D(quat)
    quat = quat./vecnorm(quat,2,2);
    % Index of dimension to be fixed positive (removes equivalent quats and reduces dims)
    posIx = 1;
    otherIx = setdiff([1:4],posIx);
    
    % Reduce dimensions
    ix = quat(:,posIx)<0;
    quat(ix,:) = -quat(ix,:);
    quatRed = quat(:,otherIx);
end