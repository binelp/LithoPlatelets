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
% Reconstructs quaternions (4D) from reduced 3D "quaternions". It is
% based on the fact that for any quaternion q, |q|=1 and q = -q.
%
% Input arguments:
% - quats234:   [N x 3] Reduced 3D "quaternion"
%
% Output arguments:
% - quats:      [N x 4] Reconstructed 4D quaternion
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [quats] = quat3Dto4D(quats234)
    % If slightly negative due to round off error:
    ixOver1 = vecnorm(quats234,2,2) > 1;
    quats234(ixOver1,:) = quats234(ixOver1,:)./vecnorm(quats234(ixOver1,:),2,2);
    
    cubeQuats1 = 1-(quats234(:,1).^2+quats234(:,2).^2+quats234(:,3).^2);
    % If slightly negative due to round off error:
    ixNeg = cubeQuats1 < 0;
    cubeQuats1(ixNeg) = 0;

    quat1 = sqrt(cubeQuats1);
    quats = [quat1 quats234];
end