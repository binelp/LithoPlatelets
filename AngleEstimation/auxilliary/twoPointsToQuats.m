%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ETH Zurich, Switzerland
% Separation Processes Laboratory
%
% Project:  Lithoplatelets
% Year:     2021
% MATLAB:   R2019b, Windows 64bit
% Authors:  Anna Jaeggi (AJ)
%
% Purpose:
% Quaternion describing rotation from point P1 to P2 (both P1 and P2 should be 3D unit vectors)
%
% Input arguments:
% - P1, P2:     Two 3D points
%
% Output arguments:
% - quats:   [n x 4] array containing new quaternions 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [quats] = twoPointsToQuats(P1,P2)
axRot = cross(P1,P2,2);
axRot = axRot./norm(axRot);
alpha = acos(dot(P1,P2,2));
quats = [cos(alpha), sin(alpha).*axRot];
quats = quats./vecnorm(quats,2,2);
end