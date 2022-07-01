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
% Since the camera angle is untrustworthy, the camera angle distribution is 
% set to a uniform distribution. This is done by transforming the normal
% vector to cylindrical coordinates, setting theta randomly and then
% extracting the quaternions again.
%
% Input arguments:
% - quats:      [n x 4] array containing quaternions describing orientations
% - HM:         Cell containing H and M matrix describing particle polytope
% - L:          Particle lengths
%
% Output arguments:
% - newQuats:   [n x 4] array containing new quaternions 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [newQuats] = setCamAngleUniform(quats) 
newQuats = zeros(size(quats));
% Unrotated H-matrix
A = eye(3);
for ii=1:size(quats,1)
    % Rotate with given quaternion
    q = quaternion(quats(ii,:));
    p = quaternion([[0;0;0], A]);
    [~,x,y,z] = parts(q*p*q');
    % Transform to cylindrical coordinates
    [theta, rho]=cart2pol(x, y);
    % Uniformly random theta
    addTheta = rand * 2 * pi;
    theta = theta - min(theta);
    theta = theta + addTheta;
    % Transform modified cylindrical coordinates back to cartesian coordinates
    [x, y] = pol2cart(theta, rho);
    % Transform rotated points into axis/angle representation
    [alpha, w] = getAngleFromPoints(A, [x,y,z]);
    % Quaternion representation
    newQuats(ii,:) = [cos(-alpha/2) sin(-alpha/2)*w];
end

end
