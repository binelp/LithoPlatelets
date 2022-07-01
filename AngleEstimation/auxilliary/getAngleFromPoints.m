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
% Get the axis/ angle representation or a quaternion from two sets of
% points P1, P2, P3 and Q1, Q2, Q3, where Qi are the rotated versions of
% Pi, meaning Pi*R = Qi. Helpful links:
% https://math.stackexchange.com/questions/222113/given-3-points-of-a-rigid-body-in-space-how-do-i-find-the-corresponding-orienta
% https://math.stackexchange.com/questions/893984/conversion-of-rotation-matrix-to-quaternion
% https://en.wikipedia.org/wiki/Axis%E2%80%93angle_representation
% 
%
% Input arguments:
% - P:      [3 x 3] array P1,P2,P3 in the ROWS
% - Q:      [3 x 3] array Q1,Q2,Q3 in the ROWS
%                                  P/Q = [x1, y1, z1;
%                                         x2, y2, z2;
%                                         z3, y3, z3]
% Output arguments:
% - alpha:  double rotation angle 
% - w:      [1 x 3] array containing axis of rotation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [alpha, w] = getAngleFromPoints(P,Q)
% Rotation matrix (since P*R = Q)
R = Q*inv(P); 
% Rotation angle 
alpha = acos((trace(R)-1)/2);
% Axis of rotation
% To avoid numerical issues, split cases
nearPi = abs(alpha - pi)<10e-6;
nearZero = abs(alpha - 0)<10e-6;
if ~nearPi && ~nearZero
    w = 1/(2*sin(alpha))*[R(3,2)-R(2,3), R(1,3)-R(3,1), R(2,1)-R(1,2)];
else
    B = 1/2* (R + eye(3));
    wSign = sign(1/(2*sin(alpha))*[R(3,2)-R(2,3), R(1,3)-R(3,1), R(2,1)-R(1,2)]);
    w = [sqrt(B(1,1)) sqrt(B(2,2)) sqrt(B(3,3))];
    w = wSign.* w;
end
end