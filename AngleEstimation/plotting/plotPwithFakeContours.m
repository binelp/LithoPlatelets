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
% Plot a polyhedron with it's contours/projections. 
%
% Input arguments:
% - P:          Polyhedron object
% - colors:     Color
% - alpha:      Transparency of projections
% - dist:       Distance of contours from centerg
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = plotPwithFakeContours(P, colors, alpha, dist)
    xz = convhull(P.V(:,1),P.V(:,3));
    yz = convhull(P.V(:,2),P.V(:,3));
    hold on    ;
    P.plot("color", colors)
    fill3(P.V(xz,1),dist*ones(size(xz,1),1),P.V(xz,3),colors,"facealpha", alpha);
    fill3(dist*ones(size(yz,1),1),P.V(yz,2),P.V(yz,3),colors,"facealpha", alpha);
    axis equal
end