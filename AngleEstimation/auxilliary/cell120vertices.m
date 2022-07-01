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
% Calculates all the vertices of the 120-cell, which is the 4D equivalent of 
% a dodecahedron: https://en.wikipedia.org/wiki/120-cell, for the purpose
% of systematically sampling quaternions from the surface of the 4D sphere
% as initial points for the orientation-optimization.
%
% Input arguments: None
%
% Output arguments:
% - vert:      [600 x 4] 120-cell vertices
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vert] = cell120vertices() 
phi = GoldenRatio;
% All permutations of:         
allPermsOf = [0, 0, 2, 2;
               1, 1, 1, sqrt(5);
               phi^(-2), phi, phi, phi;
               phi^(-1), phi^(-1), phi^(-1), phi^(2)];

% only even permutations of:  
evenPermsOf = [0, phi^(-2), 1, phi^(2);
                0, phi^(-1), phi, sqrt(5);
                phi^(-1), 1, phi, 2];

[vert] = verticesGeneration4D(allPermsOf, evenPermsOf);
% normalize to radius 1
vert = vert./vecnorm(vert,2,2);
end
