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
% Calculate cross-moments of 3D angle/ quaternion distribution (Adapted from 
% 3D moment calculation KIDDO or Shape)
%
% Input arguments:
% - AD:         Angle distribution struct containing density (F) and grid
%               boundaries and mid-points (grid. boundaries and grid.y)
% - maxMom:     Highest order moment to be calculated
%
% Output arguments:
% - muMat:      Matrix containing cross moments (i.e. mu_000 = muMat(1,1,1))
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function muMat = calcMomentsAD(AD, maxMom)
f = permute(AD.F,[3,2,1]);
% Moments
iMax = maxMom(1);
jMax = maxMom(2);
kMax = maxMom(3);
muMat = zeros(iMax+1,jMax+1,kMax+1);
delL1 = AD.grid(1).boundaries(2)-AD.grid(1).boundaries(1);
delL2 = AD.grid(2).boundaries(2)-AD.grid(2).boundaries(1);
delL3 = AD.grid(3).boundaries(2)-AD.grid(3).boundaries(1);

L1 = AD.grid(1).y;
L1=reshape(L1,1,1,length(L1)); % L1 needs to be a vector pointing to the 3rd dimension
L2 = AD.grid(2).y;% L2 needs to be a ROW vector!
L3 = AD.grid(3).y;
L3 = L3.'; % L3 needs to be a COLUMN vector!

for i=0:iMax
    for j=0:jMax
        for k=0:kMax
            % The line below requires >R2016b (for element-wise multiplication of
            % matrices and vectors, which corresponds to column- or row-wise multiplication).
            % For older releases, have a look at the function bsxfun(.), which does the same.
            muMat(i+1,j+1,k+1) = sum(sum(sum((((f.*(L3.^k)).*(L2.^j)).*L1.^i).*delL1,3).*delL2,2).*delL3);
        end
    end
end
end