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
% Calculates all the vertices of the 600-cell (which is the 4D equivalent of 
% an icosahedron: https://en.wikipedia.org/wiki/600-cell) for the purpose
% of systematically sampling quaternions from the surface of the 4D sphere
% as initial points for the orientation-optimization. The original 600-cell
% consist of tetrahedrons which can be further subdivided into smaller
% tetrahedrons to generate more points (similarly to triangle-subdivision of 
% an icosahedron for a geodesic grid covering a 3D sphere surface)
%
% Input arguments: 
% - subdivision:    How often are the tetrahedron cells subdivided into
%                   smaller tetrahedrons to generate more vertices
%
% Output arguments:
% - vert:      [n x 4] 600-cell vertices 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vert] = cell600vertices(subdivisions) 

phi = GoldenRatio;
% All permutations of:         
allPermsOf = [0,0,0,1;
               0.5,0.5,0.5,0.5];

% only even permutations of:  
evenPermsOf = [phi/2, 0.5, (phi^-1)/2, 0];

[vert] = verticesGeneration4D(allPermsOf, evenPermsOf);
if subdivisions
    % subdivide
    for subdivIx=1:subdivisions
        % find tetrahedrons composing 600-cell through distance

        % initialize array of indexes marking tetrahedrons
        tetIxs = zeros([0,4]);
        tetEdgeLen = (1/phi)/subdivIx;
        % loop through all points
        for ix1=1:length(vert)
            dists = vecnorm(vert(ix1,:)-vert,2,2);
            otherIx = setdiff(1:size(vert,1), ix1);
            dists = dists(otherIx);
            minIx1 = find(abs(dists - tetEdgeLen) < 10e-5);
            minIx1 = otherIx(minIx1);

            % ii = 1;
            for ii = 1:length(minIx1)
                ix2 = minIx1(ii);
                otherIx2 = setdiff(minIx1, ix2);
                dists2 = vecnorm(vert(ix2,:)-vert(otherIx2,:),2,2);
                minIx2 = find(abs(dists2 - tetEdgeLen) < 10e-5);
                minIx2 = otherIx2(minIx2);

                for jj = 1:  length(minIx2)
                    ix3 = minIx2(jj);
                    otherIx3 = setdiff(minIx2, ix3);
                    dists3 = vecnorm(vert(ix3,:)-vert(otherIx3,:),2,2);
                    minIx3 = find(abs(dists3 - tetEdgeLen) < 10e-5);
                    minIx3 = otherIx3(minIx3);

                    % kk = 2;
                    for kk = 1:length(minIx3)
                        ix4 = minIx3(kk);
                        % save
                        tetIxs(end+1,:) = [ix1,ix2,ix3,ix4];

                        tet = vert(tetIxs(ii,:),:);
                        [isTet] = checkIfTet(tet);
                        if ~isTet
                            break
                        end

                    end
                end
            end
        end
        % remove duplicates (order doesnt matter)
        tetIxs = sort(tetIxs,2);
        tetIxs = unique(tetIxs,"rows");
        % original 600-cell should consist of 600 tetrahedrons
        if subdivIx ==1 && length(tetIxs)~= 600
            break
        end
        % calculate new vertices of sub-tetrahedra for each tetrahedron
        allNewVerts = zeros([0,4]);
        for ii=1:length(tetIxs)
             tet = vert(tetIxs(ii,:),:);
            [newVerts] = subdivideTet(tet);
            allNewVerts(end+1:end+6,:) = newVerts;
        end
        % remove duplicates (order matters)
        allNewVerts = unique(allNewVerts,"rows");

        vert = cat(1,vert,allNewVerts);
    end
end
end
function [isTet] = checkIfTet(tet)
edgeLen = norm(tet(1,:) - tet(2,:));
isTet = true;
for a = 1:4
    for b = 1:4
        if a ~=b
            if abs(norm(tet(a,:) - tet(b,:)) - edgeLen)>10e-5
                isTet = false;
                break
            end   
        end
    end
end
    
end
function [newVerts] = subdivideTet(tet)
edgeLen = norm(tet(1,:) - tet(2,:));
% subdivide (get all half points of tetrahedron edges)
newVerts = zeros([0,4]);
for ii=1:4
    otherIx = setdiff(1:4,ii);
    for jj = 1:3
        newVerts(end+1,:) = tet(ii,:) + (tet(otherIx(jj),:)-tet(ii,:))*1/2;
    end
end
newVerts = unique(newVerts,"rows");
end