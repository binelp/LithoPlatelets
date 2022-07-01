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
% Calculates vertices of 4D polytopes according to given permutations (can be found on wikipedia)
%
% Input arguments: 
% - allPermsOf:     [n x 4] Array containing points of which to take all 
%                   permutations
% - evenPermsOf:    [n x 4] Array containing points of which to take only even
%                   permutations
%
% Output arguments:
% - vert:      [n x 4] polytope vertices 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vert] = verticesGeneration4D(allPermsOf, evenPermsOf)
% All positive/ negative conbinations of 4 elements
posNegPerms = [[1 1 1 1];
                unique(perms([1,1,1,-1]),"rows");
                unique(perms([1,1,-1,-1]),"rows");
                unique(perms([1,-1,-1,-1]),"rows");
                [-1,-1,-1,-1]];   
            
% All permutations:      
vert1=zeros([0,4]);        
for ii=1:size(allPermsOf,1)
    % Generate all possible combinations of +/- vec
    posNegCombs = unique(allPermsOf(ii,:).* posNegPerms,"rows");
    for jj=1:size(posNegCombs,1)
        % Generate all possible permutations of one vec
        vert1 = cat(1, vert1, unique(perms(posNegCombs(jj,:)),"rows"));
    end
end
% Only even permutations:  
allPerms = perms([1 2 3 4]);
for ii=1:length(allPerms)
    I4 = eye(4);
    a = allPerms(ii,:);
    sign = det(I4(:,a));
    even(ii) = sign > 0;
end
evenPerms = allPerms(even,:);

vert2=zeros([0,4]);
for ii=1:size(evenPermsOf,1)
    % Generate all possible combinations of +/- vec
    posNegCombs = unique(evenPermsOf(ii,:).* posNegPerms,"rows");
    for jj=1:size(posNegCombs,1)
        % Generate all even permutations of one vec
        bla2 = posNegCombs(jj,:);
        vert2 = cat(1,vert2,bla2(evenPerms));
    end
end
% Combine all points
vert = cat(1, vert1, vert2);
% Remove duplicates
vert = unique(vert,"rows");
end