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
% Get the rotation with minimal angle and any axis to get from orientation 1 
% to orientation 2, for a cuboid with given lengths.
% 
% Input arguments:
% - L:          [1 x 3] array containing three particle lengths
% - quats1:     [n x 4] array containing quaternions describing orientation 1
% - quats2:     [n x 4] array containing quaternions describing orientation 2
% - reduceFlag: Wether to calculate reduced quaternion or just return angle
% - considerL:  Wether to make use symmetry caused by L1=L2 or L2=L3 or
%               only consider symmetry applicable to any cuboid
%
% Output arguments:
% - minAlpha:   Minimal angle for rotation
% - minW:       Axis corresponding to minAlpha
% - quatsRed:   Quaternion describing rotation 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [minAlpha, minW, quatsRed] = minAngleBetweenOrientations(L, quats2, quats1, reduceFlag, considerL)
load("CrystalData.mat");
HM = {Crystals(15).H,Crystals(15).M};

% Create cuboid with given lengths
PolyOriginal = Polyhedron(HM{1},HM{2}*(L./2)');
p = quaternion([[0;0;0;0;0;0], PolyOriginal.A]);

for ii=1:length(quats1)
    % Rotate cuboid to obtain orientation 1
    q1 = quaternion(quats1(ii,:));
    [~,x,y,z] = parts(q1*p*q1');
    Arot1 = [x,y,z];

    % Rotate cuboid to obtain orientation 2
    q2 = quaternion(quats2(ii,:));
    [~,x,y,z] = parts(q2*p*q2');
    Arot2 = [x,y,z];
    
    % Vectors describing cuboid orientations 1 and 2
    vecs1or = Arot1([1,3,5],:);
    vecs2or = Arot2([1,3,5],:);
    % Get angle between two closest permutations of +/- vectors
    alpha = [];
    w = zeros([0,3]);
    if ~considerL
        % Maintain chirality of particle coordinate system
        posNegPerms = [1 1 1; unique(perms([-1,-1,1]),"rows")];
        for kk =1:length(posNegPerms)
            vecs2(1,:) = posNegPerms(kk,1).*vecs2or(1,:);
            vecs2(2,:) = posNegPerms(kk,2).*vecs2or(2,:);
            vecs2(3,:) = posNegPerms(kk,3).*vecs2or(3,:);
            [alpha(end+1),  w(end+1,:)] = getAngleFromPoints(vecs1or',vecs2');
        end
    else
        % Use information about equal lengths to detect more equivalent combinations 
        % If needle
        if L(2)==L(3)
            flips = [[0 0]; perms([0,1]); [1 1]];
            for ll=1:4
                vecs2fl = vecs2or;
                vecs1fl = vecs1or;
                if flips(ll,1)
                    vecs2fl(2:3,:) = vecs2or([3,2],:);
                end
                if flips(ll,2)
                    vecs1fl(2:3,:) = vecs1or([3,2],:);
                end
                % Maintain chirality 
                if mod(sum(flips(ll,:)),2)==0
                    posNegPerms = [1 1 1; unique(perms([-1,-1,1]),"rows")];
                else
                    posNegPerms = [unique(perms([-1,1,1]),"rows"); -1 -1 -1];
                end
                for kk =1:length(posNegPerms)
                    vecs2(1,:) = posNegPerms(kk,1)*vecs2fl(1,:);
                    vecs2(2,:) = posNegPerms(kk,2)*vecs2fl(2,:);
                    vecs2(3,:) = posNegPerms(kk,3)*vecs2fl(3,:);
                    [alpha(end+1), w(end+1,:)] = getAngleFromPoints(vecs1fl',vecs2');
                end
            end
        % If equant platelet 
        elseif  L(1)==L(2)
            flips = [[0 0]; perms([0,1]); [1 1]];
            for ll=1:4
                vecs2fl = vecs2or;
                vecs1fl = vecs1or;
                if flips(ll,1)
                    vecs2fl(1:2,:) = vecs2or([2,1],:);
                end
                if flips(ll,2)
                    vecs1fl(1:2,:) = vecs1or([2,1],:);
                end
                % Maintain chirality of particle coordinate system
                if mod(sum(flips(ll,:)),2)==0
                    posNegPerms = [1 1 1; unique(perms([-1,-1,1]),"rows")];
                else
                    posNegPerms = [unique(perms([-1,1,1]),"rows"); -1 -1 -1];
                end
                for kk =1:length(posNegPerms)
                    vecs2(1,:) = posNegPerms(kk,1).*vecs2fl(1,:);
                    vecs2(2,:) = posNegPerms(kk,2).*vecs2fl(2,:);
                    vecs2(3,:) = posNegPerms(kk,3).*vecs2fl(3,:);
                    [alpha(end+1),  w(end+1,:)] = getAngleFromPoints(vecs1fl', vecs2');
                end
            end
        else
            % Maintain chirality of particle coordinate system
            posNegPerms = [1 1 1; unique(perms([-1,-1,1]),"rows")];
            for kk =1:length(posNegPerms)
                vecs2(1,:) = posNegPerms(kk,1).*vecs2or(1,:);
                vecs2(2,:) = posNegPerms(kk,2).*vecs2or(2,:);
                vecs2(3,:) = posNegPerms(kk,3).*vecs2or(3,:);
                [alpha(end+1),  w(end+1,:)] = getAngleFromPoints(vecs1or',vecs2');
            end
        end
    end
    % Find smallest angle
    [minAlpha(ii), ixMin] = min(alpha);
    % If requested, claculate correstonding quaternion, if not, don't to
    % save time
    if reduceFlag
        % Axis corresponding to smallest alpha
       minW(ii,:) = w(ixMin,:);
       % Quaternion representation
       quatsRed(ii,:) = [cos(minAlpha(ii)/2) sin(minAlpha(ii)/2)*minW(ii,:)];
       
    else
       minW(ii,:) = [nan,nan,nan];
       quatsRed(ii,:) = [nan,nan,nan,nan];
    end
end
end