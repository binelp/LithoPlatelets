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
% Calculates some angles with a physical meaning based on the
% quaternion representation of particle orientation in order to analyze
% more abstract quaternion distributions.
%
% Input arguments: 
% - quats:          [n x 4] Array containing quaternions
% - L:              [1 x 3] Particle lengths
% - savePath:       Path where results should be saved (can be false/zero)
%
% Output arguments:
% - physicalAngles: Struct containing physicalAngles
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [physicalAngles] = getPhysicalAngles(quats, L, savePath)
% Number of orientations 
N = size(quats,1);
% Load cuboidal morphology
load("CrystalData.mat");
HM = {Crystals(15).H,Crystals(15).M};
% Only plot particles if low number of orientations
if N<4
    plotFlag = true;
else
    plotFlag = false;
end
% Generate unrotated polytope
PolyOriginal = Polyhedron(HM{1},HM{2}*(L./2)');
% Initialize physical angles to be calculated
flowAnglePlane = zeros([1,N]);
flowAngleVec = zeros([1,N]);
camAngle = zeros([1,N]);

for ii=1:N
    % If quaternions are nan, return nan
    if any(isnan(quats(ii,:)))
        flowAnglePlane(ii) = nan;
        flowAngleVec(ii) = nan;
        camAngle(ii) = nan;
        continue
    end
    % Rotate
    q = quaternion(quats(ii,:));
    p = quaternion([[0;0;0;0;0;0], PolyOriginal.A]);
    [~,x,y,z] = parts(q*p*q');
    Arot = [x,y,z];
    Poly = Polyhedron(Arot,PolyOriginal.b);
    % In case of singularity
    while length(Poly.V)<8
        quatsBla = quats(ii,:) + sign(rand([1,4])-0.5).*[10e-7, 10e-7, 10e-7, 10e-7];
        q = quaternion(quatsBla);
        [~,x,y,z] = parts(q*p*q');
        Arot = [x,y,z];
        Poly = Polyhedron(Arot,PolyOriginal.b);
    end
    % Plot oriented polygon 
    if plotFlag
        figure; hold on;

        Poly.plot("color",[0.5,0.5,0.5])
        xlabel("x"); ylabel("y"); zlabel("z");
        view(3)
        axis equal;
    end
    % Vectors normal to faces of cuboid
    nL1L2 = Poly.A(5,:);
    nL1L3 = Poly.A(3,:);
    nL2L3 = Poly.A(1,:);
    % Alignment with Flow Angles
    n_xy = [0,0,1];
    if L(1)==L(2)
        % Equant plates
        flowAnglePlane(ii) = abs(acos(nL1L2*n_xy'/norm(nL1L2))*180/pi);
        % Take smallest angle of n1/ n2
        needleVec1 = nL2L3;
        needleVec1 = needleVec1./norm(needleVec1);
        % Flip to be positive in z dir
        if needleVec1(3)<0
            needleVec1 = -needleVec1;
        end
        % Take smallest angle of n1/ n2
        needleVec2 = nL1L3;
        needleVec2 = needleVec2./norm(needleVec2);
        % Flip to be positive in z dir
        if needleVec2(3)<0
            needleVec2 = -needleVec2;
        end
        % Two equivalent angles, because L1 == L2
        theta1 = acos(needleVec1*n_xy'/norm(needleVec1))*180/pi;
        theta2 = acos(needleVec2*n_xy'/norm(needleVec2))*180/pi;
        % Choose one of them randomly
        bla = [theta1,theta2];
        ix = (rand<0.5)+1;
        flowAngleVec(ii) = bla(ix);
        if theta1 < theta2
            needleVec = needleVec1;
        else
            needleVec = needleVec2;
        end
            
    elseif L(2)==L(3)
        % Needles 
        needleVec = nL2L3;
        needleVec = needleVec./norm(needleVec);
        % Flip to be positive in z dir
        if needleVec(3)<0
            needleVec = -needleVec;
        end
        flowAngleVec(ii) = acos(needleVec*n_xy'/norm(needleVec))*180/pi;

        % Two equivalent angles, because L1 == L2
        flowAnglePlane1b = abs(acos(nL1L2*n_xy'/norm(nL1L2))*180/pi);
        flowAnglePlane2b = abs(acos(nL1L3*n_xy'/norm(nL1L3))*180/pi);
        % Choose one of them randomly
        bla = [flowAnglePlane1b,flowAnglePlane2b];
        ix = (rand<0.5)+1;
        flowAnglePlane(ii) = bla(ix);
    else
        % Unequant plates
        needleVec = nL2L3;
        needleVec = needleVec./norm(needleVec);
        % Flip to be positive in z dir
        if needleVec(3)<0
            needleVec = -needleVec;
        end
        flowAngleVec(ii) = acos(needleVec*n_xy'/norm(needleVec))*180/pi;
        flowAnglePlane(ii) = abs(acos(nL1L2*n_xy'/norm(nL1L2))*180/pi);

    end

    % Alginment with camera imaging planes
    if L(2)~=L(3)
        xyProj = [nL1L2(1:2) 0];
    else
        % Randomly pick one of the two equivalent faces
        r = rand;
        if r > 0.5
            xyProj = [nL1L2(1:2) 0];
        else
            xyProj = [nL1L3(1:2) 0];
        end
    end
    xyProj = xyProj./norm(xyProj);
    if xyProj(2)<0
        xyProj = -xyProj;
    end
    camAngle(ii) = acos(xyProj(1))*180/pi;
    
    % Plot particle wioth relevant vectors
    if plotFlag
        figure; hold on;
        % Rotate
        Poly.plot("color",[0.5,0.5,0.5])
        xlabel("x"); ylabel("y"); zlabel("z");
        view(3)
        axis equal;
        fill3(L3rectangle(:,1),L3rectangle(:,2),L3rectangle(:,3),[0,0,1],"facealpha", 0.4)
        % Camera planes
        dist = -200;
        hw = 200;
        fill3([-1 -1 1 1 -1].*hw,[dist dist dist dist dist],[-1 1 1 -1 -1].*hw,[0.3,0.3,0.3],"facealpha", 0.4)
        fill3([dist dist dist dist dist],[-1 -1 1 1 -1].*hw,[-1 1 1 -1 -1].*hw,[0.3,0.3,0.3],"facealpha", 0.4)
    end
end
physicalAngles.flowAngleVec = flowAngleVec;
physicalAngles.flowAnglePlane = flowAnglePlane;
physicalAngles.camAngle = camAngle;

if isstring(savePath)
    save(savePath, "physicalAngles")
end
end