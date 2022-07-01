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
% Estimates orientation of a particle based on its contours and lengths. 
%
% Input arguments:
% - A, B:        Arrays containing points of contour for each camera
% - L:          Particle lengths
% - HM:         Cell containing H and M matrix describing particle polytope
% - options:    Struct containing options (plotFlag for plotting, moveFlag 
%               to move contours to same height, scaleExpFlag to scale 
%               experimental contours to same size, scaleFlag to allow for
%               scaling during optimization (helps with blur issues))
%
% Output arguments:
% - quats:     Quaternion array describing particle orientation. Size(quats)
%               = [1 x 4] 
% - fOpt:      Optimal cost function value (can indicate how well the 
%              orientation matches contours)
% - focusDist: Distance of object from focal plane
% - lz:        Length of contours in z/flow direction
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [quats, fOpt, focusDist, lz] = getAngleCuboidQuat(A, B, L, HM, options)
    id ='MATLAB:polyshape:repairedBySimplify';
    warning('off',id);
    
    % Different convention than in technique_IA, z last here
    A = A(:,[2,1]);
    B = B(:,[2,1]);
    % Get position of particle
    xzPos = min(A) + (max(A)-min(A))/2;
    yzPos = min(B) + (max(B)-min(B))/2;
    % Get distance from focus plane (located at 1000um)
    focusDistA = (xzPos(1) - 1000);
    focusDistB = (yzPos(1) - 1000);
    focusDist = [focusDistA focusDistB];
    % Move contours to origin to simplify rotations
    A = A - xzPos;
    B = B - yzPos;
    % Length (along z-axis) of contours (should be same in both)
    lA = max(A(:,2))-min(A(:,2));
    lB = max(B(:,2))-min(B(:,2));
    lz = [lA lB];
    lDiff = lB - lA;

    % Scale larger contour to the size of the smaller one
    if options.scaleExpFlag
        if lDiff < 0
            % lA > lB
            A(:,2) = A(:,2).*lB/lA;
        else
            % lA < lB
            B(:,2) = B(:,2).*lA/lB;
        end
    end
    % Get areas of both contours
    xzConvHullVert = convhull(A(:,1),A(:,2));
    xzOriginal = polyshape(A(xzConvHullVert(1:end-1),:));
    xzOriginalArea = area (xzOriginal);
    
    yzConvHullVert = convhull(B(:,1),B(:,2));
    yzOriginal = polyshape(B(yzConvHullVert(1:end-1),:));
    yzOriginalArea = area (yzOriginal);
    
    % Ignore bubbles
    xzCircularity = (perimeter (xzOriginal)^ 2) / (4 * pi * xzOriginalArea );
    yzCircularity = (perimeter (yzOriginal)^ 2) / (4 * pi * yzOriginalArea );
    % Ignore nonconvex particles (empty inside due to threshholding)
    xzConvexity = area(polyshape(A))/xzOriginalArea;
    yzConvexity = area(polyshape(B))/yzOriginalArea;
    if (xzCircularity < options.circLimit && yzCircularity < options.circLimit ) || (xzConvexity < options.convLimit || yzConvexity < options.convLimit)
       quats = [nan nan nan nan];
       fOpt = nan;
       lz = nan;
       focusDist = nan;
       return
    end
    %% Fit particle
    % Make cuboid polyhedron object with same dimensions as real lithoplatelet
    P = Polyhedron(HM{1},HM{2}*L');
    % Get initial guess
    if options.noScreens > 0
        quats0 = screening(options.noScreens, P, xzOriginal, xzOriginalArea, yzOriginal, yzOriginalArea);
    else
        ax = [1,1,1];
        ax = ax./norm(ax);
        angle = pi/4;
        quats0 = [cos(angle) sin(angle).*ax];
        quats0 = quats0./vecnorm(quats0,2,2);
    end
    % Plot initial guess of orientation
    if options.plotFlag
        dist = 200;
        figure;
        set(gcf, 'Units', 'Inches', 'Position', [3, 3, 3.5, 3]);
        hold on;
        % Rotate cuboid with initial orientation quaternion
        q = quaternion(quats0);
        p = quaternion([[0;0;0;0;0;0], P.A]);
        [~,x,y,z] = parts(q*p*q');
        Arot = [x,y,z];
        P0 = Polyhedron(Arot,P.b);

        colors = {[0, 83, 143]./255, [34, 245, 179]./255};
        % Plot contours from image
        plotContours(A, B, yzConvHullVert, xzConvHullVert, colors{1}, dist)
        % Plot cuboid and it's projections
        plotPwithFakeContours(P0, colors{2}, 0.4, dist - 3)
        xlabel("x"); ylabel("y"); zlabel("z");
        axis equal;box on; grid off;
        ax = gca ;
        ax.XTick = [];
        ax.YTick = [];
        ax.ZTick = [];
        title("Initial Fit")
    end
    % Define cost function for optimization
    fun = @(quats)getLonelyAreas(quats, P, xzOriginal, xzOriginalArea, yzOriginal, yzOriginalArea, options.allowScaling);
         
    % Define lower and upper bounds of quaternion for optimization
    lb = [-1,-1,-1,-1];
    ub = [1,1,1,1];
    % If plotting activated, also print information about optimization,
    % otherwise don't 
    if options.plotFlag
        fminconOptions = optimoptions('fmincon','Display','iter','Algorithm','sqp');
    else
        fminconOptions = optimoptions('fmincon','Algorithm','sqp','Display','off');
    end
    % Optimization
    [quats,fOpt] = fmincon(fun, quats0,[],[],[],[],lb,ub,[],fminconOptions);

    % Normalize quaternion (to maintain size)
    quats = quats./vecnorm(quats,2,2);
    % Plot final optimized orientation
    if options.plotFlag
        figure
        dist = 150;
        f=0.7;
        set(gcf, 'Units', 'Inches', 'Position', [3, 3, 3.5*f, 3*f]);
        hold on;
        % Rotate cuboid with final orientation quaternion
        q = quaternion(quats);
        p = quaternion([[0;0;0;0;0;0], P.A]);
        [~,x,y,z] = parts(q*p*q');
        Arot = [x,y,z];
        Popt = Polyhedron(Arot,P.b);
        % Plot contours from image
        plotContours(A, B, yzConvHullVert, xzConvHullVert,colors{1}, dist)
        % Plot cuboid and it's projections
        plotPwithFakeContours(Popt, colors{2}, 0.4, dist - 3)
        xlabel("x");ylabel("y");zlabel("z");
        view(3)
        axis equal;box on; grid off;
        ax = gca ;
        ax.XTick = [];
        ax.YTick = [];
        ax.ZTick = [];
        title("Final Fit")
    end
end
% Function for screening various orientations to find a good initial
% orientation
function [quats0] = screening(n, P, xzOriginal, xzOriginalArea, yzOriginal, yzOriginalArea)
if n == 1
    [vertices] = cell600vertices(0);
elseif n == 2
    [vertices] = cell120vertices();
elseif n == 3
    load("cell600_subdiv1.mat");
    vertices = vert;
elseif n == 4
    load("cell600_subdiv2.mat");
    vertices = vert;
end
% Remove equivalent quaternions (-q same as q)
ixNeg = vertices(:,1)<0;
vertices(ixNeg,:) = -vertices(ixNeg,:);
vertices = unique(vertices,"rows");
% Add a random rotation to everything to not bias optimization to same
% inital guesses
qref = quaternion(uniformSampledQuats(1));
vertices = qref.*quaternion(vertices(:,:));
[w,x,y,z] = parts(vertices);
pts = [w,x,y,z];
% Evaluate all options
metric = zeros([length(pts),1]);
for ii=1:length(pts)
    quats0 = pts(ii,:)./norm(pts(ii,:));
    metric(ii) = getLonelyAreas(quats0, P, xzOriginal, xzOriginalArea, yzOriginal, yzOriginalArea, 0);
end
% Pick best one 
minVal = min(metric,[],"all");
ix = find(metric==minVal);
if length(ix)>1
    ix = ix(1);
end
quats0 = pts(ix,:);
end
% Function to calculate non-overlapping areas between real contours and
% projections of cuboid, to evaluate goodness of fit
function [lonelyAreas] = getLonelyAreas(quats, P, xzOriginal, xzOriginalArea, yzOriginal, yzOriginalArea, allowScaling)
% Normalize quaternion to avoid scaling
if ~allowScaling
    quats = quats./vecnorm(quats,2,2);
end
% Rotate
q = quaternion(quats);
p = quaternion([[0;0;0;0;0;0], P.A]);

[~,x,y,z] = parts(q*p*q');
Arot = [x,y,z];
Px = Polyhedron(Arot,P.b);

% Compute projections 
try
    yz = convhull(Px.V(:,2),Px.V(:,3));
catch
    % In case of failure, set very high cost function value so optimization doesnt settle here 
    lonelyAreas = 10e10;
    return
end
try
    xz = convhull(Px.V(:,1),Px.V(:,3));
catch
    % In case of failure, set very high cost function value so optimization doesnt settle here 
    lonelyAreas = 10e10;
    return
end
% Get intersecting area between real and fake contours for xz contour
xzCurrent = polyshape(Px.V(xz(1:end-1),[1,3]));
xzCurrentArea = area(xzCurrent);
polyout = intersect(xzOriginal, xzCurrent);
xzIntersectArea = area(polyout);
% Calculate non-intersecting area (normalized)
xzLonelyAreas = (xzCurrentArea + xzOriginalArea- 2*xzIntersectArea)/xzOriginalArea;

% Get intersecting area between real and fake contours for yz contour
yzCurrent = polyshape(Px.V(yz(1:end-1),[2,3]));
yzCurrentArea = area(yzCurrent);
polyout = intersect(yzOriginal, yzCurrent);
yzIntersectArea = area(polyout);
% Calculate non-intersecting area (normalized)
yzLonelyAreas = (yzCurrentArea + yzOriginalArea - 2*yzIntersectArea)/yzOriginalArea;

% Total non-intersecting area for both contours
lonelyAreas = xzLonelyAreas + yzLonelyAreas;
end
function [] = plotContours(A, B, yzConvHullVert, xzConvHullVert, colors, dist)
    fill3(A(xzConvHullVert,1),dist*ones(size(xzConvHullVert,1),1),A(xzConvHullVert,2),colors,"facealpha", 0.8);
    fill3(dist*ones(size(yzConvHullVert,1),1),B(yzConvHullVert,1),B(yzConvHullVert,2),colors,"facealpha", 0.8);
end