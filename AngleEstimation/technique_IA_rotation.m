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
% Simluates contours of virtual particles and tries to estimate the particles 
% orientation based on contours only.(Adapted from Shape/VTB/technique_IA.m, 
% branch: platelets, commitID: 1fca5a2)
%
% Input arguments:
% - P:          Polytopes
% - realQuats:  Quaternions representing orientation of particles 
% - options:    Options for simulating particles (i.e. blur) and for angle
%               estimation
%
% Output arguments:
% - out:        Struct containing results, i.e. real and measured
%               quaternions, and also options for reproducibility
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [out] = technique_IA_rotation(P, realQuats, options)
% Place, rotate and shift particles in the virtual flow cell until there are no more
% overlappings. Then, project, extract the contours and pass them to the KIDDO IA code.

% Settings
% -> Minimum number of pixels within contour for it to be processed
%This was changed from 50 (default in DISCO) to 0 because of issues when
%testing extreme aspect ratios
options.MIN_NUM_CNT_PX = 0;
% If not yet contained in the options structure, assign default values for cellDim,
% nppimg, cutOffThreshold, and magnification
options = check_IA_input(options);
% Initialize particle list object for both DP and SP

% Compute cell size in pixels [px]
ftcSize = ceil(options.cellDim(:)./options.magnification);
% Check that options.nsample is an integer multiple of options.nppimg
if ~mod(options.nsample,options.nppimg)
    nframes = options.nsample/options.nppimg;
else
    error('technique_IA:invalidSampleNumber',...
        'options.nsample needs to be an integer multiple of options.nppimg');
end

% allRealQuats = cell(nframes,1);
allMeasQuats = cell(nframes,1);
allfOpts = cell(nframes,1);

Pcell = cell(1,nframes);
for i=1:nframes
    Pcell{i} = P((i-1)*options.nppimg+1:1:i*options.nppimg);
end
if isstring(options.parforProfile)
    if strcmpi(options.parforProfile,'local')
        % Prepare cluster
        curPool = gcp('nocreate');
        % Acquire the required number of workers
        if (isempty(curPool) || (curPool.NumWorkers ~= options.parforNumCores))
            if ~isempty(curPool)
                % Delete old pool
                delete(curPool);
            end
            % Write status information to command window
            disp([datestr(now),' --> Preparing ',num2str(options.parforNumCores),...
                ' CPU cores for parallel operation. This initialization might ',...
                'require up to 20 seconds.']);
            % Initialize pool with the correct number of workers
            clusta = parcluster('local');
            parpool(clusta,options.parforNumCores);
        end
    elseif strcmpi(options.parforProfile,'eulerLSF')
        % Get an object of the EulerLSF8h cluster profile
        cluster = parcluster('EulerLSF8h');
        % Specify the time limit in hours and the requested memory in MB
        cluster.SubmitArguments = ['-W ',num2str(options.parforLSFTimeLimit),':00',...
            ' -R "rusage[mem=',num2str(options.parforMemoryRequest),']"'];
        % Open the MATLAB pool with the desired number of workers
        curPool = gcp('nocreate');
        tryCounter = 0;
        while (isempty(curPool) && tryCounter<=6)
            try
                tryCounter = tryCounter + 1;
                disp([datestr(now),' --> Preparing ',num2str(options.parforNumCores),...
                    ' CPU cores for parallel operation. This initialization might ',...
                    'require up to 20 seconds.']);
                parpool(cluster,options.parforNumCores);
            catch
                disp([datestr(now),' --> Allocating the workers failed. Will retry no. ',...
                    num2str(tryCounter),' of 7...']);
            end
            curPool = gcp('nocreate');
        end
    else
        error('technique_IA:invalidParforProfile','options.parforprofile is invalid.');
    end
    % Loop over the images: in the loop body, the counter firstid=firstidValues(idx)
    % always points to the first particle of the current image
    parfor idx = 1:nframes
        % Call loop body function: place particles within virtual FTC, project them,
        % compute contours of projection, call KIDDO IA code, and return measured particle
        % dimensions in iaOutput
        [measQuats, fOpts] =...
            loopBody(idx,ftcSize,Pcell, realQuats(idx,:), options);
%         allRealQuats{idx}=realQuats;
        allMeasQuats{idx}=measQuats;
        allfOpts{idx}=fOpts;
    end
else
    % Standard for-loop instead of parfor-loop
    for idx = 1:nframes
        % Call loop body function: place particles within virtual FTC, project them,
        % compute contours of projection, call KIDDO IA code, and return measured particle
        % dimensions in iaOutput
        [measQuats, fOpts] =...
            loopBody(idx,ftcSize,Pcell, realQuats(idx,:), options);
%         allRealQuats{idx}=realQuats;
        allMeasQuats{idx}=measQuats;
        allfOpts{idx}=fOpts;
    end
end

% out.realQuats = cell2mat(allRealQuats);
out.measQuats = cell2mat(allMeasQuats);
out.fOpts = cell2mat(allfOpts);

end

function [measQuats,fOpts] =loopBody(idx,ftcSize,Pcell,realQuat, options)
% for- or parfor-loop body for generating options.nppimg particles per frame, extracting
% there contours and characterizing them using the KIDDO IA pipeline

imgxz = false(ftcSize(3),ftcSize(1));
imgyz = false(ftcSize(3),ftcSize(2));
img45z = false(ftcSize(3),round(ftcSize(1)*sqrt(2)));

% BELOW THIS LINE: EVERYTHING HAS THE UNIT OF [um]

% Now, place the correct number of particles per image (options.nppimg) within the
% virtual flow cell. These should be inside of the flow cell, and must not overlap with
% other particles!
subloopP = Pcell{idx};
L = subloopP.b([1,3,5])';

measQuats = zeros([0,4]);
fOpts = zeros([0,1]);
for curid = 1:length(subloopP)
    % Initialize booleans
    inside = false;
    % Continue random rotation (R) and translation (vtrans) until
    % particle is inside the cell and does not clip=overlap with other particles
    while ~inside
        % Random translation vector, bounded by cell dimensions
        vtrans = rand([3,1]).*sign(rand([3,1])-0.5).*options.cellDim(:)./2;

        % Rotation with quaternions
        q = quaternion(realQuat);
        p = quaternion([[0;0;0;0;0;0], subloopP(curid).A]);

        [~,x,y,z] = parts(q*p*q');
        Arot = [x,y,z];
        curP = Polyhedron(Arot,subloopP(curid).b + subloopP(curid).A*vtrans);

        % Get vertices of current particle
        verticesCurP = curP.V;
        
    end % while loop
    
    % Now that we have a feasible rotation and translation, assign curP to original
    % set of polytopes
    subloopP(curid) = curP;

    % Project particle on xz-plane
    xzProjVert = [verticesCurP(:,1),verticesCurP(:,3)];
    % Compute convex hull of projection. MATLAB doc: "The convex hull K is expressed
    % in terms of a vector of point indices arranged in a counterclockwise cycle
    % around the hull."
    try
        xzConvHullVert = convhull(verticesCurP(:,1),verticesCurP(:,3));
    catch
        % In some situations, the particle appears like a plate with zero
        % thickness and this will not allow the computation of the convex
        % hull of the particle, as the particle will be composed of only 2
        % point. Print a warning for such a case and continue to the next
        % iteration of the for loop
        warning on;
        warning('techniqueIA:notEnoughPointsConvexHull',...
            'Error computing the convex hull. Not enough unique points specified.');
        warning off;
        continue;
    end
    % BELOW THIS LINE: EVERYTHING HAS THE UNIT OF pixels
    
    % Translate result into pixel coordinate frame that has its origin at a
    % corner of the xz-frame
    xzProjVert = xzProjVert./options.magnification;
    xzProjVert(:,1) = xzProjVert(:,1) + ftcSize(1)/2;
    xzProjVert(:,2) = xzProjVert(:,2) + ftcSize(3)/2;
    % Get bounding box of the current xz convex hull and round to integer pixel coordinate
    xzMaxX = ceil(max(xzProjVert(xzConvHullVert,1)));
    xzMinX = floor(min(xzProjVert(xzConvHullVert,1)));
    xzMaxZ = ceil(max(xzProjVert(xzConvHullVert,2)));
    xzMinZ = floor(min(xzProjVert(xzConvHullVert,2)));
    % Set of points to be tested if inside polygon -> inpoly
    [xzMeshX,xzMeshZ] = meshgrid(xzMinX:xzMaxX,xzMinZ:xzMaxZ);
    xzPoints = [xzMeshX(:),xzMeshZ(:)];
    % Vertices of convex hull -> inpoly
    xzNodes = [xzProjVert(xzConvHullVert,1),xzProjVert(xzConvHullVert,2)];
    % inpoly(.) test
    xzIn = inpoly(xzPoints,xzNodes);
    %xzIn = inpolygon(xzPoints(:,1),xzPoints(:,2),xzNodes(:,1),xzNodes(:,2));
    % Update binary image after inpoly(.) test
    xzInsidePoints = xzPoints(xzIn,:);
    for pxCount = 1:size(xzInsidePoints,1)
        imgxz(xzInsidePoints(pxCount,2),xzInsidePoints(pxCount,1)) = true;
    end
    % Repeat everything for yz
    yzProjVert = [verticesCurP(:,2),verticesCurP(:,3)];
    try
        yzConvHullVert = convhull(verticesCurP(:,2),verticesCurP(:,3));
    catch
        % In some situations, the particle appears like a plate with zero
        % thickness and this will not allow the computation of the convex
        % hull of the particle, as the particle will be composed of only 2
        % point. Print a warning for such a case and continue to the next
        % iteration of the for loop
        warning on;
        warning('techniqueIA:notEnoughPointsConvexHull',...
            'Error computing the convex hull. Not enough unique points specified.');
        warning off;
        continue;
        
    end
    yzProjVert = yzProjVert./options.magnification;
    yzProjVert(:,1) = yzProjVert(:,1) + ftcSize(2)/2;
    yzProjVert(:,2) = yzProjVert(:,2) + ftcSize(3)/2;
    yzMaxY = ceil(max(yzProjVert(yzConvHullVert,1)));
    yzMinY = floor(min(yzProjVert(yzConvHullVert,1)));
    yzMaxZ = ceil(max(yzProjVert(yzConvHullVert,2)));
    yzMinZ = floor(min(yzProjVert(yzConvHullVert,2)));
    [yzMeshY,yzMeshZ] = meshgrid(yzMinY:yzMaxY,yzMinZ:yzMaxZ);
    yzPoints = [yzMeshY(:),yzMeshZ(:)];
    yzNodes = [yzProjVert(yzConvHullVert,1),yzProjVert(yzConvHullVert,2)];
    yzIn = inpoly(yzPoints,yzNodes);
    yzInsidePoints = yzPoints(yzIn,:);
    for pxCount = 1:size(yzInsidePoints,1)
        imgyz(yzInsidePoints(pxCount,2),yzInsidePoints(pxCount,1)) = true;
    end   

    %avoid error later on (Undefined function or variable 'contoursxz'.) 
    %which happens when tempcontours empty. Just redoes everything
    % Create contours
    tempcontoursxz = bwboundaries(imgxz,'noholes'); % Returns a cell array where each element of the cell
    tempcontoursyz = bwboundaries(imgyz,'noholes');
    
    counter = 1;
    for jj=1:length(tempcontoursxz)
        if size(tempcontoursxz{jj},1) > options.MIN_NUM_CNT_PX
            contoursxz{counter} = tempcontoursxz{jj}; %#ok<*AGROW> % [px]
            counter = counter + 1;
        end 
    end
    counter = 1;
    for jj=1:length(tempcontoursyz)
        if size(tempcontoursyz{jj},1) > options.MIN_NUM_CNT_PX
            contoursyz{counter} = tempcontoursyz{jj}; % [px]
            counter = counter + 1;
        end 
    end
end % for all polytopes in a single frame
% Plot only first frame if plots enabled
if options.doPlots && idx==1
    
    figure; 
    x0=600; y0=300; width=1800; height=400;
    set(gcf,'position',[x0,y0,width,height]);
    hold on;
    
    subplot(1,3,1); hold on;
    for kk=1:length(subloopP)
        % Plot original particle
        subloopP(kk).plot;
    end
    xlabel(['x [',char(181),'m]'],'interpreter','tex');
    ylabel(['y [',char(181),'m]'],'interpreter','tex');
    zlabel(['z [',char(181),'m]'],'interpreter','tex');
    hold on; axis equal; box off;
    axis([-options.cellDim(1)/2 options.cellDim(1)/2 -options.cellDim(2)/2 options.cellDim(2)/2 -options.cellDim(3)/2 options.cellDim(3)/2]);
    set(gcf,'color','w')
   
    subplot(1,3,2);hold on;
    for kk=1:length(subloopP)
        subloopP(kk).plot;
    end
    xlabel(['x [',char(181),'m]'],'interpreter','tex');
    ylabel(['y [',char(181),'m]'],'interpreter','tex');
    zlabel(['z [',char(181),'m]'],'interpreter','tex');
    hold on; axis equal; box on;
    axis([-options.cellDim(1)/2 options.cellDim(1)/2 -options.cellDim(2)/2 options.cellDim(2)/2 -options.cellDim(3)/2 options.cellDim(3)/2]);
    set(gcf,'color','w')
     view(90,0);
     
     subplot(1,3,3);hold on;
    for kk=1:length(subloopP)
        subloopP(kk).plot;
    end
    xlabel(['x [',char(181),'m]'],'interpreter','tex');
    ylabel(['y [',char(181),'m]'],'interpreter','tex');
    zlabel(['z [',char(181),'m]'],'interpreter','tex');
    hold on; axis equal; box on;
    axis([-options.cellDim(1)/2 options.cellDim(1)/2 -options.cellDim(2)/2 options.cellDim(2)/2 -options.cellDim(3)/2 options.cellDim(3)/2]);
    set(gcf,'color','w')
     view(0,0);
    if options.saveFlag
        export_fig(['SimulationResults/PlateletAngles_',datestr(now,'ddmmyyyy_HHMMss'),'.png'])
    end
    
    figure; hold on; imshow(flipud(imgxz)); xlabel('x'); ylabel('z');axis equal;
    figure; hold on; imshow(flipud(imgyz)); xlabel('y'); ylabel('z');axis equal;
    
    if options.thirdCam
        figure; hold on; imshow(flipud(img45z)); xlabel('x/y'); ylabel('z');axis equal;
    end
    % Delete subloopP since not needed anymore
    clear subloopP;
end

% Round contours to simulate blur
SE1 = strel("disk",options.SEsizes(1));
SE2 = strel("disk",options.SEsizes(2));
er = imerode(imgxz, SE1);
dil = imdilate(er, SE2);
imgxz = dil;

er = imerode(imgyz, SE1);
dil = imdilate(er, SE2);
imgyz = dil;

% Create contours
tempcontoursxz = bwboundaries(imgxz,'noholes'); % Returns a cell array where each element of the cell
tempcontoursyz = bwboundaries(imgyz,'noholes');

% Remove the particles where one of the contours contains less than MIN_NUM_CNT_PX pixels
counter = 1;
for jj=1:length(tempcontoursxz)
    if size(tempcontoursxz{jj},1) > options.MIN_NUM_CNT_PX
        contoursxz{counter} = tempcontoursxz{jj}; %#ok<*AGROW> % [px]
        counter = counter + 1;
    end 
end
counter = 1;
for jj=1:length(tempcontoursyz)
    if size(tempcontoursyz{jj},1) > options.MIN_NUM_CNT_PX
        contoursyz{counter} = tempcontoursyz{jj}; % [px]
        counter = counter + 1;
    end 
end

% Process the simulated particles
A = contoursxz{1};
B = contoursyz{1};

% Estimate orientation based on generated contours
[quat, fOpt] = getAngleCuboidQuat(A, B, L, {options.A,options.M}, options.optOpts);

% Append Results
measQuats(end+1,:) = quat;
fOpts(end+1) = fOpt;
end

function options = check_IA_input(options)
% If not yet contained in the options structure, assign default values for options.cellDim,
% options.nppimg, options.cutOffThreshold, and options.magnification

% Camera magnification in [um] per pixel
if ~isfield(options,'magnification')
    options.magnification = 0.9922; % new isotropic scaling factor [um/pixel]
end
% Photographed flow cell dimensions in [um]
if ~isfield(options,'cellDim') 
    options.cellDim = floor([2048, 2048, 2448].*options.magnification);
    % FOV: 2032 x 2032 x 2428 [um] in "real" world coordinates
end
% Number of particles per image
if ~isfield(options,'nppimg') 
    options.nppimg = 10;
end
% Cut off threshold
if ~isfield(options,'cutOffThreshold') 
    options.cutOffThreshold = 15; % [px] (for matching)
end

end % function
