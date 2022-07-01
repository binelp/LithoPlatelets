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
% Simulates contours of a cuboid at many random orientations, and then estimates
% its orientation based only on the simulated contours
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all;
% Get git commit ID (reproducibility)
[status,cmdout] = system('git rev-parse HEAD');
if status == 0
    % Command was successful
    commitId = cmdout(1:7);
else
    commitId = [];
end

% Options for generating projections
options.calcMoments = false;
options.thirdCam = false;
options.collectData = false;
options.instruments = {'DP'}; % Dual Projection
options.crystalHeteronym = "platelet";
options.saveFlag = false; % Saving is done outside of function

options.doPlots = 0; % Plotting
options.parforProfile = false; % Parallelization
options.parforNumCores  = 4; % Cores to use in parallel
options.nsample = 15000; % Number of samples
options.nppimg = 1; % Particles per Image
options.SEsizes =  [10,15]; % Structural element sizes for erosion and blurring

% Angle optimization options
optOpts.noScreens = 2; % 0 = no sampling, 1 = 4D icosahedron (75 points), 2 = 4D dodecahedron (330 points), 3 = 4D icosahedron subdivided (471 points), 4 =  4D icosahedron subdivided twice (3111 points)
optOpts.plotFlag = 0;
optOpts.moveFlag = 1;
optOpts.scaleExpFlag = 1;
optOpts.circLimit = 1.01; % circularity limit for excluding bubbles
optOpts.convLimit = 0.5; % convexity limit for excluding badly threshholded platelets
optOpts.allowScaling = 1; % allow scaling during the optimization (helps 45deg angle bias)

options.optOpts = optOpts;
% Lithoplatelet lengths
Ls = {[300,100,100],[100,100,30],[300,100,30],[200,100,66],[200,100,100],[300,100,66],[200,100,30]};

% Reproducibility
rng(0)
% Generate Quaternions for all runs (better comparision between different
% settings / lengths (can compare exactly same rotations))
realQuats = uniformSampledQuats(options.nsample);

simTime = datestr(now,'ddmmyyyy_HHMMss');

% Loop through different Lithoplatelet populations
for ii = 1:3
    fprintf('\n%s -> Simulating lengths %i of %i.\n',datestr(now),ii,length(Ls));
    tic
    [out, outOptions] = VTB_main_rotation(Ls{ii}, realQuats, options);
    elapsed = toc;
    % Save results
    outputStruct.timeElapsed = toc;
    outputStruct.simOptions = options;
    outputStruct.vtbOptions = outOptions;
    outputStruct.results = out;
    outputStruct.L = Ls{ii};
    outputStruct.realQuats = realQuats;

    simulationFileName = join(['SimulationResults_s',num2str(jj),'_L',num2str(ii),'_',simTime,'_',commitId],"");
    simulationFilePath = join(['D:/SimResultsCuboidAngles/',filesep,simulationFileName],"");
    save(simulationFilePath,'outputStruct');
end
