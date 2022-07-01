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
% Loads and pre-processes images taken by the DISCO for the purpose of
% estimating particle orientations based on their contour. 
%
% Input arguments:
% - myDir:      Path to folder containing "data_camA" and "data_camB"
%               subfolders.
% - HM:         Cell containing H and M matrix describing particle polytope
% - L:          Particle lengths
% - options:    Struct containing options (plotFlag, moveFlag, scaleFlag)
%
% Output arguments:
% - allQuats:   [n x 4] array containing quaternions for each particle, describing 
%               their orientation. 
% - fOpts:      [n x 1] array containing optimal cost function value for each 
%               particle (can indicate how well the orientation matches contours)
% - noP:        Number of particles 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [allQuats, fOpts, noP, focusDists, lzs] = angleEstimReal(myDir, HM, L, options)
%% Load Contours
myFilesA = dir(join([myDir,filesep,'data_camA',filesep,'A*.csv'],""));
myFilesB = dir(join([myDir,filesep,'data_camB',filesep,'B*.csv'],""));

% Parse out time information from the *.csv filenames
timeInfoA = myFilesA(1).name(3:end-6); timeInfoB = myFilesB(1).name(3:end-6);

% Convert filenames into datetime class objects (case separtation)
timeDataA = datetime(timeInfoA,'inputformat','ddMMyyyy-HHmmss');

% Parse out file counter of camera A from the *.csv filenames
for i = 1:length(myFilesA)
    fileCounterA(i,1) = str2double(myFilesA(i).name(19:end-4));
end

% Parse out file counter of camera B from the *.csv filenames
for i = 1:length(myFilesB)
    fileCounterB(i,1) = str2double(myFilesB(i).name(19:end-4));
end

% Perform time matching
counterDifference = zeros(length(fileCounterA),1); % memory pre-allocation
Imatch(:,1) = 1:length(fileCounterA); 
for i = 1:length(fileCounterA)
    [counterDifference(i), Imatch(i,2)] = min(abs(fileCounterA(i)-fileCounterB)); % find closest match
end
Imatch(counterDifference>0,:) = []; % if counter difference is larger than zero then discard

nMatchingFiles = size(Imatch,1); % number of matching files after time matching

DeltaZthreshold = 35;
cutOffThreshold = [-DeltaZthreshold DeltaZthreshold -0.6 0.6 -5 5 0]; % [um]

allQuats = zeros([0,4]);
fOpts = zeros([0]);

focusDists = zeros([0,2]);
lzs = zeros([0,2]);

noP = 0;
% Loop over all time-matched contours from both the cameras
for i=1:nMatchingFiles
    %% Image analysis routine
    % Generate a list of *.csv files to be read
    contFileA = join([myDir,filesep,'data_camA',filesep,'A-',timeInfoA,'-',num2str(fileCounterA(i)),'.csv'],"");
    contFileB = join([myDir,filesep,'data_camB',filesep,'B-',timeInfoB,'-',num2str(fileCounterB(i)),'.csv'],"");

    % Read contour data of images captured by cameras A/B
    [A, tbA] = readContourFile(contFileA);
    [B, tbB] = readContourFile(contFileB);
    
    % Match the particles observed by both the cameras based on their
    % z-coordinates and their centroids
    [matches, numberOfPairs] = matchProjections(tbA,tbB,cutOffThreshold);
    % Go through all pairs and estimate orientations for each
    for n=1:numberOfPairs
        noP = noP + 1;
        contourA = A{matches(n,1)}(:,1:2); contourB = B{matches(n,2)}(:,1:2);
        [quats, fOpt, focusDist, lz] = getAngleCuboidQuat(contourA, contourB, L, HM, options);
        allQuats(end+1,:) = quats;
        fOpts(end+1)=fOpt;
        focusDists(end+1,:) = focusDist;
        lzs(end+1,:)= lz;
    end
end
end