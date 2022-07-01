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
% Goes through all images of all lithoplatelet populations obtained by the DISCO
% and tries to estimate the orientation of the particles, then saves that information. 
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
% Make useful saving file name
simTime = datestr(now,'ddmmyyyy_HHMMss');
saveName = join(['resultQuatsRealSize_',simTime,'_',commitId],"");

% Load cuboidal morphology
load("CrystalData.mat");
HM = {Crystals(15).H,Crystals(15).M};
% Path to folder containing different lithoplatelet populations and corresponding lengths 
mainDir ="";
% Lithoplatelet lengths
Ls = {[205,105,66],[205,105,105],[306,106,66],[306,110,106],[101,101,31],[202,102,32],[300,102,32]};
% Directory names for corresponding lithoplatelet population
dirs = {"20210812_200x100x66_TH33","20210812_200x100x100_TH33","20210812_300x100x66_TH33",...
    "20210812_300x100x100_TH33","20210903_100x100x30_TH33","20210903_200x100x30_TH33","20210903_300x100x30_TH33"};

% Angle optimization options
options.noScreens = 2; % 0 = no sampling, 1 = 4D icosahedron (75 points), 2 = 4D dodecahedron (330 points), 3 = 4D icosahedron subdivided (471 points), 4 =  4D icosahedron subdivided twice (3111 points)
options.plotFlag = 0;
options.scaleExpFlag = 1; % scale experimental contours to even length
options.allowScaling = 1; % allow scaling during the optimization (helps 45deg angle bias)
options.circLimit = 1.01; % circularity limit for excluding bubbles
options.convLimit = 0.5; % convexity limit for excluding badly threshholded platelets

% Loop through different lithoplatelet populations
for ii = 3
	fprintf('\n%s -> Population %i with L = [%i %i %i].\n',datestr(now),ii,Ls{ii}(1),Ls{ii}(2),Ls{ii}(3));
    L = Ls{ii}./2;
    subDir = join([mainDir,dirs{ii}],"");

    % Get a list of all files and folders in this folder.
    files = dir(subDir);
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subsubFolders = files(dirFlags);
    subsubFolders = subsubFolders(3:end);

    allAngles = {};
    allfOpts = {};
    noP = {};
    lzs = {};
    focusDists = {};
    clear tempSave;
    tempSave(length(subsubFolders)) = struct();
    parfor k = 1: length(subsubFolders)
        fprintf('\n%s -> Subfolder %i of %i: %s\n', datestr(now), k, length(subsubFolders), subsubFolders(k).name);
        subFolderPath = join([subDir,"/",subsubFolders(k).name],"");
        tic
        % Estimate Angles of all particles in subfloder
        [allAngles{k}, allfOpts{k}, noP{k}, focusDists{k}, lzs{k}] = angleEstimReal(subFolderPath, HM, L, options);
        elapsed = toc;
        fprintf('         Measured %i Particles in %d.2f minutes. sec/particle: %d.2f \n', noP{k}, elapsed/60, elapsed/noP{k});
        
        % Save subsubfolder results in case of crash
        tempSave(k).angles = allAngles{k};
        tempSave(k).fOpts = allfOpts{k};
        tempSave(k).path = subFolderPath;
        tempSave(k).noP = noP{k};
        tempSave(k).focusDists = focusDists{k};
        tempSave(k).lzs = lzs{k};
        saveParfor(join([subFolderPath,"/tempSave.mat"],""),tempSave(k));
    end
    % Save final results
    results.allAngles = allAngles;
    results.fOpts = allfOpts;
    results.path = subDir;
    results.noP = noP;
    results.lzs = lzs;
    results.focusDists = focusDists;
    results.L = Ls{ii};
    saveParfor(join([subDir,"/",saveName,".mat"],""),results);
end

function []=saveParfor(name,var)
save(name,'var')
end