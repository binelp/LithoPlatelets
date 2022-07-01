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
% Post-process and plot the orientation distributions obtained from
% lithoplatelet images. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all;

% ID of saved results to be used
runID = "";
% Directory where results were saved
mainDir ="";
% Directories of images 
dirs = ["20210812_200x100x66_TH33","20210812_200x100x100_TH33","20210812_300x100x66_TH33",...
    "20210812_300x100x100_TH33","20210903_100x100x30_TH33","20210903_200x100x30_TH33","20210903_300x100x30_TH33"];
% Load cuboid geometry information
load("CrystalData.mat");
HM = {Crystals(15).H,Crystals(15).M};
physAnglesCell = {};

% Colors for lithoplatelets
colors = { [11 171 224]./255,...
[0 114 201]./255,...
[11 51 212]./255,...
[214 100 11]./255,...
[227 40 11]./255,...
[148 227 79]./255,...
[84 214 75]./255};
% Parameters
removalThreshhold = 1.5; 
for kk = 1:7
    % Load results
    load(join([mainDir,dirs(kk),"/resultQuatsRealSize",runID,".mat"],""))
    results = var;
    quats = vertcat(results.allAngles{:});
    fOpts = [results.fOpts{:}];
    N = sum([results.noP{:}]);
    L = results.L;
    
    %%%%%%%%%%%%%%%%%%%%%%    PREPROCESSING    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot fOpts & removal threshhold
    figure; hold on;
    set(gcf, 'Units', 'Inches', 'Position', [3, 3, 6, 6/3]);
    title(sprintf("L = [%i %i %i]",L(1),L(2),L(3)));
    numBins = 100;
    binLims = [0,2.5];
    transparency=0.4;
    hold on; grid on; box on;
    histogram(fOpts,'Facecolor',colors{kk},'Facealpha',transparency,'edgealpha',0,'NumBins',numBins,'BinLimits',binLims);
    histogram(fOpts,'Edgecolor',colors{kk},'Facealpha',0,'edgealpha',1,'NumBins',numBins, 'DisplayStyle', 'stairs','linewidth',1.5,'BinLimits',binLims);

    % Remove worst fits
    okIx = fOpts < removalThreshhold;
    fOpts = fOpts(okIx);
    quats = quats(okIx,:);
    N = length(fOpts);
   
    disp(join([sprintf("L = [%i %i %i]",L(1),L(2),L(3)), " - Percentage removed [%] : ", 100*(1-sum(okIx)/length(okIx))]))
    
    % Set camera angle to be uniform
    savePath = join([mainDir,dirs(kk),"\camAngleCorr",runID,".mat"],"");
    try
        load(savePath)
    catch
        [quats] = setCamAngleUniform(quats);
        save(savePath, "quats")
    end
    
    % Reduce equivalent quaternions
    savePath = join([mainDir,dirs(kk),"\reducedEquiv",runID,".mat"],"");
    try
        load(savePath)
    catch
        noRot = [zeros([length(quats),1]), ones([length(quats),1]), zeros([length(quats),2])];
        [~, ~, quats] = minAngleBetweenOrientations(L, quats, noRot, 1, 0);
        save(savePath, "quats")
    end
    
    %%%%%%%%%%%%%%%%%%%% ANGLE DISTRIBUTION PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Reduce number of dimensions
    quatsRed = quat4Dto3D(quats);
    % Reconstruct reduced quaternion distribution
    lims =  [-0.9,0.9;-0.9,0.9;-0.9,0.9];
    nbins = [20 20 20];
    [QD] = angleListToDist(quatsRed, lims, nbins, 1);
    % Plot reduced quaternion distribution
    figure; hold on;
    suptitle(sprintf("L = [%i %i %i]",L(1),L(2),L(3)));
    subplot(1,2,1)
    f = 0.8;
    set(gcf, 'Units', 'Inches', 'Position', [3, 3, 3.5*f, 3*f]);
    plotAD(QD,[1 0 0], 0,0.4,2,8)
    % Mirror distribution 
    quatsRed(quatsRed(:,2)<0,:) = [-quatsRed(quatsRed(:,2)<0,1) -quatsRed(quatsRed(:,2)<0,2) quatsRed(quatsRed(:,2)<0,3)];
    
    % Reconstruct reduced quaternion distribution
    lims =  [-0.9,0.9;-0,0.9;-0.9,0.9];
    nbins = [20 10 20];
    [QD] = angleListToDist(quatsRed, lims, nbins, 1);
    % Plot reduced quaternion distribution
    subplot(1,2,2)
    figure; hold on;
    set(gcf, 'Units', 'Inches', 'Position', [3, 3, 3.5*f, 3*f]);
    plotAD(QD,[1 0 0], 0,0.4,2,8)
    
    % Save quaternion distribution
    savePath = join([mainDir,dirs(kk),"\QD_",runID,".mat"],"");
    save(savePath, "QD")

    %%%%%%%%%%%%%%%%%%%%%%% PHYSICAL ANGLE DIST %%%%%%%%%%%%%%%%%%%%%%%%%%%
    savePath = join([mainDir,dirs(kk),"\physicalAngles_",runID,".mat"],"");

    % Calculate and save physical angles
    [physicalAngles] = getPhysicalAngles(quats, L, savePath);
    physAnglesCell{kk} = physicalAngles;
end