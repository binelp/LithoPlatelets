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
% Combines all lithoplatelet orientation distributions to be interpolated
% in the VTB (simulation of length measurement).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;

%% Save file containing all angle distributions and corresponding lengths for interpolation
% ID of saved results to be used
runID = "";
% Directory where data and results stored
mainDir ="";
% Subdirectory for each population
dirs = ["20210812_200x100x66_TH33","20210812_200x100x100_TH33","20210812_300x100x66_TH33",...
    "20210812_300x100x100_TH33","20210903_100x100x30_TH33","20210903_200x100x30_TH33","20210903_300x100x30_TH33"];
% Lithoplatelet lengths
Ls = [205,105,66; 
    205,105,105; 
    306,106,66; 
    306,110,106; 
    101,101,31; 
    202,102,32; 
    300,102,32];
% Load all the individual orientation distributions from each lihtoplatelet population 
for kk=1:7
    path = join([mainDir,dirs(kk),"\QD_",runID,".mat"],"");
    load(path, "QD")
    allQDs(kk).QD = QD;
    allQDs(kk).AS = [Ls(kk,1)./Ls(kk,2) Ls(kk,2)./Ls(kk,3)];
end

% Add uniform dist at AS = [1, 1] (equant, therefore all oreintations equally likely)
allQDs(8).AS = [1 1];
% Reduce equivalent quaternions
quats = uniformSampledQuats(15000);
noRot = [zeros([length(quats),1]), ones([length(quats),1]), zeros([length(quats),2])];
[~, ~, quats] = minAngleBetweenOrientations([300 200 100], quats, noRot, 1, 0);
% Reduce number of dimensions
quatsRed = quat4Dto3D(quats);
% Mirror distribution 
quatsRed(quatsRed(:,2)<0,:) = [-quatsRed(quatsRed(:,2)<0,1) -quatsRed(quatsRed(:,2)<0,2) quatsRed(quatsRed(:,2)<0,3)];
% Reconstruct reduced quaternion distribution
lims =  [-0.9,0.9;-0,0.9;-0.9,0.9];
nbins = [20 10 20];
[QD] = angleListToDist(quatsRed, lims, nbins, 1);
allQDs(8).QD = QD;

% Interpolate QD at [1, 1.5]
AS = cat(1,allQDs.AS);
denom = AS(8,2) - AS(5,2);
I = (AS(8,2)-1.5)/denom*allQDs(5).QD.F + (1.5-AS(5,2))/denom*allQDs(8).QD.F;
allQDs(9).QD.F = I;
allQDs(9).AS = [1, 1.5];
allQDs(9).QD.grid = allQDs(1).QD.grid;

% Copy distribution to edges (maximum aspect ratio)
endAS = 5000;
ixs = [5,6,7,9,1,3,8,2,4];
allQDs = allQDs(ixs);

% Put aspect ratios and orientation into one struct to save
allQDs(10).AS = [1,endAS];
allQDs(10).QD = allQDs(1).QD;
allQDs(11).AS = [2,endAS];
allQDs(11).QD = allQDs(2).QD;
allQDs(12).AS = [3,endAS];
allQDs(12).QD = allQDs(3).QD;
allQDs(13).AS = [endAS,endAS];
allQDs(13).QD = allQDs(3).QD;
allQDs(14).AS = [endAS,3.33333];
allQDs(14).QD = allQDs(3).QD;
allQDs(15).AS = [endAS,1.5];
allQDs(15).QD = allQDs(6).QD;
allQDs(16).AS = [endAS,1];
allQDs(16).QD = allQDs(9).QD;

% Save QD struct for usage in VTB
save( join([mainDir,"allQDs_",runID,".mat"],""),"allQDs")
AS = cat(1,allQDs.AS);
QDs = {allQDs.QD};

% Plot Quaternion Distributions
fontsize = 14;
figure; 
set(gcf, 'Units', 'Inches', 'Position', [3, 0, 10, 9]);
subplot(3,3,1);
plotAD(allQDs(1).QD, [0.8,0,0], 0, 0.4, 2, fontsize );
subplot(3,3,2);
plotAD(allQDs(2).QD, [0.8,0,0], 0, 0.4, 2, fontsize );
subplot(3,3,3);
plotAD(allQDs(3).QD, [0.8,0,0], 0, 0.4, 2, fontsize );

subplot(3,3,4);
plotAD(allQDs(4).QD, [0.8,0,0], 0, 0.4, 2, fontsize );
subplot(3,3,5);
plotAD(allQDs(5).QD, [0.8,0,0], 0, 0.4, 2, fontsize );
subplot(3,3,6);
plotAD(allQDs(6).QD, [0.8,0,0], 0, 0.4, 2, fontsize );

subplot(3,3,7);
plotAD(allQDs(7).QD, [0.8,0,0], 0, 0.4, 2, fontsize );
subplot(3,3,8);
plotAD(allQDs(8).QD, [0.8,0,0], 0, 0.4, 2, fontsize );
subplot(3,3,9);
plotAD(allQDs(9).QD, [0.8,0,0], 0, 0.4, 2, fontsize );