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
% To analyze the simulations of the orientation estimation. To asses and
% compare methodologies. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all;
% Directory where results were saved
mainDir = "";
% Files containing results to consider
files = {};
% Plotting size factor and colors
sizeFactor = 2;
colors = [[0.9843 0.7294 0.4471]  [0.7922 0.3255 0.0627]  [0.7333 0.3020 0]  [0.5608 0.1451 0.0471]  [0.4118 0.1176 0.0235]];

% Loop through files
for ii = 1:length(files)
    file=files{ii};
    load(join([mainDir,file],""))

    measQuats{ii} = outputStruct.results.measQuats;
    realQuats{ii} = outputStruct.realQuats;
    
    fOpts{ii} = outputStruct.results.fOpts;
    L{ii} = outputStruct.L;
    
    % Find smallest rotation needed to go from real to estimated
    % orientation to compare the two
    [minAngles{ii}, ~, ~] = minAngleBetweenOrientations(L{ii}, measQuats{ii}, realQuats{ii}, 0, 0);
    minAngles{ii} = minAngles{ii}./pi.*180;
    length(find(minAngles{ii} < 10))/length(realQuats{ii})
    figure
    hist(minAngles{ii},50)
    
    % look at good examples
    [~,sIx] = sort(minAngles{ii});
    for jj = 1:10
        plotPartComp(L{ii},realQuats{ii}(sIx(jj),:),measQuats{ii}(sIx(jj),:))
    end
    % look at medium examples
    for jj = 1:10
        plotPartComp(L{ii},realQuats{ii}(sIx(11500+jj),:),measQuats{ii}(sIx(11500+jj),:))
    end
    % look at bad examples
    for jj = 1:10
        plotPartComp(L{ii},realQuats{ii}(sIx(end-jj),:),measQuats{ii}(sIx(end-jj),:))
    end
    
    % Reduce to 3D
    measQuatsRed = quat4Dto3D(measQuats{ii});
    realQuatsRed = quat4Dto3D(realQuats{ii});
    % Construct measured AngleDist
    lims =  [-1.25,1.25;-1.25,1.25;-1.25,1.25];
    nbins = [15 15 15];
    [ADmeas{ii}] = angleListToDist(measQuatsRed, lims, nbins, 1);
    [ADreal] = angleListToDist(realQuatsRed, lims, nbins, 1);
    
    % Plot orientation distributions
    colors = {[0, 83, 143]./255, [34, 245, 179]./255};
    figure
    f = 1;
    set(gcf, 'Units', 'Inches', 'Position', [3, 3, 3.5*f, 3*f]);
    plotAD(ADreal, colors{1}, 0,0.4,2,8)
    
    figure
    f = 1;
    set(gcf, 'Units', 'Inches', 'Position', [3, 3, 3.5*f, 3*f]);
    plotAD(ADmeas{ii}, colors{2}, 0,0.4,2,8)
end
%% Plots
realPhysAngles = realPhysAngles.physicalAngles;
physAngles{1} = physAngles{1}.physicalAngles;
physAngles{2} = physAngles{2}.physicalAngles;
physAngles{3} = physAngles{3}.physicalAngles;
% Plot physical angles
plotPhysAngles(physAngles, realPhysAngles, colors, sizeFactor);
% Plot distribution of optimized cost function values
h = plotFoptDist(fOpts, colors, sizeFactor);
%% Functions
function [] = plotPhysAngles(physAngles, realPhysAngles, colors, sizeFactor)
% Plot Physical Angle Distribution
figure; hold on;
set(gcf, 'Units', 'Inches', 'Position', [3, 0, 7*sizeFactor, 3*3.3*sizeFactor]);
numBins = 40;
binLims = [0,90];
transparency=0.4;
subplot(3,2,1);
hold on; grid on; box on;
histogram(realPhysAngles.flowAngleVec,'Facecolor',colors{1},'Facealpha',transparency,'edgealpha',0,'NumBins',numBins,'Normalization','pdf','BinLimits',binLims);
histogram(realPhysAngles.flowAngleVec,'Edgecolor',colors{1},'Facealpha',0,'edgealpha',1,'NumBins',numBins,'Normalization','pdf', 'DisplayStyle', 'stairs','linewidth',1.5,'BinLimits',binLims)
xlabel("flow angle vec")
ylim([0,0.025])
subplot(3,2,2)
hold on; grid on; box on;
for jj=1:3
h(jj)=histogram(physAngles{jj}.flowAngleVec,'Facecolor',colors{jj+1},'Facealpha',transparency,'edgealpha',0,'NumBins',numBins,'Normalization','pdf','BinLimits',binLims);
histogram(physAngles{jj}.flowAngleVec,'Edgecolor',colors{jj+1},'Facealpha',0,'edgealpha',1,'NumBins',numBins,'Normalization','pdf', 'DisplayStyle', 'stairs','linewidth',1.5,'BinLimits',binLims)
end
xlabel("flow angle vec")
ylim([0,0.025])

subplot(3,2,3);
hold on; grid on; box on;
histogram(realPhysAngles.flowAnglePlane,'Facecolor',colors{1},'Facealpha',transparency,'edgealpha',0,'NumBins',numBins,'Normalization','pdf','BinLimits',binLims);
histogram(realPhysAngles.flowAnglePlane,'Edgecolor',colors{1},'Facealpha',0,'edgealpha',1,'NumBins',numBins,'Normalization','pdf', 'DisplayStyle', 'stairs','linewidth',1.5,'BinLimits',binLims)
xlabel("flowAnglePlane")
ylim([0,0.035])
subplot(3,2,4)
hold on; grid on; box on;
for jj=1:3
h(jj)=histogram(physAngles{jj}.flowAnglePlane,'Facecolor',colors{jj+1},'Facealpha',transparency,'edgealpha',0,'NumBins',numBins,'Normalization','pdf','BinLimits',binLims);
histogram(physAngles{jj}.flowAnglePlane,'Edgecolor',colors{jj+1},'Facealpha',0,'edgealpha',1,'NumBins',numBins,'Normalization','pdf', 'DisplayStyle', 'stairs','linewidth',1.5,'BinLimits',binLims)
end
xlabel("flowAnglePlane")
ylim([0,0.035])

binLims = [0,180];
subplot(3,2,5);
hold on; grid on; box on;
histogram(realPhysAngles.camAngle,'Facecolor', colors{1},'Facealpha',transparency,'edgealpha',0,'NumBins',numBins,'Normalization','pdf','BinLimits',binLims);
histogram(realPhysAngles.camAngle,'Edgecolor', colors{1},'Facealpha',0,'edgealpha',1,'NumBins',numBins,'Normalization','pdf', 'DisplayStyle', 'stairs','linewidth',1.5,'BinLimits',binLims)
xlabel("camAngle")
ylim([0,0.015])
subplot(3,2,6)
hold on; grid on; box on;
for jj=1:3
h(jj)=histogram(physAngles{jj}.camAngle,'Facecolor',colors{jj+1},'Facealpha',transparency,'edgealpha',0,'NumBins',numBins,'Normalization','pdf','BinLimits',binLims);
histogram(physAngles{jj}.camAngle,'Edgecolor',colors{jj+1},'Facealpha',0,'edgealpha',1,'NumBins',numBins,'Normalization','pdf', 'DisplayStyle', 'stairs','linewidth',1.5,'BinLimits',binLims)
end
xlabel("camAngle")
ylim([0,0.015])


binLims = [0,90];
figure; hold on;
set(gcf, 'Units', 'Inches', 'Position', [3, 3, 6, 3.5]);
subplot(2,1,1);
hold on; grid on; box on;
histogram(realPhysAngles.camAngle,'Facecolor', colors{1},'Facealpha',transparency,'edgealpha',0,'NumBins',numBins,'Normalization','pdf','BinLimits',binLims);
histogram(realPhysAngles.camAngle,'Edgecolor', colors{1},'Facealpha',0,'edgealpha',1,'NumBins',numBins,'Normalization','pdf', 'DisplayStyle', 'stairs','linewidth',1.5,'BinLimits',binLims)
ylim([0,0.03])
ylabel("$p$", 'interpreter', 'latex')
ax = gca;
ax.FontSize =  14; 
xlim([0,90])

subplot(2,1,2); hold on; grid on; box on;
for jj=1:3
h(jj)=histogram(physAngles{jj}.camAngle,'Facecolor',colors{jj+1},'Facealpha',transparency,'edgealpha',0,'NumBins',numBins,'Normalization','pdf','BinLimits',binLims);
histogram(physAngles{jj}.camAngle,'Edgecolor',colors{jj+1},'Facealpha',0,'edgealpha',1,'NumBins',numBins,'Normalization','pdf', 'DisplayStyle', 'stairs','linewidth',1.5,'BinLimits',binLims)
end
xlabel("$\psi$ [$^\circ$]", 'interpreter', 'latex')
ylabel("$p$", 'interpreter', 'latex')
ylim([0,0.03])
xlim([0,90])
ax = gca;
ax.FontSize =  14; 

end
% Plots distributions of optimized cost function values
function [h] = plotFoptDist(fOpts, colors, sizeFactor)
% Plot fOpts Distribtuion
figure; hold on;
set(gcf, 'Units', 'Inches', 'Position', [3, 0, 7*sizeFactor, 3.3*sizeFactor]);
numBins=50;
binLims = [0,1.5];
transparency=0.4;
hold on; grid on; box on;
for jj=1:3
h(jj)=histogram(fOpts{jj},'Facecolor',colors{jj+1},'Facealpha',transparency,'edgealpha',0,'NumBins',numBins,'Normalization','pdf','BinLimits',binLims);
histogram(fOpts{jj},'Edgecolor',colors{jj+1},'Facealpha',0,'edgealpha',1,'NumBins',numBins,'Normalization','pdf', 'DisplayStyle', 'stairs','linewidth',1.5,'BinLimits',binLims)
end
xlabel("$f_{\mathrm{opt}}$",'Interpreter','Latex','Fontsize',8*sizeFactor)
% legend(h,{sprintf("L = [%i %i %i]",L{1}(1),L{1}(2),L{1}(3)),sprintf("L = [%i %i %i]",L{2}(1),L{2}(2),L{2}(3)),sprintf("L = [%i %i %i]",L{3}(1),L{3}(2),L{2}(3))},'Fontsize',8*sizeFactor)
legend(h,{"Needle","Platelet","Needle-platelet"},'Fontsize',8*sizeFactor)
% 
end
% Plot a comparision between two different orientations and the same
% particle
function [] = plotPartComp(L,quats1,quats2)
load("CrystalData.mat");
HM = {Crystals(15).H,Crystals(15).M};

PolyOriginal = Polyhedron(HM{1},HM{2}*(L./2)');
p = quaternion([[0;0;0;0;0;0], PolyOriginal.A]);

% Rotate
q1 = quaternion(quats1);
[~,x,y,z] = parts(q1*p*q1');
Arot1 = [x,y,z];
Poly1 = Polyhedron(Arot1,PolyOriginal.b);

% Rotate 
q2 = quaternion(quats2);
[~,x,y,z] = parts(q2*p*q2');
Arot2 = [x,y,z];
Poly2 = Polyhedron(Arot2,PolyOriginal.b);

figure
f = 1.5
set(gcf, 'Units', 'Inches', 'Position', [3, 3, 3.5*f, 3*f]);
dist = 300;
colors = {[0, 83, 143]./255, [34, 245, 179]./255};
plotPwithFakeContours(Poly1, colors{1}, 0.8, dist+3)
plotPwithFakeContours(Poly2, colors{2}, 0.4, dist)
xlabel("x", 'fontsize',16); ylabel("y", 'fontsize',16); zlabel("z", 'fontsize',16);
        
ax = gca ;
ax.FontSize = 14;
zlim(ax.ZLim*1.05)
xlim([ax.XLim(1)*1.05 ax.XLim(2)])
ylim([ax.YLim(1)*1.05 ax.YLim(2)])
ax.XTick = [];
ax.YTick = [];
ax.ZTick = [];
end