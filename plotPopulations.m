close all
clear

expNames = {    'Exp_100x100x30',...
    'Exp_200x100x30',...
    'Exp_300x100x30',...
    'Exp_200x100x66',...
    'Exp_300x100x66',...
    'Exp_200x100x100',...
    'Exp_300x100x100'};

realSizes = [101 101 31.3;
    201.5 101.9 31.7;
    300.4 102.0 32.4;
    205.0 105.0 65.5;
    305.9 106.1 66.0;
    204.8 105.2 104.8;
    305.6 109.8 105.5];


colors = {  [11 171 224],...
    [0 114 201],...
    [11 51 212],...
    [214 100 11],...
    [227 40 11],...
    [148 227 79],...
    [84 214 75]         };

colors = cellfun(@(x)x./255,colors,'UniformOutput',false);

N = length(expNames);

for kk=1:N
    load([expNames{kk},'.mat']);
    sizeOBB{kk} = catsizeObb;
    sizeML{kk} = catsize;
end

averagesML = cell2mat(cellfun(@(x)mean(x),sizeML(:),'UniformOutput',false));
stdsML = cell2mat(cellfun(@(x)std(x),sizeML(:),'UniformOutput',false));
averagesOBB = cell2mat(cellfun(@(x)mean(x),sizeOBB(:),'UniformOutput',false));
stdsOBB = cell2mat(cellfun(@(x)std(x),sizeOBB(:),'UniformOutput',false));

%% MAKE 3D PSSD PLOTS
for kk=1:N
    if any([kk==1, kk==4, kk==6])
        figure;
        ax = axes(gcf); %#ok<LAXES>
        hold(ax,'on');
    end
    
    % Load data
    load([expNames{kk},'.mat']);
    
        PSSD_3D = PSSD_3DObb;
        lvl = lvlObb;
        lvlMarginal = lvlMarginalObb;

    plot3DPSSD(PSSD_3D,colors{kk},0.3,false,lvl);
    alpha = 0.3;

    for ignoredDim=1:3
        dims = find([1,2,3]~=ignoredDim);
        q_m  = {squeeze(sum(PSSD_3D.PSSD.F,ignoredDim))};
        plot3DPSSD_marginal(dims, q_m,...
            PSSD_3D.PSSD.grid(dims(1)).y, PSSD_3D.PSSD.grid(dims(2)).y,...
            colors{kk}, alpha, 2, [], [], lvlMarginal{ignoredDim});
        hold on;
    end
    
    xlim(ax,[0 400]);
    ylim(ax,[0 200]);
    zlim(ax,[0 200]);
    set(ax,'FontSize',8)
    set(ax,'XTick',[0   100   200   300   400],...
        'YTick',[0   100   200],...
        'ZTick',[0 50 100 150 200]);
end

% Find all figures and axes
fig = findobj('Type','figure');
ax = get(fig,'children');
ax = [ax{:}];

% Change size and appearance 
for kk=1:length(fig)
    view(ax(kk),151.4324,24.0921);
    set(fig(kk),'Renderer','painters')
    typeSetPlot(fig(kk),'single','padding',0);
    typeSetPlot(fig(kk),'single','padding',0);
end

%% MAKE PARITY PLOTS
% 
% for kk=1:3
%     figB(kk) = figure;
%     axB(kk) = axes(figB(kk));
%     line(gca,[0 500],[0 500],'Color','k','LineWidth',0.75,'LineStyle','-');
%     hold on
%     scatter(realSizes(:,kk),averagesML(:,kk),37,vertcat(colors{:}),'filled','MarkerFaceAlpha',0.7)
%     scatter(realSizes(:,kk),averagesOBB(:,kk),30,vertcat(colors{:}))
%     box on
%     xlabel(['True{\it L}_',num2str(kk),' [',char(181),'m]'],'FontSize',8)
%     ylabel(['Measured{\it L}_',num2str(kk),' [',char(181),'m]'],'FontSize',8)
%     axis equal   
%     grid on 
%     set(gca,'FontSize',8)
% end
% 
% axis(axB(1),[0 400 0 400])
% set(axB(1),'XTick',[0 100 200 300 400]);
% set(figB(1),'Position',[2069         613         237         162]);
% 
% axis(axB(2),[0 200 0 200])
% set(axB(2),'XTick',[0 50 100 150 200]);
% set(figB(2),'Position',[2069         613         237         162]);
% 
% axis(axB(3),[0 150 0 150])
% set(axB(3),'XTick',[0 50 100 150]);
% set(figB(3),'Position',[2069         613         237         162]);


%% MAKE VIOLIN PLOTS

for kk=1:N
    % pull from https://github.com/bastibe/Violinplot-Matlab
    violinFig(kk) = makeViolinPlot(sizeML{kk},sizeOBB{kk},colors{kk},realSizes(kk,:));
end

% return
% close all
% makeboxplot2(sizeML{1},sizeOBB{1},colors{1},[100,33]);
% % export_fig(gcf,'1.png','-painters','-transparent','-r600')
% makeboxplot2(sizeML{2},sizeOBB{2},colors{2},[200,100,33]);
% % export_fig(gcf,'2.png','-painters','-transparent','-r600')
% makeboxplot2(sizeML{3},sizeOBB{3},colors{3},[300,100,33]);
% % export_fig(gcf,'3.png','-painters','-transparent','-r600')
% makeboxplot2(sizeML{4},sizeOBB{4},colors{4},[200,100,66]);
% % export_fig(gcf,'4.png','-painters','-transparent','-r600')
% makeboxplot2(sizeML{5},sizeOBB{5},colors{5},[300,100,66]);
% % export_fig(gcf,'5.png','-painters','-transparent','-r600')
% makeboxplot2(sizeML{6},sizeOBB{6},colors{6},[200,100]);
% % export_fig(gcf,'6.png','-painters','-transparent','-r600')
% makeboxplot2(sizeML{7},sizeOBB{7},colors{7},[300,100]);
% % export_fig(gcf,'7.png','-painters','-transparent','-r600')
% spreadFig;

%% Create ;-separated strings for excel to make Latex table
tbReal = [];
tbML = [];
tbOBB = [];
separator = ';';
realSizes = round(realSizes);

reshapedRealSizes = reshape(realSizes',1,3*N);
reshapedAveragesML = reshape(averagesML',1,3*N);
reshapedAveragesOBB = reshape(averagesOBB',1,3*N);
reshapedStdsML = reshape(stdsML',1,3*N);
reshapedStdsOBB = reshape(stdsOBB',1,3*N);

for kk=1:length(reshapedAveragesML)
    tbReal = strcat(tbReal,...
         num2str(reshapedRealSizes(kk)), ';');
    tbML = strcat(tbML,...
        sprintf('%.1f',round(reshapedAveragesML(kk),1)), '(', sprintf('%.0f',10*round(reshapedStdsML(kk),1)), ');');  
    tbOBB = strcat(tbOBB,...
        sprintf('%.1f',round(reshapedAveragesOBB(kk),1)), '(', sprintf('%.0f',10*round(reshapedStdsOBB(kk),1)), ');');
end

%%
spreadFig;


klist=1:5;%the number of clusters you want to try
myfunc = @(X,K)(kmeans(X, K));
eva = evalclusters(sizeOBB{5},myfunc,'silhouette','klist',klist)
classes=kmeans(sizeOBB{5},eva.OptimalK);

return
%% Export all figures
fig = findobj('Type','figure');
for kk=1:length(fig)
    print(fig(kk),['plot',num2str(kk),'.pdf'],'-dpdf');
end