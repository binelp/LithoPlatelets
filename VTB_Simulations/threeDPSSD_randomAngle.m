close all
clear

expNames = {    'Sim_100x100x30',...
    'Sim_200x100x30',...
    'Sim_300x100x30',...
    'Sim_200x100x66',...
    'Sim_300x100x66',...
    'Sim_200x100x100',...
    'Sim_300x100x100'};

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
    sizeML{kk} = outputStruct.results.DP.measuredParticleList_ML'; %#ok<SAGROW>
end

%% MAKE 3D PSSD PLOTS
for kk=1:N
    if any([kk==1, kk==4, kk==6])
        figure;
        ax = axes(gcf); %#ok<LAXES>
        hold(ax,'on');
    end
    
    % Load data
    load([expNames{kk},'.mat']);

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
    zlim(ax,[0 150]);
    set(ax,'FontSize',8)
    set(ax,'XTick',[0   100   200   300   400],...
        'YTick',[0   100   200],...
        'ZTick',[0 50 100 150]);
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

%%
spreadFig;

return
%% Export all figures
fig = findobj('Type','figure');
for kk=1:length(fig)
    print(fig(kk),['plot',num2str(kk),'.pdf'],'-dpdf');
end