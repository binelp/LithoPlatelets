close all
clear

expNames = {    'Exp_100x100x30',...
    'Exp_200x100x30',...
    'Exp_300x100x30',...
    'Exp_200x100x66',...
    'Exp_300x100x66',...
    'Exp_200x100x100',...
    'Exp_300x100x100'};

colors = {  [11 171 224],...
    [0 114 201],...
    [11 51 212],...
    [214 100 11],...
    [227 40 11],...
    [148 227 79],...
    [84 214 75]         };

realSizes = [101 101 31.3;
    201.5 101.9 31.7;
    300.4 102.0 32.4;
    205.0 105.0 65.5;
    305.9 106.1 66.0;
    204.8 105.2 104.8;
    305.6 109.8 105.5];

realVolumes = realSizes(:,1).*realSizes(:,2).*realSizes(:,3);

colors = cellfun(@(x)x./255,colors,'UniformOutput',false);

N = length(expNames);

for kk=1:N
    load([expNames{kk},'.mat']);
    sizeOBB{kk} = catsizeObb;
    volumeOBB{kk} = catsizeObb(:,1).*catsizeObb(:,2).*catsizeObb(:,3);
    sizeML{kk} = catsize;
    volumeML{kk} = catsize(:,1).*catsize(:,2).*catsize(:,3);
end


for kk=1:N
    fig(kk) = figure;
    ax(kk) = axes(fig(kk));
    histogram(volumeML{kk}./realVolumes(kk),'FaceColor',colors{kk},'FaceAlpha',0.7,'EdgeColor','none','Normalization','pdf');
    hold(ax(kk),'on');
    histogram(volumeML{kk}./realVolumes(kk),'FaceColor','none','EdgeColor',[0.15 0.15 0.15],'EdgeAlpha',0.9,'Normalization','pdf','DisplayStyle','stairs');
    histogram(volumeOBB{kk}./realVolumes(kk),'FaceColor',[0.35 0.35 0.35],'FaceAlpha',0.2,'EdgeColor','none','Normalization','pdf');
    histogram(volumeOBB{kk}./realVolumes(kk),'FaceColor','none','EdgeColor',[0.15 0.15 0.15],'EdgeAlpha',0.9,'Normalization','pdf','DisplayStyle','stairs');
    xlabel('Normalized volume');
    ylabel('Probability density');
    xlim([0 6]);
    set(ax(kk),'XTick',[0 1 2 3 4 5 6]);
    typeSetPlot(fig(kk),'single','padding',0.15);
    set(fig(kk),'Position',[734.4000  459.0000  312.0000  245.0000]);
    set(ax(kk),'PlotBoxAspectRatio',[1.0000    0.9617    0.9617]);
    ax(kk).YAxis.Exponent = 1;
end

spreadFig;

figure;
insubplot(3,3,1,fig(1));
insubplot(3,3,2,fig(2));
insubplot(3,3,3,fig(3));
insubplot(3,3,4,fig(4));
insubplot(3,3,5,fig(5));
insubplot(3,3,7,fig(6));
insubplot(3,3,8,fig(7));
set(gcf,'Position',[474   233   787   797]);


print(gcf,'fig_name.pdf','-dpdf');

