function fig = makeViolinPlot(x,y,color,targetsize)

fig = figure;
ax = axes(fig);
hold(ax,'on');

LineStyle = {':','--','-.'};
% targetsize = unique(targetsize);
for kk=1:length(targetsize)
    line([0,5],[targetsize(kk), targetsize(kk)],'Color','k','LineStyle',LineStyle{kk},'LineWidth',0.5)
end

violinplot({x,y},[1,2,3],'Width',0.37,'ShowData',false,'ViolinColor',{repmat(color,size(x,2),1),[0.85 0.85 0.85]},...
    'BoxWidth',0.13,'ViolinAlpha',0.85,'MarkerSize',5,'BoxColor',color,'LineWidth',0.25,...
    'ShowBox',false,'ShowMedian',false,'ShowWhiskers',false,'EdgeColor',[0.3 0.3 0.3]);


% chil = get(ax,'Children');
% chil = chil(1).Children;
% for kk=1:length(chil)
%     if strcmpi(chil(kk).Tag,'box')
%         chil(kk).LineWidth = 8;
%     end
% end
% for kk=1:length(chil)
%     if strcmpi(chil(kk).Tag,'median')
%         chil(kk).LineWidth = 0.75;
%     end
% end
% for kk=1:length(chil)
%     if strcmpi(chil(kk).Tag,'whisker')
%         chil(kk).LineWidth = 0.75;
%     end
% end


ylim(ax,[0 400])
xlim(ax,[0.5 3.5])
typeSetPlot(fig,'half','padding',0.2)
ylabel(ax,['Size [',char(181),'m]'])
set(ax,'TickLabelInterpreter', 'tex');
set(ax,'XTick',[1 2 3],'XTickLabel',{'{\it L}_1','{\it L}_2','{\it L}_3'})
set(ax,'YTick',[0 100 200 300 400])
set(ax,'LineWidth',0.75)
set(fig,'Position',[559   549   156   320],'Renderer','painters');

typeSetPlot(fig,'half','padding',0.15)

end