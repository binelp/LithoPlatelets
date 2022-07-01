% Generate and saves the contour levers used to plot

expNames = {    '100x100x30',...
                '200x100x30',...
                '300x100x30',...
                '200x100x66',...
                '300x100x66',...
                '200x100x100',...
                '300x100x100'};
  
            
for kk=1:length(expNames)
    load(['Exp_',expNames{kk},'.mat']);
    [x,y] = computeContourLevels(PSSD_3D);
    lvl = x;
    lvlMarginal = y;
    [x,y] = computeContourLevels(PSSD_3DObb);
    lvlObb = x;
    lvlMarginalObb = y;
	save(['Exp_',expNames{kk},'.mat'],'lvl','lvlObb','lvlMarginal','lvlMarginalObb','-append');
end