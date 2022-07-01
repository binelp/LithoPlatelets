function processPng(treeHome,expName,outDir,L)

%/home/binelp/Desktop/08_Code_PopulationControl/KIDDO/BURST/bin/staticAnalysis -dir /imageData/binelp/20210415_20_65_90/15042021-173149 -o /imageData/binelp/20210415_20_65_90_TH33 -ts 15042021-173149 -th 33 -n 400
% treeHome = '/imageData/binelp';
% expName = '20210811_100x100x33';
% 
% outDir = '/imageData/binelp/20210811_100x100x33_TH33';  %without last /

L1 = L(1);
L2 = L(2);
L3 = L(3);

if strcmp(outDir(end),filesep)
    outDir(end) = [];
end
if strcmp(treeHome(end),filesep)
    treeHome(end) = [];
end

expDir = dir([treeHome,filesep,expName]);
subExps = {expDir.name};
subExps(~[expDir.isdir]) = [];
subExps(1:2) = [];      % remove . and ..

% Based on 10.1098/rspa.2000.0667
mincntpx = 2*round(L2+L3);
if (L2^2+L3^2)>=L1
    maxntpx = 1.2*round(2*sqrt(2)*sqrt(L1^2+L2^2+L3^2));
else 
    maxntpx = 1.2*round(2*(sqrt(L3^2+L2^2)+L1));
end

for kk=1:length(subExps)
    targetCmd = [which('staticAnalysis.'),...
        ' -dir ',[treeHome,filesep,expName],...
        ' -o ',outDir,...
        ' -ts ',subExps{kk},...
        ' -cntpx ',num2str(mincntpx),...
        ' -maxcntpx ',num2str(maxntpx),...
        ' -th 33',...
        ' -n 300'];
    system(targetCmd);
end
end