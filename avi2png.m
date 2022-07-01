%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ETH Zurich, Switzerland
% Separation Processes Laboratory
%
% Project:  Roche 2
% Year:     2021
% MATLAB:   R2019b, Linux 64bit
% Authors:  Pietro Binel (PB)
%
% Purpose:          Extract an .avi uncompressed file into grayscale .png
%                   (lossless) images for improved storage
% Input arguments:
% - treeHome:       path of the experiment e.g. /imageData/binelp
% - expName:        name of the experiment e.g. 20210811_100x100x33
%
% Last modified:
% - 2021-10-05, PB: Initial creation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ######################    ONLY WORKS FOR 300 IMAGES AS SUCH

function avi2png(treeHome,expName)

if strcmp(treeHome(end),filesep)
    treeHome(end) = [];
end

expDir = dir([treeHome,filesep,expName]);
subExps = {expDir.name};
subExps(~[expDir.isdir]) = [];
subExps(1:2) = [];      % remove . and ..

%length(subExps)
for kk=1:length(subExps)
    baseVideoNameA = ['A-',subExps{kk},'-000'];
    baseNameA = ['A-',subExps{kk},'-'];
    extractAvi({[baseVideoNameA,'0'],[baseVideoNameA,'1'],[baseVideoNameA,'2']},...
        'inputPath',[treeHome,filesep,expName,filesep,subExps{kk}],...
        'outputPath',[treeHome,filesep,expName,filesep,subExps{kk},filesep,'data_camA'],...
        'basename',{baseNameA,baseNameA,baseNameA},...
        'q',100,'startn',[0, 143, 286]);
    
    baseVideoNameB = ['B-',subExps{kk},'-000'];
    baseNameB = ['B-',subExps{kk},'-'];
    extractAvi({[baseVideoNameB,'0'],[baseVideoNameB,'1'],[baseVideoNameB,'2']},...
        'inputPath',[treeHome,filesep,expName,filesep,subExps{kk}],...
        'outputPath',[treeHome,filesep,expName,filesep,subExps{kk},filesep,'data_camB'],...
        'basename',{baseNameB,baseNameB,baseNameB},...
        'q',100,'startn',[0, 143, 286]);
    
end
end