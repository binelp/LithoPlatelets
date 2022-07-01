% This function is useful to pull a property from a set of experiments
% For instance if expPath = imageData/binelp/expName_PROCESSED
% and property = sizedataNeedle then the function will return the
% concatenated vector of the sizedataNeedle properties of all .mat files in
% the _PROCESSED folder

function R = getProp(expPath,property)

if nargin==0
    disp('sizedataSphere');
    disp('sizedataNeedle');
    disp('sizedataCuboid');
    disp('sizedataPlatelet');
    disp('sizedataCuboidImposedNeedle');
    disp('sizedataAgglo');
    disp('numSphere');
    disp('numNeedle');
    disp('numCuboid');
    disp('numPlatelet');
    disp('numAgglo');
    disp('totalParticles');
    disp('VHVolumesSphere');
    disp('VHVolumesNeedle');
    disp('VHVolumesCuboid');
    disp('VHVolumesPlatelet');
    disp('VHVolumesAgglo');
    return;
end

singleexpPaths = dir([expPath,filesep,'*.mat']);

R = [];

for k = 1:size(singleexpPaths,1)
    temporaryData = load(strcat(expPath,filesep,singleexpPaths(k).name));
    if ~isfield(temporaryData.outputMatrix,property)
        continue;
    end

    R = vertcat(R,temporaryData.outputMatrix.(property));
end


end