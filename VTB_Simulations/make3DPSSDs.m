clear;
close all;

% expNames = {    'Sim_100x100x30',...
%     'Sim_200x100x30',...
%     'Sim_300x100x30',...
%     'Sim_200x100x66',...
%     'Sim_300x100x66',...
%     'Sim_200x100x100',...
%     'Sim_300x100x100'};

expNames = {    'Sim_LithoAngles_100x100x30',...
    'Sim_LithoAngles_200x100x30',...
    'Sim_LithoAngles_300x100x30',...
    'Sim_LithoAngles_200x100x66',...
    'Sim_LithoAngles_300x100x66',...
    'Sim_LithoAngles_200x100x100',...
    'Sim_LithoAngles_300x100x100'};

N = length(expNames);

for kk=1:N
    load([expNames{kk},'.mat']);
    sizeOBB{kk} = outputStruct.results.DP.measuredParticleList;
    sizeML{kk} = outputStruct.results.DP.measuredParticleList_ML';
    
    
    PSSD_3D = generatePSSD(sizeML{kk} );
    dL = 4;
    maxL = [PSSD_3D.PSSD.grid(1).boundaries(end),...
        PSSD_3D.PSSD.grid(2).boundaries(end),...
        PSSD_3D.PSSD.grid(3).boundaries(end)];
    PSSD_3D = generatePSSD(sizeML{kk},ceil(maxL./dL));
    
    PSSD_3DObb = generatePSSD(sizeOBB{kk});
    dL = 4;
    maxL = [PSSD_3DObb.PSSD.grid(1).boundaries(end),...
        PSSD_3DObb.PSSD.grid(2).boundaries(end),...
        PSSD_3DObb.PSSD.grid(3).boundaries(end)];
    PSSD_3DObb = generatePSSD(sizeOBB{kk},ceil(maxL./dL));
    
    [x,y] = computeContourLevels(PSSD_3D);
    lvl = x;
    lvlMarginal = y;
    [x,y] = computeContourLevels(PSSD_3DObb);
    lvlObb = x;
    lvlMarginalObb = y;
     
    save([expNames{kk},'.mat'],'PSSD_3D','PSSD_3DObb',...
        'lvl','lvlObb','lvlMarginal','lvlMarginalObb','-append');
    
end