treeHome = '/imageData/binelp';


close all; 
pStruct = inspectStereoImage('/imageData/binelp/20210903_200x100x30_TH33/03092021-160516/data_camB/B-03092021-160516-7.csv',...
        'imgdir','/imageData/binelp/20210903_200x100x30/03092021-160516/data_camB');
    
close all; 
pStruct = inspectStereoImage('/imageData/binelp/20210903_300x100x30_TH33/03092021-163358/data_camA/A-03092021-163358-60.csv',...
        'imgdir','/imageData/binelp/20210903_300x100x30/03092021-163358/data_camB');
    

close all; 
pStruct = inspectStereoImage('/imageData/binelp/20210903_100x100x30_TH33/03092021-145039/data_camB/B-03092021-145039-34.csv',...
        'imgdir','/imageData/binelp/20210903_100x100x30/03092021-145039/data_camA');

