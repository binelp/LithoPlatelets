clear
close all

myPath = '/imageData/binelp';
suffix = '_TH33_PROCESSED';

expNames = {    '100x100x30',...
                '200x100x30',...
                '300x100x30',...
                '200x100x66',...
                '300x100x66',...
                '200x100x100',...
                '300x100x100'};
            
expDates = {    '20210903',...
                '20210903',...
                '20210903',...
                '20210812',...
                '20210812',...
                '20210812',...
                '20210812'};

                        
VHLowThreshold  =   [1E5,   0.4E6,  5.5E5,  0.6E6,  1E6,    1E6,    2E6     ];
VHHighThreshold =   [1E6,   1.5E6,  2.2E6,  3E6,    4E6,    4.4E6,  6.9E6   ];

L1LowThreshold  =   [70,    150,    255,    150,    250,    170,    250     ];
L1HighThreshold =   [140,   260,    350,    250,    350,    260,    360     ];
    
L2LowThreshold  =   [60,    60,     85,     85,      80,    90,     90      ];
L2HighThreshold =   [120,   140,    127,    150,    170,    180,    200     ];

L3LowThreshold  =   [7,     0,      0,      32,      24,    50,     50      ];
L3HighThreshold =   [Inf,   80,     70,     Inf,    100,    140,    130     ];


for kk=1:length(expNames)
    fullExpPath = [myPath,filesep,expDates{kk},'_',expNames{kk},suffix];
    
    sizedataCuboid = getProp(fullExpPath,'sizedataCuboid');
    sizedataCuboidObb = getProp(fullExpPath,'sizedataCuboidObb');
    sizedataNeedle = getProp(fullExpPath,'sizedataNeedleML');
    sizedataNeedleObb = getProp(fullExpPath,'sizedataNeedleObb');
    sizedataPlatelet = getProp(fullExpPath,'sizedataPlatelet');
    sizedataPlateletObb = getProp(fullExpPath,'sizedataPlateletObb');
      
    checkOutlier = @(data)any( [data(:,1)<L1LowThreshold(kk),data(:,1)>L1HighThreshold(kk),...
                                data(:,2)<L2LowThreshold(kk),data(:,2)>L2HighThreshold(kk),...
                                data(:,3)<L3LowThreshold(kk),data(:,3)>L3HighThreshold(kk),...
                                prod(data,2)<VHLowThreshold(kk),prod(data,2)>VHHighThreshold(kk)],2);
    
    
    outlierCuboid = checkOutlier(sizedataCuboid);
    outlierNeedle =  checkOutlier(sizedataNeedle);
    outlierPlatelet =  checkOutlier(sizedataPlatelet);
    
    nOutlierCuboid = sum(outlierCuboid);
    nOutlierNeedle = sum(outlierNeedle);
    nOutlierPlatelet = sum(outlierPlatelet);
    
    sizedataCuboid(outlierCuboid,:) = [];
    sizedataNeedle(outlierNeedle,:) = [];
    sizedataPlatelet(outlierPlatelet,:) = [];
    
    sizedataCuboidObb(outlierCuboid,:) = [];
    sizedataNeedleObb(outlierNeedle,:) = [];
    sizedataPlateletObb(outlierPlatelet,:) = [];
    
    catsize = [sizedataCuboid;sizedataNeedle;sizedataPlatelet];
    catsizeObb = [sizedataCuboidObb;sizedataNeedleObb;sizedataPlateletObb];
    allMean = mean(catsize);
    allStd = std(catsize);
    
    
    PSSD_3D = generatePSSD(catsize);
    dL = 4;
    maxL = [PSSD_3D.PSSD.grid(1).boundaries(end),...
        PSSD_3D.PSSD.grid(2).boundaries(end),...
        PSSD_3D.PSSD.grid(3).boundaries(end)];
    PSSD_3D = generatePSSD(catsize,ceil(maxL./dL));
    
%     lvl = computeContourLevels(PSSD_3D);
    
    PSSD_3DObb = generatePSSD(catsizeObb);
    dL = 4;
    maxL = [PSSD_3DObb.PSSD.grid(1).boundaries(end),...
        PSSD_3DObb.PSSD.grid(2).boundaries(end),...
        PSSD_3DObb.PSSD.grid(3).boundaries(end)];
    PSSD_3DObb = generatePSSD(catsizeObb,ceil(maxL./dL));
    
%     lvlObb = computeContourLevels(PSSD_3DObb);
    
    save(['Exp_',expNames{kk},'.mat'],'sizedataCuboid','sizedataCuboidObb',...
        'sizedataNeedle','sizedataNeedleObb',...
        'sizedataPlatelet','sizedataPlateletObb',...
        'nOutlierCuboid','nOutlierNeedle','nOutlierPlatelet',...
        'allMean','allStd','PSSD_3D','PSSD_3DObb','catsize','catsizeObb'); 
end 
