%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ETH Zurich, Switzerland
% Separation Processes Laboratory
%
% Project:  Lithoplatelets
% Year:     2021
% MATLAB:   R2019b, Windows 64bit
% Authors:  Anna Jaeggi (AJ)
%
% Purpose:
% Plot 2D distributions of physical angles as contourplots from both
% experimental and simulational angle distributions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all;

% ID of saved results to be used
runID = "";
% Directory where data and results stored
mainDir ="";
% Subdirectory for each population
dirs = ["20210812_200x100x66_TH33","20210812_200x100x100_TH33","20210812_300x100x66_TH33",...
    "20210812_300x100x100_TH33","20210903_100x100x30_TH33","20210903_200x100x30_TH33","20210903_300x100x30_TH33"];
% Lithoplatelet lengths
Ls = {[205,105,66],[205,105,105],[306,106,66],[306,110,106],[101,101,31],[202,102,32],[300,102,32]};
% Load cuboid geometry information
load("CrystalData.mat");
HM = {Crystals(15).H,Crystals(15).M};
physAnglesCell = {};

% Colors for lithoplatelets
colors = { [214 100 11]./255,...
[148 227 79]./255,...
[227 40 11]./255,...
[84 214 75]./255,...
[11 171 224]./255,...
[0 114 201]./255,...
[11 51 212]./255};

% Generate physical angles from uniformly distributed quaternions
q = uniformSampledQuats(100000);
randPhysAngles = getPhysicalAngles(q, [120,100,30], "uniformlyRandomPhysAngles.mat");

% Load physical angles from data and results directory
for kk = 1:7
    L = Ls{kk};
    savePath = join([mainDir,dirs(kk),"\physicalAngles_",runID,".mat"],"");
    load(savePath)
    physAnglesCell{kk} = physicalAngles;
end
%% Contourplot experimetnal main part 
% Binning
nbins = 20;
x = randPhysAngles.flowAngleVec;
y = randPhysAngles.flowAnglePlane;
% Obtain counts
Xedges=linspace(0,90,nbins+1);
Yedges=linspace(0,90,nbins+1);
Nrand = histcounts2(x,y,Xedges,Yedges);
Nrand = Nrand/nansum(Nrand,'all');

% Index for corresponding lihtoplatelet lengths
parts = [6,5,4, 1, 2, 3, 7];

levels = 0:0.5:2.6;
for ii = 1:7
    % Binning
    x = physAnglesCell{parts(ii)}.flowAngleVec;
    y = physAnglesCell{parts(ii)}.flowAnglePlane;
    % Create bins and obtain counts
    Xedges=linspace(0,90,nbins+1);
    Yedges=linspace(0,90,nbins+1);
    N = histcounts2(x,y,Xedges,Yedges);

    [X,Y] = meshgrid(linspace(0,90,nbins),...
        linspace(0,90,nbins));
    
    % Tranform to make uniform uniform
    N = N/nansum(N,'all');

    N = N./Nrand;
    figs{ii} = figure;
    hold on;
    title(sprintf("L = [%i %i %i]",Ls{parts(ii)}(1),Ls{parts(ii)}(2),Ls{parts(ii)}(3)))
    c = pink(100);
    colormap(c(15:80,:));
    contourf(X,Y,N',levels,"ShowText","off")
    caxis([min(levels),max(levels)])
    tx = [2, 88,2,2];
    ty = [88,2,2,88];
    fill(tx,ty,'w', 'EdgeAlpha', 0 )
    xlim([0,90])
    ylim([0,90])
    xlabel("\theta")
    ylabel("\phi")
    ax = gca;
    if ii == 1
        typeSetPlot(figs{ii}, 'single', 'noPADDING')
        f = 0.85;
        set(ax,'PlotBoxAspectRatio',[ 1.0679    1.0000    1.0000]);
        set(figs{ii}, 'Position', [680 558 312*f 420.0000*f]);
        
    else
        typeSetPlot(figs{ii}, 'half', 'noPADDING')
        set(ax,'PlotBoxAspectRatio',[ 1.0000    0.9364    0.9364]);
        set(figs{ii},'Position',[377 441 156 157]);
    end
    print(gcf,sprintf('angleDistContour%i.pdf',ii),'-dpdf');
    
    h = colorbar;
    h.Limits = [min(levels),max(levels)];
    
    set(figs{ii}, 'Position', [680 558 312*f 420.0000*f]);
    print(gcf,sprintf('angleDistContourColorbar%i.pdf',ii),'-dpdf');
end

%% Contourplot SI Experimental
% Binning
nbins = 20;
x = randPhysAngles.flowAngleVec;
y = randPhysAngles.flowAnglePlane;
% Obtain counts
Xedges=linspace(0,90,nbins+1);
Yedges=linspace(0,90,nbins+1);
Nrand = histcounts2(x,y,Xedges,Yedges);
    Nrand = Nrand/nansum(Nrand,'all');

% Index for corresponding lihtoplatelet lengths
parts = [6,5,4, 1, 2, 3, 7];

levels = 0:0.5:2.6;
for ii = 1:7
    % Binning
    x = physAnglesCell{parts(ii)}.flowAngleVec;
    y = physAnglesCell{parts(ii)}.flowAnglePlane;
    % Create bins and obtain counts
    Xedges=linspace(0,90,nbins+1);
    Yedges=linspace(0,90,nbins+1);
    N = histcounts2(x,y,Xedges,Yedges);
    [X,Y] = meshgrid(linspace(0,90,nbins),...
        linspace(0,90,nbins));
    
    % Tranform to make uniform uniform
    N = N/nansum(N,'all');

    N = N./Nrand;
    max(N,[],'all')
    sum(N>10,'all')/length(N)

    % Create countor, transponsing N so rows correspond to Y and columns correspond to X
    figs{ii} = figure;
    hold on;

    title(sprintf("L = [%i %i %i]",Ls{parts(ii)}(1),Ls{parts(ii)}(2),Ls{parts(ii)}(3)))
    c = pink(100);
    colormap(c(15:80,:));
    contourf(X,Y,N',levels,"ShowText","off")
    caxis([min(levels),max(levels)])
    tx = [2, 88,2,2];
    ty = [88,2,2,88];
    fill(tx,ty,'w', 'EdgeAlpha', 0 )
    xlim([0,90])
    ylim([0,90])
    xlabel("\theta")
    ylabel("\phi")
    ax = gca;
    if ii == 1
        typeSetPlot(figs{ii}, 'single', 'noPADDING')
        f = 0.85;
        set(ax,'PlotBoxAspectRatio',[ 1.0679    1.0000    1.0000]);
        set(figs{ii}, 'Position', [680 558 312*f 420.0000*f]);
        
    else
        typeSetPlot(figs{ii}, 'half', 'noPADDING')
        set(ax,'PlotBoxAspectRatio',[ 1.0000    0.9364    0.9364]);
        set(figs{ii},'Position',[377 441 156 157]);
    end
    newColormap()
    print(gcf,sprintf('angleDistContour%i.pdf',ii),'-dpdf');
    
    h = colorbar;
    h.Limits = [min(levels),max(levels)];
    newColormap()
    set(figs{ii}, 'Position', [680 558 312*f 420.0000*f]);
    print(gcf,sprintf('angleDistContourColorbar%i.pdf',ii),'-dpdf');
end

%% Contourplot from Simulations Supplementary Information
% Paths to files obtained from the angle estimation simulation
files = {};
% Index for corresponding lihtoplatelet lengths
parts = [4,5,7];
for ii = 1:length(files)
    file = files{ii};
    % Load physical angles, or calculate them if not already done
    try
        load(join(["Results/",file(1:end-4),"_measPhysical.mat"],""));
        measAngles = physicalAngles;
    catch
        load(join(["D:\SimResultsCuboidAngles\",file],""))

        measQuats = outputStruct.results.measQuats;
        realQuats = outputStruct.realQuats;
        L = outputStruct.L;
        fOpts = outputStruct.results.fOpts;

        % Remove worst fits
        removalThreshhold = 1.5;
        okIx = fOpts < removalThreshhold;
        fOpts = fOpts(okIx);
        measQuats = measQuats(okIx,:);
        realQuats = realQuats(okIx,:);
        N = length(fOpts);

        getPhysicalAngles(measQuats, L, join(["Results/",file(1:end-4),"_measPhysical.mat"],""))
    end
    % Contourplot for each case
    for phA = [measAngles]
        % Binning
        x = phA.flowAngleVec;
        y = phA.flowAnglePlane;
        % Create bins and obtain counts
        Xedges = linspace(0,90,nbins+1);
        Yedges = linspace(0,90,nbins+1);
        N = histcounts2(x,y,Xedges,Yedges);

        [X,Y] = meshgrid(linspace(0,90,nbins),...
            linspace(0,90,nbins));

        % Tranform to make uniform uniform
        N = N/nansum(N,'all');
        N = N./Nrand;

        figure;
        title(sprintf("L = [%i %i %i]",Ls{parts(ii)}(1),Ls{parts(ii)}(2),Ls{parts(ii)}(3)))
        hold on;
        c = pink(100);
        colormap(c(15:80,:));
        contourf(X,Y,N',levels,"ShowText","off")
        caxis([min(levels),max(levels)])
        tx = [2, 88,2,2];
        ty = [88,2,2,88];
        fill(tx,ty,'w', 'EdgeAlpha', 0 )
        xlim([0,90])
        ylim([0,90])
        xlabel("\theta")
        ylabel("\phi")

        ax = gca;

        typeSetPlot(gcf, 'half', 'noPADDING')
        set(ax,'PlotBoxAspectRatio',[ 1.0000    0.9364    0.9364]);
        set(gcf,'Position',[377 441 156 157]);
        newColormap()
        print(gcf,sprintf('anglevtb_10_20_%i.pdf',ii),'-dpdf');
        colorbar
        print(gcf,sprintf('anglevtb_10_20_%i_COLORBAR.pdf',ii),'-dpdf');
    
    end
end

%% Color palette used in plots
function []= newColormap()
colormap([0.7490    0.3373    0.4745
    0.7496    0.3397    0.4750
    0.7502    0.3422    0.4756
    0.7509    0.3446    0.4761
    0.7515    0.3471    0.4767
    0.7521    0.3496    0.4772
    0.7527    0.3520    0.4777
    0.7533    0.3545    0.4783
    0.7539    0.3569    0.4788
    0.7546    0.3594    0.4794
    0.7552    0.3619    0.4799
    0.7558    0.3643    0.4804
    0.7564    0.3668    0.4810
    0.7570    0.3692    0.4815
    0.7576    0.3717    0.4820
    0.7582    0.3742    0.4826
    0.7589    0.3766    0.4831
    0.7595    0.3791    0.4837
    0.7601    0.3815    0.4842
    0.7607    0.3840    0.4847
    0.7613    0.3865    0.4853
    0.7619    0.3889    0.4858
    0.7626    0.3914    0.4864
    0.7632    0.3938    0.4869
    0.7638    0.3963    0.4874
    0.7644    0.3988    0.4880
    0.7650    0.4012    0.4885
    0.7656    0.4037    0.4890
    0.7662    0.4062    0.4896
    0.7669    0.4086    0.4901
    0.7675    0.4111    0.4907
    0.7681    0.4135    0.4912
    0.7687    0.4160    0.4917
    0.7693    0.4185    0.4923
    0.7699    0.4209    0.4928
    0.7705    0.4234    0.4933
    0.7712    0.4258    0.4939
    0.7718    0.4283    0.4944
    0.7724    0.4308    0.4950
    0.7730    0.4332    0.4955
    0.7736    0.4357    0.4960
    0.7742    0.4381    0.4966
    0.7749    0.4406    0.4971
    0.7755    0.4431    0.4977
    0.7761    0.4455    0.4982
    0.7767    0.4480    0.4987
    0.7773    0.4504    0.4993
    0.7779    0.4529    0.4998
    0.7785    0.4554    0.5003
    0.7792    0.4578    0.5009
    0.7798    0.4603    0.5014
    0.7804    0.4627    0.5020
    0.7822    0.4656    0.5033
    0.7839    0.4684    0.5047
    0.7857    0.4713    0.5061
    0.7875    0.4741    0.5075
    0.7892    0.4770    0.5089
    0.7910    0.4798    0.5103
    0.7928    0.4827    0.5116
    0.7945    0.4855    0.5130
    0.7963    0.4884    0.5144
    0.7981    0.4912    0.5158
    0.7998    0.4940    0.5172
    0.8016    0.4969    0.5186
    0.8034    0.4997    0.5200
    0.8052    0.5026    0.5213
    0.8069    0.5054    0.5227
    0.8087    0.5083    0.5241
    0.8105    0.5111    0.5255
    0.8122    0.5140    0.5269
    0.8140    0.5168    0.5283
    0.8158    0.5196    0.5296
    0.8175    0.5225    0.5310
    0.8193    0.5253    0.5324
    0.8211    0.5282    0.5338
    0.8228    0.5310    0.5352
    0.8246    0.5339    0.5366
    0.8264    0.5367    0.5379
    0.8281    0.5396    0.5393
    0.8299    0.5424    0.5407
    0.8317    0.5453    0.5421
    0.8334    0.5481    0.5435
    0.8352    0.5509    0.5449
    0.8370    0.5538    0.5463
    0.8388    0.5566    0.5476
    0.8405    0.5595    0.5490
    0.8423    0.5623    0.5504
    0.8441    0.5652    0.5518
    0.8458    0.5680    0.5532
    0.8476    0.5709    0.5546
    0.8494    0.5737    0.5559
    0.8511    0.5765    0.5573
    0.8529    0.5794    0.5587
    0.8547    0.5822    0.5601
    0.8564    0.5851    0.5615
    0.8582    0.5879    0.5629
    0.8600    0.5908    0.5642
    0.8617    0.5936    0.5656
    0.8635    0.5965    0.5670
    0.8653    0.5993    0.5684
    0.8671    0.6022    0.5698
    0.8688    0.6050    0.5712
    0.8706    0.6078    0.5725
    0.8712    0.6112    0.5741
    0.8718    0.6146    0.5756
    0.8724    0.6180    0.5772
    0.8730    0.6214    0.5787
    0.8737    0.6248    0.5802
    0.8743    0.6281    0.5818
    0.8749    0.6315    0.5833
    0.8755    0.6349    0.5849
    0.8761    0.6383    0.5864
    0.8767    0.6417    0.5879
    0.8774    0.6451    0.5895
    0.8780    0.6484    0.5910
    0.8786    0.6518    0.5925
    0.8792    0.6552    0.5941
    0.8798    0.6586    0.5956
    0.8804    0.6620    0.5972
    0.8810    0.6654    0.5987
    0.8817    0.6687    0.6002
    0.8823    0.6721    0.6018
    0.8829    0.6755    0.6033
    0.8835    0.6789    0.6048
    0.8841    0.6823    0.6064
    0.8847    0.6857    0.6079
    0.8854    0.6890    0.6095
    0.8860    0.6924    0.6110
    0.8866    0.6958    0.6125
    0.8872    0.6992    0.6141
    0.8878    0.7026    0.6156
    0.8884    0.7060    0.6171
    0.8890    0.7093    0.6187
    0.8897    0.7127    0.6202
    0.8903    0.7161    0.6218
    0.8909    0.7195    0.6233
    0.8915    0.7229    0.6248
    0.8921    0.7263    0.6264
    0.8927    0.7296    0.6279
    0.8933    0.7330    0.6295
    0.8940    0.7364    0.6310
    0.8946    0.7398    0.6325
    0.8952    0.7432    0.6341
    0.8958    0.7466    0.6356
    0.8964    0.7499    0.6371
    0.8970    0.7533    0.6387
    0.8977    0.7567    0.6402
    0.8983    0.7601    0.6418
    0.8989    0.7635    0.6433
    0.8995    0.7669    0.6448
    0.9001    0.7702    0.6464
    0.9007    0.7736    0.6479
    0.9013    0.7770    0.6494
    0.9020    0.7804    0.6510
    0.9029    0.7825    0.6528
    0.9039    0.7847    0.6547
    0.9048    0.7869    0.6565
    0.9058    0.7890    0.6584
    0.9068    0.7912    0.6602
    0.9077    0.7933    0.6621
    0.9087    0.7955    0.6639
    0.9097    0.7976    0.6657
    0.9106    0.7998    0.6676
    0.9116    0.8019    0.6694
    0.9125    0.8041    0.6713
    0.9135    0.8062    0.6731
    0.9145    0.8084    0.6750
    0.9154    0.8105    0.6768
    0.9164    0.8127    0.6787
    0.9173    0.8148    0.6805
    0.9183    0.8170    0.6824
    0.9193    0.8191    0.6842
    0.9202    0.8213    0.6860
    0.9212    0.8235    0.6879
    0.9221    0.8256    0.6897
    0.9231    0.8278    0.6916
    0.9241    0.8299    0.6934
    0.9250    0.8321    0.6953
    0.9260    0.8342    0.6971
    0.9270    0.8364    0.6990
    0.9279    0.8385    0.7008
    0.9289    0.8407    0.7027
    0.9298    0.8428    0.7045
    0.9308    0.8450    0.7063
    0.9318    0.8471    0.7082
    0.9327    0.8493    0.7100
    0.9337    0.8514    0.7119
    0.9346    0.8536    0.7137
    0.9356    0.8557    0.7156
    0.9366    0.8579    0.7174
    0.9375    0.8601    0.7193
    0.9385    0.8622    0.7211
    0.9394    0.8644    0.7230
    0.9404    0.8665    0.7248
    0.9414    0.8687    0.7266
    0.9423    0.8708    0.7285
    0.9433    0.8730    0.7303
    0.9443    0.8751    0.7322
    0.9452    0.8773    0.7340
    0.9462    0.8794    0.7359
    0.9471    0.8816    0.7377
    0.9481    0.8837    0.7396
    0.9491    0.8859    0.7414
    0.9500    0.8880    0.7433
    0.9490    0.9216    0.7176
    0.9500    0.9231    0.7200
    0.9510    0.9246    0.7224
    0.9520    0.9262    0.7248
    0.9530    0.9277    0.7272
    0.9540    0.9293    0.7296
    0.9550    0.9308    0.7319
    0.9560    0.9323    0.7343
    0.9570    0.9339    0.7367
    0.9580    0.9354    0.7391
    0.9590    0.9369    0.7415
    0.9600    0.9385    0.7439
    0.9610    0.9400    0.7463
    0.9620    0.9416    0.7486
    0.9630    0.9431    0.7510
    0.9640    0.9446    0.7534
    0.9650    0.9462    0.7558
    0.9660    0.9477    0.7582
    0.9670    0.9493    0.7606
    0.9680    0.9508    0.7629
    0.9690    0.9523    0.7653
    0.9700    0.9539    0.7677
    0.9710    0.9554    0.7701
    0.9720    0.9569    0.7725
    0.9730    0.9585    0.7749
    0.9740    0.9600    0.7772
    0.9750    0.9616    0.7796
    0.9760    0.9631    0.7820
    0.9770    0.9646    0.7844
    0.9780    0.9662    0.7868
    0.9790    0.9677    0.7892
    0.9800    0.9692    0.7915
    0.9810    0.9708    0.7939
    0.9820    0.9723    0.7963
    0.9830    0.9739    0.7987
    0.9840    0.9754    0.8011
    0.9850    0.9769    0.8035
    0.9860    0.9785    0.8058
    0.9870    0.9800    0.8082
    0.9880    0.9815    0.8106
    0.9890    0.9831    0.8130
    0.9900    0.9846    0.8154
    0.9910    0.9862    0.8178
    0.9920    0.9877    0.8201
    0.9930    0.9892    0.8225
    0.9940    0.9908    0.8249
    0.9950    0.9923    0.8273
    0.9960    0.9938    0.8297
    0.9970    0.9954    0.8321
    0.9980    0.9969    0.8344
    0.9990    0.9985    0.8368
    1.0000    1.0000    0.8392])
end