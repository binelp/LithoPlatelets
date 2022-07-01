%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ETH Zurich, Switzerland
% Separation Process Laboratory
%
% Project:  CrystOCAM 2.0
% Year:     2018
% MATLAB:   R2017a, Windows 64bit
% Authors:  Ashwin Kumar Rajagopalan (AK)
%           Stefan Boetschi (SB)
%
% Purpose:
% Plots a color coded contour plot for the initial and final distribution
%
% Last modified:
% - 2022-03-15, PB: Cleaned, modified, adapted
% - 2019-02-13, SB: Minor plotting bugfix (empty average dimensions)
% - 2019-01-31, SB: Average L1 lines become a black disk marker
% - 2019-01-28, SB: Average L1 lines become thin solid black
% - 2019-01-21, SB: char(956) -> char(181)
% - 2019-01-14, SB: Add average L1 lines to the plots
% - 2018-10-03, SB: Introduce optional input argument defaultFigureStyle
% - 2018-07-26, SB: Return subplot handles
% - 2018-06-13, SB: Adapt contour levels
% - 2018-04-05, AK: Minor cosmetic changes
% - 2018-03-09, SB: L1 and L2 inputs can be cell arrays
% - 2018-03-05, AK: Initial creation
%
% Input arguments:
% - q_cell:             Normalized, number-weighted PSSD,
%                       size(qV) = [length(L2) length(L1]),
%                       units are [um^-2]
% - L1:                 L1 pivot vector [um] or cell array of vectors
% - L2:                 L2 pivot vector [um] or cell array of vectors
% - colorPSSDs:         Color vector for seeds and products
% - defaultFigureStyle: Boolean. Optional input argument. If true, no
%                       figure styles will be set in this script.
%                       Default: false
% - numAvgL_cell:       Optional input argument. Cell array with the same
%                       length as q_cell and qV_cell. Each cell contains a
%                       1x2 vector with the number-weighted average
%                       dimensions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ax, fig] =  plot3DPSSD(PSSD, colorPSSDs, alpha,...
    defaultFigureStyle, contourLvlsNum)

[L1, L2, L3] = PSSD.PSSD.grid.y;
q_cell = {PSSD.PSSD.F};

% Set defaultFigureStyle
if nargin < 4
    defaultFigureStyle = false;
end

% Make sure that L1 and L2 are cell arrays
if ~iscell(L1) || ~iscell(L2) || ~iscell(L3)
    L1Vec = L1;
    L2Vec = L2;
    L3Vec = L3;
    L1 = cell(1,length(q_cell));
    L2 = cell(1,length(q_cell));
    L3 = cell(1,length(q_cell));
    for ii=1:length(q_cell)
        L1{ii} = [L1Vec(1)-1E5*eps L1Vec L1Vec(end)+1E5*eps];
        L2{ii} = [L2Vec(1)-1E5*eps L2Vec L2Vec(end)+1E5*eps];
        L3{ii} = [L3Vec(1)-1E5*eps L3Vec L3Vec(end)+1E5*eps];
    end
end

% Contour levels and axis limits
% contourLvlsVol = [0.2,0.5,0.8];
if ~exist('contourLvlsNum','var')
    contourLvlsNum = [0.05,0.2,0.6];
end
numAxisLimits = [0 400 0 150 0 150];
% volAxisLimits = [0 200 0 200 0 200];

% Plot number weighted distributions
%     numWeightHandle = subplot(1,1,1);

ax = gca;
hold(ax,'on');

for ii = 1:length(q_cell)
    %%
    %% Here we permute dimension 1 and 2 to correct L1 and L2
    
    temp_Cell   = cell(1,1);
    temp_Cell{1}= permute(q_cell{1},[2 1 3]);
    q_cell      = temp_Cell;

    % Number-weighted PSSD
    for j = 1:length(contourLvlsNum)
        [x, y, z] = meshgrid(L1{ii},L2{ii},L3{ii});

        p = patch(ax, isosurface(x,y,z,padarray(q_cell{ii}./max(q_cell{ii}(:)),[1 1 1]), contourLvlsNum(j)));

%         normqcell = q_cell{ii}./max(q_cell{ii}(:));
%         logicalq = normqcell>contourLvlsNum(j);
%         subPSSD.F = q_cell{1};
%         subPSSD.F(~logicalq) = 0;
%         subPSSD.grid = PSSD.PSSD.grid;
%         subPSSD.F = permute(subPSSD.F,[2 1 3]);
%         submom(j) = compute3DCrossMoments(subPSSD); %#ok<AGROW>
        
        set(p,'facecolor',colorPSSDs(ii,:)+(j>1)*(1-colorPSSDs(ii,:))*(1/length(contourLvlsNum))^(j-1),'facealpha',alpha,'Edgecolor','none');
    end
    view(90+50,30)
end

% fullmom = compute3DCrossMoments(PSSD.PSSD);
% [submom.mu000]./fullmom.mu000

xlabel(ax,['{\it L}_1 [',char(181),'m]']);
ylabel(ax,['{\it L}_2 [',char(181),'m]']);
zlabel(ax,['{\it L}_3 [',char(181),'m]']);
axis(ax,numAxisLimits);
box(ax,'on');
grid(ax,'on');
daspect(ax,[1 1 1]);

if ~defaultFigureStyle
    ax.LineWidth = 1;
    ax.FontName = 'Helvetica';
    ax.FontSize = 8;
    ax.FontWeight = 'bold';
end

end