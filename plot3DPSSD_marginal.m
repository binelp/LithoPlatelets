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
%           Anna Jaeggi (AJ)
%
% Purpose:
% Plots a color coded contour plot for the initial and final distribution
%
% Last modified:
% - 2022-03-15, PB: Cleaned, modified, adapted
% - 2019-10-11, AJ: Adapted for 3D marginal distribution
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
% - plane:              1x2 vector indicating which plane to plot the
%                       marginal distribution in ([1,2]->xy plane,
%                       [1,3] -> xz plane ect...)
% - q_cell:             Normalized, number-weighted PSSD,
%                       size(qV) = [length(L2) length(L1]),
%                       units are [um^-2]
% - qV_cell:            Normalized, volume-weighted PSSD,
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
% - volAvgL_cell:       Optional input argument. Cell array with the same
%                       length as q_cell and qV_cell. Each cell contains a
%                       1x2 vector with the volume-weighted average
%                       dimensions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ax, fig] = plot3DPSSD_marginal(plane, q_cell, L1, L2, colorPSSDs,...
    opacity, offset, defaultFigureStyle, numAvgL_cell, contourLvlsNum)

% Set defaultFigureStyle
if nargin < 8
    defaultFigureStyle = false;
end
% Check optional input arguments for average sizes
if nargin < 9
    numAvgL_cell = [];
end

% Make sure that L1 and L2 are cell arrays
if ~iscell(L1) || ~iscell(L2)
    L1Vec = L1;
    L2Vec = L2;
    L1 = cell(1,length(q_cell));
    L2 = cell(1,length(q_cell));
    for ii=1:length(q_cell)
        L1{ii} = L1Vec;
        L2{ii} = L2Vec;
    end
end

% Contour levels and axis limits
%     contourLvlsVol = [0.1,0.5,0.9];
if ~exist('contourLvlsNum','var')
    contourLvlsNum = [0.05,0.2,0.6];
end
numAxisLimits = [0 300 0 300 0 150];
%     limval = max([cell2mat(L1),cell2mat(L2)]);
%     volAxisLimits = [0 limval 0 limval 0 limval];

ax = gca;
hold(ax,'on');

for ii = 1:length(q_cell)
    % Number-weighted PSSD
    for j = 1:length(contourLvlsNum)
        % Obtain contour info, padarray makes contourc output more predictable (-> always closed contours)
        C = contourc([min(L1{ii})-1e5*eps L1{ii} max(L1{ii})+1e5*eps],...
            [min(L2{ii})-1e5*eps L2{ii} max(L2{ii})+1e5*eps],...
            padarray((q_cell{ii}.')./max(q_cell{ii}(:)),[1 1]),...
            [contourLvlsNum(j) contourLvlsNum(j)]);
        % Fill the contours
        fillContours(C,opacity,colorPSSDs(ii,:)+(j>1)*(1-colorPSSDs(ii,:))*(1/5)^(j-1),colorPSSDs(ii,:),'-',plane,offset);
    end
    % Number-weighted average dimensions
    if ~isempty(numAvgL_cell) && ~(length(numAvgL_cell)<length(q_cell))
        plot(numAvgL_cell{ii}(2),numAvgL_cell{ii}(1),...
            'Marker','o','MarkerSize',4,...
            'MarkerEdgeColor',zeros(1,3),'MarkerFaceColor',zeros(1,3),...
            'LineStyle','none','LineWidth',1);
    end
end

axis(ax,numAxisLimits);
box(ax,'on');
grid(ax,'on');
daspect([1 1 1]);

xlabel(ax,['{\it L}_1 [',char(181),'m]']);
ylabel(ax,['{\it L}_2 [',char(181),'m]']);
zlabel(ax,['{\it L}_3 [',char(181),'m]']);

if ~defaultFigureStyle
    ax.LineWidth = 1;
    ax.FontName = 'Helvetica';
    ax.FontSize = 8;
end

% Subfunction of plotPSSD, responsible for parsing contourc output and
% filling the polygons.
    function fillContours(C,alpha,faceColor,lineColor,lineStyle,plane,offset)
        
        i = 1;
        
        while i < size(C,2)
            
            % Number of vertices that are part of this contour
            np = C(2,i);
            
            % x,y coordinates of these vertices
            vertices = C(:,i+1:i+np);
            
            % if not closed
            if ~all(vertices(:,1) == vertices(:,end))
                vertices = [vertices [min(vertices(1,:));min(vertices(2,:))]];    %#ok<AGROW>
            end
            
            % Plot on xy plane
            if isequal(plane,[1,2])
                hf = fill3(vertices(1,:),vertices(2,:),zeros(length(vertices(1,:)),1)+offset,faceColor);
                
                % Plot on yz plane
            elseif isequal(plane,[1,3])
                hf = fill3(vertices(1,:),zeros(length(vertices(2,:)),1)+offset,vertices(2,:),faceColor);
                
                % Plot on xz plane
            elseif isequal(plane,[2,3])
                hf = fill3(zeros(length(vertices(2,:)),1)+offset,vertices(1,:),vertices(2,:),faceColor);
            end
            
            set(hf,'facealpha',alpha,'linewidth',0.5,'edgecolor',lineColor,'linestyle',lineStyle);
            i = i+np+1;
        end
    end
end