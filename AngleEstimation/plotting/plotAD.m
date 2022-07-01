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
% Plot 3D contour plot of orientationDistributions. 
%
% Input arguments:
% - AD:         Quaternion/ angle distribution struct
% - color:      Contour surfaces color
% - offset:     Offset for marginal distributions from edge
% - alpha:      Transparency of contour surfaces
% - marginals:  Wether to plot marginal distributions or not (0 = no, 1 = 
%               only marginals, 2 = marginals and 3D dist )
% - fontsize:   fontsize
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function []=plotAD(AD, color, offset, alpha, marginals, fontsize)

AD.F=permute(AD.F,[3,2,1]);
AD.grid = AD.grid([3,2,1]);

hold on;
if marginals==0 || marginals==2
    plot3DPSSD_coloredContours({AD.F},{AD.F},AD.grid(1).y,AD.grid(2).y,AD.grid(3).y,color,alpha); hold on;
end
% Generate and plot the marginal distribution
if marginals==1 || marginals==2
    for ignoredDim=1:3
        dims=find([1,2,3]~=ignoredDim);
        AD_m={squeeze(sum(AD.F,ignoredDim))};
        plotPSSD_coloredContours_marginal(dims,AD_m,AD_m,AD.grid(dims(2)).y,AD.grid(dims(1)).y,color,alpha,AD.grid(ignoredDim).boundaries(1)-offset);
            
        box on;
        grid on;
        daspect([1 1 1]);
        view(90+50,30)
        
    end
end
xlim([AD.grid(2).boundaries(1),AD.grid(2).boundaries(end)]);
ylim([AD.grid(1).boundaries(1),AD.grid(1).boundaries(end)]);
zlim([AD.grid(3).boundaries(1),AD.grid(3).boundaries(end)]);

xlabel(sprintf('$q_%i$',3),'Interpreter','latex','fontsize',fontsize)
ylabel(sprintf('$q_%i$',4),'Interpreter','latex','fontsize',fontsize)
zlabel(sprintf('$q_%i$',2),'Interpreter','latex','fontsize',fontsize)

set(gca,'fontsize',fontsize)
end