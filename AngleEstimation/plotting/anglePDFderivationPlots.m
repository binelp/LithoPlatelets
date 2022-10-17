%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ETH Zurich, Switzerland
% Separation Processes Laboratory
%
% Project:  Lithoplatelets
% Year:     2022
% MATLAB:   R2019b, Windows 64bit
% Authors:  Anna Jaeggi (AJ)
%
% Purpose:
% Generate plots for section about derivation of distribution of phi, theta
% in the supporting information
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all;
%% Distribution of theta alone
x = 0:0.01:pi/2;
fig = figure;
plot(x, sin(x), 'color', [0 0 0])
xlim([0,pi/2])

xlabel("\theta")
ylabel("f_{\theta}(\theta)")

xticks([0 pi/4 pi/2])
xticklabels({'0','1/4 \pi','1/2 \pi'})

typeSetPlot(fig, 'half', 'noPADDING')
set(fig,'Position',[360 485 156 133]);

set(gca,'fontsize',8)

print(gcf,'deriv1.pdf','-dpdf');

%% Alpha as a function of beta
fig = figure;
hold on;
ii = 1;
colors = [150, 25, 77;
    193, 75, 46.5000;
    236, 125,16;
    145.5, 155.5, 102.5;
    55, 186, 189
]./255;
for theta = [0, pi/8, pi/4, 3/8*pi , pi/2]
    beta = 0:0.01:2*pi;
    alpha = asin(sin(theta)*sin(beta));
    h(ii) = plot(beta, alpha, 'color', colors(ii,:));
    yline(theta, 'linestyle',':', 'color', colors(ii,:))
    yline(-theta, 'linestyle',':', 'color', colors(ii,:))
    ii = ii+1;
end
xlim([0,pi*2])
xlabel("\beta")
ylabel("\alpha")

xticks([0 pi/2 pi 3/2*pi 2*pi])
xticklabels({'0','1/2 \pi','\pi','3/2 \pi','2 \pi'})

yticks([-pi/2 -pi/4 0 pi/4 pi/2 ])
yticklabels({' -1/2 \pi','-1/4 \pi','0','1/4 \pi','1/2 \pi'})

typeSetPlot(fig, 'full', 'noPADDING')
set(fig,'Position',[360 313 312 305]);

set(gca,'fontsize',8)

lgd = legend(h,"0", "1/8 \pi", "1/4 \pi", "3/8 \pi", "1/2 \pi");
lgd.Title.String = '\theta';
lgd.Title.FontSize = 8;

print(gcf,'deriv2.pdf','-dpdf');
%% Derivative of alpha with respect to beta

fig = figure;
hold on;
ii = 1;
colors = [150, 25, 77;
    193, 75, 46.5000;
    236, 125,16;
    145.5, 155.5, 102.5;
    55, 186, 189
]./255;

for theta = [pi/8, pi/4, 3/8*pi]
    phi = (pi/2 - theta + eps):0.01:pi/2;
    f = - sin(phi)./(sin(theta).*sqrt(1-cos(phi).^2./sin(theta).^2));
    h(ii) = plot(phi, f, 'color', colors(ii+1,:));
    xline(pi/2 - theta, 'linestyle',':', 'color', colors(ii+1,:))
    ii = ii+1;
end
xlim([0, pi/2])
ylim([-15,0])
xlabel("\phi")
ylabel("d\beta/d\phi")

xticks([0 pi/4 pi/2 ])
xticklabels({'0','1/4 \pi','1/2 \pi'})

typeSetPlot(fig, 'full', 'noPADDING')
set(fig,'Position',[360 313 312 305]);
set(gca,'fontsize',8)

lgd = legend(h, "1/8 \pi", "1/4 \pi", "3/8 \pi", 'Location','southwest');
lgd.Title.String = '\theta';
lgd.Title.FontSize = 8;

print(gcf,'deriv3.pdf','-dpdf');

%% Half analytical, half numerical integration for binning
n =  30;
delta = (pi/2)/n;
for ii=1:n
    for jj=1:n
        a = (ii-1)*delta;
        b = (ii)*delta;
        c = (jj-1)*delta;
        d = (jj)*delta;
        
        theta_bin(ii,jj) = (a+b)/2;
        phi_bin(ii,jj) = (c+d)/2;
        intPhi = @(theta) -asin(cos(d)./sin(theta)).*sin(theta) + asin(cos(c)./sin(theta)).*sin(theta);
        intTheta = integral(intPhi,a,b);
        intTheta = real(intTheta);
        f_bin(ii,jj) = intTheta;
        
    end
end

f_bin = f_bin/(nansum(f_bin,'all')*delta^2);
nansum(f_bin,'all')*delta^2
max(f_bin,[],'all')
f_bin(30,30)

fig = figure;
hold on;

s = surf(theta_bin, phi_bin, f_bin);
view(60,15)

xticks([0 pi/4 pi/2 ])
xticklabels({'0','1/4 \pi','1/2 \pi'})

yticks([0 pi/4 pi/2 ])
yticklabels({'0','1/4 \pi','1/2 \pi'})

xlim([0, pi/2])
ylim([0, pi/2])
zlim([0,2])
xlabel('\theta')
ylabel('\phi')
zlabel("p_{\phi, \theta}(\phi, \theta)")

set(gca,'fontsize',8)

typeSetPlot(fig, 'full', 'noPADDING')
set(fig,'Position',[360 343.6667 312 274.3333]);

export_fig(gcf, '-dpng','-transparent','-r600', "deriv7.png")
%% Sampled distribution

% Obtain distribution of theta, phi from uniform distributed orientation
% samples
q = uniformSampledQuats(1000000);
physicalAngles = getPhysicalAngles(q, [120,100,30], "uniformlyRandomPhysAngles.mat");

nbins = 30;
% Binning
x = physicalAngles.flowAngleVec;
y = physicalAngles.flowAnglePlane;

% Transform to radians
x = x/180*pi;
y = y/180*pi;

% Create bins and obtain counts
Xedges=linspace(0,pi/2,nbins+1);
Yedges=linspace(0,pi/2,nbins+1);
N = histcounts2(x,y,Xedges,Yedges);

delta = (pi/2)/nbins;
[X,Y] = meshgrid(linspace(delta/2,pi/2- delta/2,nbins),...
    linspace(delta/2,pi/2- delta/2,nbins));

% Normalize
N = N/(nansum(N,'all')*delta^2);

max(N,[],'all')
N(30,30)

fig = figure;
hold on;
surf(X,Y,N')
view(3)
xlabel("theta")
ylabel("phi")
view(60,15)

xticks([0 pi/4 pi/2 ])
xticklabels({'0','1/4 \pi','1/2 \pi'})

yticks([0 pi/4 pi/2 ])
yticklabels({'0','1/4 \pi','1/2 \pi'})

xlim([0, pi/2])
ylim([0, pi/2])
zlim([0,2])
xlabel('\theta')
ylabel('\phi')
zlabel("p_{\phi, \theta}(\phi, \theta)")

set(gca,'fontsize',8)

typeSetPlot(fig, 'full', 'noPADDING')
set(fig,'Position',[360 343.6667 312 274.3333]);
export_fig(gcf, '-dpng','-transparent','-r600', "deriv5.png")

%% Marginal dist comparison

% Marginal cumulative dist
fig = figure;
histogram(y,'Normalization','cdf','DisplayStyle','stairs','edgecolor','k')
title("")
hold on
plot(Xedges, 1 - cos(Xedges),'r-')

set(gca,'fontsize',8)
xlabel("\theta")
ylabel("F_{m,\theta}(\theta)")

xticks([0 pi/4 pi/2 ])
xlim([0,pi/2])
ylim([0,1])
xticklabels({'0','1/4 \pi','1/2 \pi'})
typeSetPlot(fig, 'full', 'noPADDING')
set(fig,'Position',[360 324.3333 312 293.6667]);

legend('Empirical CDF','Analytical CDF','Location','southeast')
print(gcf,'deriv8.pdf','-dpdf');
fig = figure;
histogram(y,'Normalization','pdf','DisplayStyle','stairs','edgecolor','k')
title("")
hold on
x_values = linspace(min(x),max(x));
plot(Xedges, sin(Xedges),'r-')
xlabel("\phi")

ylabel("F_{m,\phi}(\phi)")

xticks([0 pi/4 pi/2 ])
xlim([0,pi/2])
ylim([0,1])
xticklabels({'0','1/4 \pi','1/2 \pi'})

typeSetPlot(fig, 'full', 'noPADDING')
set(fig,'Position',[360 324.3333 312 293.6667]);

legend('Empirical CDF','Analytical CDF','Location','southeast')
print(gcf,'deriv9.pdf','-dpdf');