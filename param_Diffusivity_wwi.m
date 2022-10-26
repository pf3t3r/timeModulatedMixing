clc; clear; close all;
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [5 5 40 22]);
set(0,'defaultAxesFontSize',18);

% Add subfolders to path
addpath(genpath(pwd));

load E_wwi.mat;
load Fzt_wwi.mat;
load epsWwi.mat;
load('Mixing_parameterization_fields.mat');
% load N2v2.mat;
load IC3.mat;
% load N2-CT.mat;
load N_is.mat;
load('epsHil.mat','rhoIC3');
mooring = "IC3";

%% Diffusivity

Rf = 1/6;
D_wwi = Rf.*sixDayEpsWwiM./N2_IC3;

% Histogram of diffusivity
ax1 = figure;
histogram(D_wwi);
title('Diffusivity (before binning)');

savefig('figures/main/param/' + mooring + '_D_hist');
exportgraphics(ax1,'figures/main/param/' + mooring + '_D_hist.png');

save Matfiles/D_wwi.mat D_wwi;

%% Diffusivity: binning optimisation

% New optimisation: put values <1e-5 and values >1e-2 in one bin
% respectively
edges = [0 1e-5 3e-5 5e-5 0.7e-4 1e-4 1.4e-4 2e-4 3.5e-4 6e-4 1e-3 2e-3 4e-3 1e-2 7];

binnedDiffusivity = discretize(D_wwi,edges);

ax2 = figure;
histogram(binnedDiffusivity);
title('Diffusivity (after binning)');

savefig('figures/main/param/' + mooring + '_D_histBinned');
exportgraphics(ax2,'figures/main/param/' + mooring + '_D_histBinned.png');

%% Diffusivity: binned plot

NCbar = length(edges)-1;

ax3 = figure;
contourf(X2,Y2,binnedDiffusivity,'LineColor','none');
datetick('x','yyyy mmm','keeplimits');
c = colorbar;
colormap(flipud(cbrewer2('Spectral',NCbar)));
c.Label.String = '\epsilon_{hil} [W kg^{-1}]';
c = colorbar('YTick',1:15,'TickLabels',{num2str(edges(1)), num2str(edges(2)), num2str(edges(3)), num2str(edges(4)), ... 
    num2str(edges(5)), num2str(edges(6)), num2str(edges(7)), num2str(edges(8)), num2str(edges(9)), ...
    num2str(edges(10)), num2str(edges(11)), num2str(edges(12)), num2str(edges(13)), num2str(edges(14)), num2str(edges(15))});
xlabel('time');
ylabel('depth [m]');
zlabel('\epsilon_{hil}');
title('binned diffusivity');

savefig('figures/main/param/' + mooring + '_diffusivityBinned');
exportgraphics(ax3,'figures/main/param/' + mooring + '_diffusivityBinned.png');


%% Diffusivity: IC3 (not binned)

ax3a = figure;
contourf(X2,Y2,log10(D_wwi),-2.4:-0.1:-6,'LineColor','none');
datetick('x','yyyy mmm','keeplimits');
ylabel('Depth [m]');
% cmocean('tempo');
colormap(flipud(cbrewer2('Spectral')));
ylim([-inf -100]);
c = colorbar;
c.Label.String = 'log_{10}(D_{hil}) [log m^2 s^{-1}]';

exportgraphics(ax3a,'figures/main/param/' + mooring + '_D_Wwi.png');

%% Zoomed-in contour plot

depthToPlotFrom = 33;
t1 = 9000;
t2 = 10000;
t3 = 30000;
t4 = 31000;

ax3a = figure;
sgtitle('Diffusivity: zoomed-in');

subplot(1,2,1)
contourf(X2(depthToPlotFrom:end,t1:t2),Y2(depthToPlotFrom:end,t1:t2),binnedDiffusivity(depthToPlotFrom:end,t1:t2),'LineColor','none');
datetick('x','yyyy mmm','keeplimits');
c = colorbar;
colormap(flipud(cbrewer2('Spectral',NCbar)));
c.Label.String = '\epsilon_{hil} [W kg^{-1}]';
c = colorbar('YTick',1:15,'TickLabels',{num2str(edges(1)), num2str(edges(2)), num2str(edges(3)), num2str(edges(4)), ... 
    num2str(edges(5)), num2str(edges(6)), num2str(edges(7)), num2str(edges(8)), num2str(edges(9)), ...
    num2str(edges(10)), num2str(edges(11)), num2str(edges(12)), num2str(edges(13)), num2str(edges(14)), num2str(edges(15))});
xlabel('time');
ylabel('depth [m]');
zlabel('\epsilon_{hil}');

subplot(1,2,2)
contourf(X2(depthToPlotFrom:end,t3:t4),Y2(depthToPlotFrom:end,t3:t4),binnedDiffusivity(depthToPlotFrom:end,t3:t4),'LineColor','none');
datetick('x','yyyy mmm','keeplimits');
c = colorbar;
colormap(flipud(cbrewer2('Spectral',NCbar)));
c.Label.String = '\epsilon_{hil} [W kg^{-1}]';
c = colorbar('YTick',1:15,'TickLabels',{num2str(edges(1)), num2str(edges(2)), num2str(edges(3)), num2str(edges(4)), ... 
    num2str(edges(5)), num2str(edges(6)), num2str(edges(7)), num2str(edges(8)), num2str(edges(9)), ...
    num2str(edges(10)), num2str(edges(11)), num2str(edges(12)), num2str(edges(13)), num2str(edges(14)), num2str(edges(15))});
xlabel('time');
ylabel('depth [m]');
zlabel('\epsilon_{hil}');

clear depthToPlotFrom t1 t2 t3 t4;

% savefig('figures/main/param/' + mooring + '_diffusivityBinned_Zoom');
exportgraphics(ax3a,'figures/main/param/' + mooring + '_diffusivityBinned_Zoom.png');

%% Flux


dz_IC3 = cat(1,10*ones(5,1),25*ones(6,1),50*ones(30,1));
% We should use the original non-stabilised values for CT
load('SA_CT_interpolated.mat','CT_IC3_i');
CT_IC3 = CT_IC3_i;
cp0 = 3992;


for j=1:length(time)
    for i=1:length(dz_IC3)
        dThetadz_IC3(i,j) = (CT_IC3(i+1,j) - CT_IC3(i,j))/(dz_IC3(i));
    end
end

flux_IC3 = cp0.*rhoIC3.*D_wwi.*dThetadz_IC3;

%% Flux plot

range = [-200 -20 -1 -0.8 -0.5 -0.2 -0.1 0 0.1 0.2 0.5 0.8 1 20 200];
flux_IC3_b = discretize(flux_IC3,range);

ax4 = figure;
contourf(X2,Y2,flux_IC3_b,'LineColor','none');
datetick('x','yyyy mmm','keeplimits');
ylabel('Depth [m]');
colormap(cbrewer2('RdBu'));
c = colorbar;
c.Ticks = 1:(1-1/NCbar):length(range);
c.TickLabels = {num2str(range(1)), num2str(range(2)), num2str(range(3)), num2str(range(4)), ... 
    num2str(range(5)), num2str(range(6)), num2str(range(7)), num2str(range(8)), num2str(range(9)), ...
    num2str(range(10)), num2str(range(11)), num2str(range(12)), num2str(range(13)), num2str(range(14)), num2str(range(15))};
c.Label.String = 'Flux [W m^{-2}]';
ylim([-inf -100]);
% title('IC3: Diffusivity');
exportgraphics(ax4,'figures/main/param/' + mooring(1) + '_FluxWwi.png');

%% Optimise flux binning

% ax4 = figure;
% histogram(flux);
% title('Flux (before binning)');
% 
% savefig('figures/main/param/' + mooring + '_FluxHistogram');
% exportgraphics(ax4,'figures/main/param/' + mooring + '_FluxHistogram.png');
% 
% edges = [-10 -5 -3 -2 -1 -0.9 -0.8 -0.7 -0.6  -0.5 -0.3 -0.2 -0.1  0 0.0001];
% 
% binnedFlux = discretize(flux,edges);
% 
% ax5 = figure;
% histogram(binnedFlux);
% title('Flux (after binning)');
% savefig('figures/main/param/' + mooring + '_FluxHistogramBinned');
% exportgraphics(ax5,'figures/main/param/' + mooring + '_FluxHistogramBinned.png');

%% Plot the flux
% 
% ax6 = figure;
% %contourf(X2,Y2,flux','LineStyle','none');
% contourf(X2,Y2,binnedFlux','LineStyle','none');
% datetick('x','yyyy mmm','keeplimits');
% c = colorbar;
% colormap(flipud(cbrewer2('Spectral',NCbar)));
% c.Label.String = '\epsilon_{hil} [W kg^{-1}]';
% c = colorbar('YTick',1:15,'TickLabels',{num2str(edges(1)), num2str(edges(2)), num2str(edges(3)), num2str(edges(4)), ... 
%     num2str(edges(5)), num2str(edges(6)), num2str(edges(7)), num2str(edges(8)), num2str(edges(9)), ...
%     num2str(edges(10)), num2str(edges(11)), num2str(edges(12)), num2str(edges(13)), num2str(edges(14)), num2str(edges(15))});
% xlabel('time');
% ylabel('depth [m]');
% zlabel('\epsilon_{hil}');
% title('binned flux');
% 
% savefig('figures/main/param/' + mooring + '_FluxBinned');
% exportgraphics(ax6,'figures/main/param/' + mooring + '_FluxBinned.png');
