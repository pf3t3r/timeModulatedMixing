clc; clear; close all;
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [5 5 40 20]);
set(0,'defaultAxesFontSize',18);

% Add subfolders to path
addpath(genpath(pwd));

load Matfiles/E_hil;
load('Matfiles/Fzt.mat');
load('Matfiles/epsHil.mat');
load('Data/Mixing_parameterization_fields.mat');
load N_is.mat;
load('M1.mat','time');
mooring = ["IC3","M1"];

%% Diffusivity

Rf = 1/6;
Kp_IC3 = Rf.*sixDayEpsHM_IC3./N2_IC3;
Kp_M1 = Rf.*sixDayEpsHM_M1./N2_M1;
% Histogram of diffusivity
ax1 = figure;
histogram(Kp_IC3);
title('Diffusivity (before binning)');

savefig('figures/main/param/' + mooring(1) + '_D_hist');
exportgraphics(ax1,'figures/main/param/' + mooring(1) + '_D_hist.png');

%% Correct Diffusivity:

Kp_IC3c = Kp_IC3;
Kp_M1c = Kp_M1;


% Leave out correction just for now

% for i=1:length(Kp_IC3(:,1))
%     for j=1:length(Kp_IC3(1,:))
%         if Kp_IC3(i,j) >= 1e-2
%             Kp_IC3c(i,j) = NaN;
%         end
%     end
% end
% 
% for i=1:length(Kp_M1(:,1))
%     for j=1:length(Kp_M1(1,:))
%         if Kp_M1(i,j) >= 1e-2
%             Kp_M1c(i,j) = NaN;
%         end
%     end
% end

%% Diffusivity: IC3 (not binned)
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [5 5 35 20]);

ax1a = figure;
contourf(X2,Y2,log10(Kp_IC3c),-2.4:-0.1:-6,'LineColor','none');
datetick('x','yyyy mmm','keeplimits');
ylabel('Depth [m]');
% cmocean('tempo');
colormap(flipud(cbrewer2('Spectral')));
c = colorbar;
c.Label.String = 'log_{10}(D_{hil}) [log m^2 s^{-1}]';
ylim([-inf -100]);
% title('IC3: Diffusivity');
exportgraphics(ax1a,'figures/main/param/' + mooring(1) + '_Kp.png');

%% Diffusivity: M1 (not binned)

ax1b = figure;
contourf(X2a,Y2a,log10(Kp_M1c),-2.4:-0.1:-6,'LineColor','none');
datetick('x','yyyy mmm','keeplimits');
ylabel('Depth [m]');
% cmocean('tempo');
colormap(flipud(cbrewer2('Spectral')));
c = colorbar;
c.Label.String = 'log_{10}(D_{hil}) [log m^2 s^{-1}]';
ylim([-inf -60]);
% title('M1: Diffusivity');
exportgraphics(ax1b,'figures/main/param/' + mooring(2) + '_Kp.png');


set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [5 5 40 20]);

%% Diffusivity: binning optimisation

% Original optimisation: tuned according to histogram
%edges = [0 4e-5 8e-5 1.1e-4 1.8e-4 2.6e-4 4e-4 5.5e-4 8e-4 1.2e-3 2.6e-3 5e-3 2.5e-2 5e-2 7];

% New optimisation: put values <1e-5 and values >1e-2 in one bin
% respectively
edges = [0 1e-5 3e-5 5e-5 0.7e-4 1e-4 1.4e-4 2e-4 3.5e-4 6e-4 1e-3 2e-3 4e-3 1e-2 1e-1];

binnedDiffusivity_IC3 = discretize(Kp_IC3c,edges);
binnedDiffusivity_M1 = discretize(Kp_M1c,edges);

ax2 = figure;
histogram(binnedDiffusivity_IC3);
title('IC3: Diffusivity (after binning)');
% savefig('figures/main/param/' + mooring(1) + '_D_histBinned');
exportgraphics(ax2,'figures/main/param/' + mooring(1) + '_D_histBinned.png');

ax3 = figure;
histogram(binnedDiffusivity_M1);
title('M1: Diffusivity (after binning)');
% savefig('figures/main/param/' + mooring(2) + '_D_histBinned');
exportgraphics(ax3,'figures/main/param/' + mooring(2) + '_D_histBinned.png');

%% IC3: Diffusivity, binned

% undo thin bar temporarily
% set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [5 5 40 10]);

NCbar = length(edges)-1;

ax4 = figure;
contourf(X2,Y2,binnedDiffusivity_IC3,'LineColor','none');
datetick('x','yyyy mmm','keeplimits');
cmap = colormap(flipud(cbrewer2('Spectral',NCbar)));
c = colorbar;
c.Ticks = 1:(1-1/NCbar):length(edges);
c.TickLabels = {num2str(edges(1)), num2str(edges(2)), num2str(edges(3)), num2str(edges(4)), ... 
    num2str(edges(5)), num2str(edges(6)), num2str(edges(7)), num2str(edges(8)), num2str(edges(9)), ...
    num2str(edges(10)), num2str(edges(11)), num2str(edges(12)), num2str(edges(13)), num2str(edges(14)), num2str(edges(15))};
c.Label.String = 'K_p [m^2 s^{-1}]';
xlabel('time');
ylabel('depth [m]');
% title('IC3: Diffusivity');

% savefig('figures/main/param/' + mooring(1) + '_diffusivityBinned');
exportgraphics(ax4,'figures/main/param/' + mooring(1) + '_diffusivityBinned.png');

%% M1: Diffusivity, binned

ax5 = figure;
contourf(X2a,Y2a,binnedDiffusivity_M1,'LineColor','none');
datetick('x','yyyy mmm','keeplimits');
cmap = colormap(flipud(cbrewer2('Spectral',NCbar)));
c = colorbar;
c.Ticks = 1:(1-1/NCbar):length(edges);
c.TickLabels = {num2str(edges(1)), num2str(edges(2)), num2str(edges(3)), num2str(edges(4)), ... 
    num2str(edges(5)), num2str(edges(6)), num2str(edges(7)), num2str(edges(8)), num2str(edges(9)), ...
    num2str(edges(10)), num2str(edges(11)), num2str(edges(12)), num2str(edges(13)), num2str(edges(14)), num2str(edges(15))};
c.Label.String = 'K_p [m^2 s^{-1}]';
% % c = colorbar;
% colormap(flipud(cbrewer2('Spectral',NCbar)));
% c = colorbar('YTick',1:15,'TickLabels',{num2str(edges(1)), num2str(edges(2)), num2str(edges(3)), num2str(edges(4)), ... 
%     num2str(edges(5)), num2str(edges(6)), num2str(edges(7)), num2str(edges(8)), num2str(edges(9)), ...
%     num2str(edges(10)), num2str(edges(11)), num2str(edges(12)), num2str(edges(13)), num2str(edges(14)), num2str(edges(15))});
% c.Label.String = 'Diffusivity [m^2 s^{-1}]';
xlabel('time');
ylabel('depth [m]');
% title('M1: Diffusivity');

% savefig(ax5,'figures/main/param/'+mooring(2)+'_diffusivityBinned','compact');
exportgraphics(ax5,'figures/main/param/' + mooring(2) + '_diffusivityBinned.png');

%% Close-up of time and depth

% return to default
% set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [5 5 40 20]);

depthToPlotFrom = 30;
t1 = 1;
t2 = 26000;
t3 = 30000;
t4 = 50000;

%% IC3: Zoomed-in plot

ax6 = figure;
sgtitle('IC3: Diffusivity (zoomed-in)');

subplot(1,2,1)
contourf(X2(depthToPlotFrom:end,t1:t2),Y2(depthToPlotFrom:end,t1:t2),binnedDiffusivity_IC3(depthToPlotFrom:end,t1:t2),'LineColor','none');
cmap = colormap(flipud(cbrewer2('Spectral',NCbar)));
c = colorbar;
c.Ticks = 1:(1-1/NCbar):length(edges);
c.TickLabels = {num2str(edges(1)), num2str(edges(2)), num2str(edges(3)), num2str(edges(4)), ... 
    num2str(edges(5)), num2str(edges(6)), num2str(edges(7)), num2str(edges(8)), num2str(edges(9)), ...
    num2str(edges(10)), num2str(edges(11)), num2str(edges(12)), num2str(edges(13)), num2str(edges(14)), num2str(edges(15))};
c.Label.String = 'K_p [m^2 s^{-1}]';
datetick('x','yyyy mmm','keeplimits');
xlabel('time');
ylabel('depth [m]');

subplot(1,2,2)
contourf(X2(depthToPlotFrom:end,t3:t4),Y2(depthToPlotFrom:end,t3:t4),binnedDiffusivity_IC3(depthToPlotFrom:end,t3:t4),'LineColor','none');
datetick('x','yyyy mmm','keeplimits');
cmap = colormap(flipud(cbrewer2('Spectral',NCbar)));
c = colorbar;
c.Ticks = 1:(1-1/NCbar):length(edges);
c.TickLabels = {num2str(edges(1)), num2str(edges(2)), num2str(edges(3)), num2str(edges(4)), ... 
    num2str(edges(5)), num2str(edges(6)), num2str(edges(7)), num2str(edges(8)), num2str(edges(9)), ...
    num2str(edges(10)), num2str(edges(11)), num2str(edges(12)), num2str(edges(13)), num2str(edges(14)), num2str(edges(15))};
c.Label.String = 'K_p [m^2 s^{-1}]';
xlabel('time');
ylabel('depth [m]');

% savefig('figures/main/param/' + mooring(1) + '_diffusivityBinned_Zoom');
exportgraphics(ax6,'figures/main/param/' + mooring(1) + '_diffusivityBinned_Zoom.png');

%% M1: zoomed-in plot

ax7 = figure;
sgtitle('M1: Diffusivity (zoomed-in)');

subplot(1,2,1)
contourf(X2a(depthToPlotFrom:end,t1:t2),Y2a(depthToPlotFrom:end,t1:t2),binnedDiffusivity_M1(depthToPlotFrom:end,t1:t2),'LineColor','none');
datetick('x','yyyy mmm','keeplimits');
cmap = colormap(flipud(cbrewer2('Spectral',NCbar)));
c = colorbar;
c.Ticks = 1:(1-1/NCbar):length(edges);
c.TickLabels = {num2str(edges(1)), num2str(edges(2)), num2str(edges(3)), num2str(edges(4)), ... 
    num2str(edges(5)), num2str(edges(6)), num2str(edges(7)), num2str(edges(8)), num2str(edges(9)), ...
    num2str(edges(10)), num2str(edges(11)), num2str(edges(12)), num2str(edges(13)), num2str(edges(14)), num2str(edges(15))};
c.Label.String = 'K_p [m^2 s^{-1}]';
% c = colorbar;
% colormap(flipud(cbrewer2('Spectral',NCbar)));
% c.Label.String = '\epsilon_{hil} [W kg^{-1}]';
% c = colorbar('YTick',1:15,'TickLabels',{num2str(edges(1)), num2str(edges(2)), num2str(edges(3)), num2str(edges(4)), ... 
%     num2str(edges(5)), num2str(edges(6)), num2str(edges(7)), num2str(edges(8)), num2str(edges(9)), ...
%     num2str(edges(10)), num2str(edges(11)), num2str(edges(12)), num2str(edges(13)), num2str(edges(14)), num2str(edges(15))});
xlabel('time');
ylabel('depth [m]');

subplot(1,2,2)
contourf(X2a(depthToPlotFrom:end,t3:t4),Y2a(depthToPlotFrom:end,t3:t4),binnedDiffusivity_M1(depthToPlotFrom:end,t3:t4),'LineColor','none');
datetick('x','yyyy mmm','keeplimits');
cmap = colormap(flipud(cbrewer2('Spectral',NCbar)));
c = colorbar;
c.Ticks = 1:(1-1/NCbar):length(edges);
c.TickLabels = {num2str(edges(1)), num2str(edges(2)), num2str(edges(3)), num2str(edges(4)), ... 
    num2str(edges(5)), num2str(edges(6)), num2str(edges(7)), num2str(edges(8)), num2str(edges(9)), ...
    num2str(edges(10)), num2str(edges(11)), num2str(edges(12)), num2str(edges(13)), num2str(edges(14)), num2str(edges(15))};
c.Label.String = 'K_p [m^2 s^{-1}]';
% c = colorbar;
% colormap(flipud(cbrewer2('Spectral',NCbar)));
% c.Label.String = '\epsilon_{hil} [W kg^{-1}]';
% c = colorbar('YTick',1:15,'TickLabels',{num2str(edges(1)), num2str(edges(2)), num2str(edges(3)), num2str(edges(4)), ... 
%     num2str(edges(5)), num2str(edges(6)), num2str(edges(7)), num2str(edges(8)), num2str(edges(9)), ...
%     num2str(edges(10)), num2str(edges(11)), num2str(edges(12)), num2str(edges(13)), num2str(edges(14)), num2str(edges(15))});
xlabel('time');
ylabel('depth [m]');

% savefig('figures/main/param/' + mooring(2) + '_diffusivityBinned_Zoom');
exportgraphics(ax7,'figures/main/param/' + mooring(2) + '_diffusivityBinned_Zoom.png');

%% clear
clear depthToPlotFrom t1 t2 t3 t4;

%% Kurtosis for Kp at IC3 and M1

% Not sure if this is correct... very high kurtosis
% is is even appropriate to find kurtosis of a derived value such as this?

for i=1:length(zmid_IC3(:,1))
    kurtKp_IC3(i) = kurtosis(Kp_IC3c(i,:));
    kurtKp_M1(i) = kurtosis(Kp_M1c(i,:));
end

ax7a = figure;
subplot(1,2,1)
plot(kurtKp_IC3,zmid_IC3(:,1),'DisplayName','Kurt(K_p(z,t))');
hold on
xline(3,':','DisplayName','Kurtosis = 3');
hold off
legend('Location','best');
xlabel('Kurtosis(K_p(z,t))');
ylabel('Depth [m]');
title('IC3 (2014-2020)');

subplot(1,2,2)
plot(kurtKp_M1,zmid_IC3(:,1),'DisplayName','Kurt(K_p(z,t))');
hold on
xline(3,':','DisplayName','Kurtosis = 3');
hold off
legend('Location','best');
xlabel('Kurtosis(K_p(z,t))');
ylabel('Depth [m]');
title('M1 (2014-2020)');

savefig('figures/main/param/_kurtosisOfKpAcrossDepth');
exportgraphics(ax7a,'figures/main/param/_kurtosisOfKpAcrossDepth.png');



%% Calculate Flux(z,t) for IC3 and M1

dz_IC3 = cat(1,10*ones(5,1),25*ones(6,1),50*ones(30,1));
dz_M1 = cat(1,10*ones(5,1),25*ones(6,1),50*ones(31,1));
% load('SA_CT_interpolated.mat','CTmid_IC3_i','CTmid_M1_i','CT_M1_i');
load('SA_CT_stabilised.mat','CT_M1_is','CT_IC3_is');

% We should use the original non-stabilised values for CT
load('SA_CT_interpolated.mat','CT_IC3_i','CT_M1_i');
CT_M1 = CT_M1_i;
CT_IC3 = CT_IC3_i;

cp0 = 3992;

for j=1:length(time)
    for i=1:length(dz_M1)
        dThetadz_M1(i,j) = (CT_M1(i+1,j) - CT_M1(i,j))/(dz_M1(i)); 
    end
    for i=1:length(dz_IC3)
        dThetadz_IC3(i,j) = (CT_IC3(i+1,j) - CT_IC3(i,j))/(dz_IC3(i));
    end
end

% dThetadz_IC3(isinf(dThetadz_IC3)) = 0;
flux_IC3 = cp0.*rhoIC3.*Kp_IC3c.*dThetadz_IC3;

% dThetadz_M1(isinf(dThetadz_M1)) = 0;
flux_M1 = cp0.*rhoM1.*Kp_M1c.*dThetadz_M1;



%% figure for temp gradient only

% Not needed, just for checking

figure
contourf(X2,Y2,dThetadz_IC3,'LineColor','none');
c = colorbar;
c.Label.String = 'd\Theta/dz [Cm^{-1}]';
colormap(flipud(cbrewer2('RdBu')));
title('IC3: temperature gradient');

figure
contourf(X2a,Y2a,dThetadz_M1,'LineColor','none');
c = colorbar;
c.Label.String = 'd\Theta/dz [Cm^{-1}]';
colormap(flipud(cbrewer2('RdBu')));
title('M1: temperature gradient');


%% Conservative Temperature

load('SA_CT_interpolated.mat','CTmid_IC3_i','CTmid_M1_i');

axAA = figure;
contourf(X2,Y2,CTmid_IC3_i,'LineColor','none');
datetick('x','yyyy mmm','keeplimits');
c = colorbar;
c.Label.String = '\Theta [K]';
colormap(flipud(cbrewer2('Spectral')));
title('IC3: Conservative Temperature');
exportgraphics(axAA,'Figures/Main/Param/CT_IC3.png');

axAB = figure;
contourf(X2a,Y2a,CTmid_M1_i,'LineColor','none');
datetick('x','yyyy mmm','keeplimits');
c = colorbar;
c.Label.String = '\Theta [K]';
colormap(flipud(cbrewer2('Spectral')));
title('M1: Conservative Temperature');
exportgraphics(axAB,'Figures/Main/Param/CT_M1.png');

%% Absolute Salinity
load('SA_CT_interpolated.mat','SA_IC3_i','SAmid_IC3_i','SA_M1_i','SAmid_M1_i');

axAC = figure;
contourf(X2,Y2,SAmid_IC3_i,'LineColor','none');
datetick('x','yyyy mmm','keeplimits');
c = colorbar;
c.Label.String = 'S_A [g kg^{-1}]';
colormap(flipud(cbrewer2('Spectral')));
title('IC3: Absolute Salinity');
exportgraphics(axAC,'Figures/Main/Param/SA_IC3.png');

axAD =figure;
contourf(X2a,Y2a,SAmid_M1_i,'LineColor','none');
datetick('x','yyyy mmm','keeplimits');
c = colorbar;
c.Label.String = 'S_A [g kg^{-1}]';
colormap(flipud(cbrewer2('Spectral')));
title('M1: Absolute Salinity');
exportgraphics(axAD,'Figures/Main/Param/SA_M1.png');

%% Salinity Gradient

for j=1:length(time)
    for i=1:length(dz_M1)
        dSAdz_M1(i,j) = (SA_M1_i(i+1,j) - SA_M1_i(i,j))/(dz_M1(i)); 
    end
    for i=1:length(dz_IC3)
        dSAdz_IC3(i,j) = (SA_IC3_i(i+1,j) - SA_IC3_i(i,j))/(dz_IC3(i));
    end
end

axAE = figure;
contourf(X2,Y2,1000*dSAdz_IC3,'LineColor','none');
datetick('x','yyyy mmm','keeplimits');
c = colorbar;
c.Label.String = '\Delta S_A [mg kg^{-1}]';
colormap(flipud(cbrewer2('RdBu')));
title('IC3: \Delta S_A');
exportgraphics(axAE,'Figures/Main/Param/deltaSA_IC3.png');

axAF =figure;
contourf(X2a,Y2a,1000*dSAdz_M1,'LineColor','none');
datetick('x','yyyy mmm','keeplimits');
c = colorbar;
c.Label.String = '\Delta S_A [mg kg^{-1}]';
colormap(flipud(cbrewer2('RdBu')));
title('M1: \Delta S_A');
exportgraphics(axAF,'Figures/Main/Param/deltaSA_M1.png');

%% Density

load('epsHil.mat','rhoIC3','rhoM1');

axAG = figure;
contourf(X2,Y2,rhoIC3,'LineColor','none');
datetick('x','yyyy mmm','keeplimits');
c = colorbar;
c.Label.String = '\rho [kg m^{-3}]';
colormap(flipud(cbrewer2('Spectral')));
title('IC3: Density');
exportgraphics(axAG,'Figures/Main/Param/rho_IC3.png');

axAH = figure;
contourf(X2a,Y2a,rhoM1,'LineColor','none');
datetick('x','yyyy mmm','keeplimits');
c = colorbar;
c.Label.String = '\rho [kg m^{-3}]';
colormap(flipud(cbrewer2('Spectral')));
title('M1: Density');
exportgraphics(axAH,'Figures/Main/Param/rho_M1.png');


%% Flux, no bin (IC3)
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [5 5 35 20]);

figure
test = [-200 -50 -10:0.1:10 50 200];
ax7b = figure;
contourf(X2,Y2,flux_IC3,test,'LineColor','none');
datetick('x','yy mmm','keeplimits','keepticks');
ylabel('Depth [m]');
% cmocean('tempo');
colormap(flipud(cbrewer2('RdBu')));
c = colorbar;
c.Label.String = 'Flux [W m^{-2}]';
% title('IC3: Diffusivity');
exportgraphics(ax7b,'figures/main/param/' + mooring(1) + '_Flux.png');


%% Optimise flux binning

ax8 = figure;
histogram(flux_IC3);
title('IC3: Flux (before binning)');

savefig(ax8,'figures/main/param/' + mooring(1) + '_FluxHistogram');
exportgraphics(ax8,'figures/main/param/' + mooring(1) + '_FluxHistogram.png');

ax8a = figure;
histogram(flux_M1);
title('M1: Flux (before binning)');

% edges_IC3 = [-200 -10 -3 -2 -1.6 -1.3 -1.1 -0.85 -0.7 -0.6 -0.45 -0.3 0 1 max(max(flux_IC3))];
edges_IC3 = [-200 -20 -1 -0.8 -0.5 -0.2 -0.1 0 0.1 0.2 0.5 0.8 1 20 200];

% edges_M1 = [-50 -10 -3 -2 -1.6 -1.3 -1.1 -0.85 -0.7 -0.6  -0.45 -0.3 0 1 max(max(flux_M1))];
edges_M1 = edges_IC3;

binnedFlux_IC3 = discretize(flux_IC3,edges_IC3);
binnedFlux_M1 = discretize(flux_M1,edges_M1);

ax9 = figure;
histogram(binnedFlux_IC3);
title('IC3: Flux (after binning)');
savefig(ax9,'figures/main/param/' + mooring(1) + '_FluxHistogramBinned');
exportgraphics(ax9,'figures/main/param/' + mooring(1) + '_FluxHistogramBinned.png');

ax9a = figure;
histogram(binnedFlux_M1);
title('M1: Flux (after binning)');

%% IC3: Flux(z,t)

ax10 = figure;
contourf(X2,Y2,binnedFlux_IC3,'LineStyle','none');
datetick('x','yyyy mmm','keeplimits');
cmap = colormap(cbrewer2('RdBu',NCbar));
c = colorbar;
c.Ticks = 1:(1-1/NCbar):length(edges_IC3);
c.TickLabels = {num2str(edges_IC3(1)), num2str(edges_IC3(2)), num2str(edges_IC3(3)), num2str(edges_IC3(4)), ... 
    num2str(edges_IC3(5)), num2str(edges_IC3(6)), num2str(edges_IC3(7)), num2str(edges_IC3(8)), num2str(edges_IC3(9)), ...
    num2str(edges_IC3(10)), num2str(edges_IC3(11)), num2str(edges_IC3(12)), num2str(edges_IC3(13)), num2str(edges_IC3(14)), num2str(edges_IC3(15))};
% colormap(cbrewer2('RdBu',NCbar));
% c = colorbar('YTick',1:15,'TickLabels',{num2str(edges_IC3(1)), num2str(edges_IC3(2)), num2str(edges_IC3(3)), num2str(edges_IC3(4)), ... 
%     num2str(edges_IC3(5)), num2str(edges_IC3(6)), num2str(edges_IC3(7)), num2str(edges_IC3(8)), num2str(edges_IC3(9)), ...
%     num2str(edges_IC3(10)), num2str(edges_IC3(11)), num2str(edges_IC3(12)), num2str(edges_IC3(13)), num2str(edges_IC3(14)), num2str(edges_IC3(15))});
c.Label.String = 'Flux [W m^{-2}]';
% xlabel('time');
ylabel('Depth [m]');
ylim([-inf -100]);
% title('IC3: flux');

% savefig(ax10,'figures/main/param/' + mooring(1) + '_FluxBinned');
exportgraphics(ax10,'figures/main/param/' + mooring(1) + '_FluxBinned.png');


%% M1: Flux(z,t)

ax11 = figure;
contourf(X2a,Y2a,binnedFlux_M1,'LineStyle','none');
datetick('x','yyyy mmm','keeplimits');
cmap = colormap(cbrewer2('RdBu',NCbar));
c = colorbar;
c.Ticks = 1:(1-1/NCbar):length(edges_M1);
c.TickLabels = {num2str(edges_M1(1)), num2str(edges_M1(2)), num2str(edges_M1(3)), num2str(edges_M1(4)), ... 
    num2str(edges_M1(5)), num2str(edges_M1(6)), num2str(edges_M1(7)), num2str(edges_M1(8)), num2str(edges_M1(9)), ...
    num2str(edges_M1(10)), num2str(edges_M1(11)), num2str(edges_M1(12)), num2str(edges_M1(13)), num2str(edges_M1(14)), num2str(edges_M1(15))};
% colormap(cbrewer2('RdBu',NCbar));
% c = colorbar('YTick',1:15,'TickLabels',{num2str(edges_M1(1)), num2str(edges_M1(2)), num2str(edges_M1(3)), num2str(edges_M1(4)), ... 
%     num2str(edges_M1(5)), num2str(edges_M1(6)), num2str(edges_M1(7)), num2str(edges_M1(8)), num2str(edges_M1(9)), ...
%     num2str(edges_M1(10)), num2str(edges_M1(11)), num2str(edges_M1(12)), num2str(edges_M1(13)), num2str(edges_M1(14)), num2str(edges_M1(15))});
c.Label.String = 'Flux [W m^{-2}]';
% xlabel('time');
ylabel('Depth [m]');
ylim([-inf -60]);
% title('M1: Flux');

% savefig(ax11,'figures/main/param/' + mooring(2) + '_FluxBinned');
exportgraphics(ax11,'figures/main/param/' + mooring(2) + '_FluxBinned.png');
