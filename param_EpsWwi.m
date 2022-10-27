clc; clear; close all;
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [5 5 40 15]);
set(0,'defaultAxesFontSize',18);

% Add subfolders to path
addpath(genpath(pwd));

load E_wwi.mat;
% load('E_hil.mat','zq');
load('epsHil.mat','rhoIC3');
load Fzt_wwi.mat;
load('Mixing_parameterization_fields.mat');
% load N2v2.mat;
load N_is.mat;
load IC3.mat;
mooring = "IC3";

%% Compute eps_wwi(z,t)

rho = rhoIC3;

for i=1:length(time)
    eps_wwi(:,i) = E_wwi_IC3_mod_6dm(i).*verticalStructureIC3Wwi(:,i)./rho(:,i);
end
% mean_eps_wwi2_IC3 = mean(eps_wwi,'omitnan');

sixDayEpsWwiM = movmean(eps_wwi,149,2);

%% kurtosis: eps_wwi

for i=1:length(zmid_IC3(:,1))
    kurtEpsWwi(i) = kurtosis(sixDayEpsWwiM(i,:));
end

ax0 = figure;
plot(kurtEpsWwi,zmid_IC3(:,1),'DisplayName','Kurt(\epsilon_{wwi}(z,t))');
hold on
xline(3,':','DisplayName','Kurtosis = 3');
hold off
legend();
xlabel('Kurtosis');
ylabel('Depth [m]');
title('Kurt(\epsilon_{wwi}(z,t)): 2014-2020');

% savefig('figures/main/param/' + mooring + '_kurtosisOfEpsWwiAcrossDepth');
exportgraphics(ax0,'figures/main/param/' + mooring + '_kurtosisOfEpsWwiAcrossDepth.png');

%% eps_wwi: distribution of values

ax1 = figure;
histogram(sixDayEpsWwiM);
title('\epsilon_{wwi} (6dm): before binning');

% savefig('figures/main/param/' + mooring + '_epsWwi_distribution');
exportgraphics(ax1,'figures/main/param/' + mooring + '_epsWwi_distribution.png');

%% Hovmoeller Diagrams of eps_wwi

depthsIC3 = zmid_IC3(:,1);
time2 = time';
timeNum = datenum(time2);

[X2,Y2] = meshgrid(timeNum,depthsIC3);

%% eps_wwi: optimise the binning

% New optimisation: put all values < 1e-9 in one bin
% edges = [0 1e-9 1.3e-9 1.6e-9 1.9e-9 2.2e-9 2.5e-9 3.0e-9 3.5e-9 4.5e-9 5.5e-9 0.8e-8 1.7e-8 2.5e-8 14e-8];
edges = [0 1e-10 5e-10 1e-9 1.3e-9 1.6e-9 1.9e-9 3.0e-9 3.5e-9 4.5e-9 5.5e-9 0.8e-8 1.7e-8 2.5e-8 14e-8];

[epsWwiBinned] = discretize(sixDayEpsWwiM,edges);

% Distribution after binning
ax2 = figure;
histogram(epsWwiBinned);
title('\epsilon_{wwi} (6dm): after binning');
% savefig('figures/main/param/' + mooring + '_epsWwi_binnedDistribution');
exportgraphics(ax2,'figures/main/param/' + mooring + '_epsWwi_binnedDistribution.png');

%% eps_wwi: binned plot
NCbar = length(edges)-1;

ax3 = figure;
contourf(X2,Y2,epsWwiBinned,'LineColor','none');
datetick('x','yyyy mmm','keeplimits');
c = colorbar;
colormap(flipud(cbrewer2('Spectral',NCbar)));
c.Label.String = '\epsilon_{wwi} [W kg^{-1}]';
c = colorbar('YTick',1:15,'TickLabels',{num2str(edges(1)), num2str(edges(2)), num2str(edges(3)), num2str(edges(4)), ... 
    num2str(edges(5)), num2str(edges(6)), num2str(edges(7)), num2str(edges(8)), num2str(edges(9)), ...
    num2str(edges(10)), num2str(edges(11)), num2str(edges(12)), num2str(edges(13)), num2str(edges(14)), num2str(edges(15))});
xlabel('time');
ylabel('depth [m]');
zlabel('\epsilon_{wwi}');
title('binned epsilon');

% savefig('figures/main/param/' + mooring + '_epsWwi_Hovmoeller_Binned');
exportgraphics(ax3,'figures/main/param/' + mooring + '_epsWwi_Hovmoeller_Binned.png');

%% epsWwi IC3: not binned

ax3a = figure;
contourf(X2,Y2,log10(eps_wwi),-7.7:-0.1:-11,'LineColor','none');
datetick('x','yyyy mmm','keeplimits');
ylabel('Depth [m]');
colormap(flipud(cbrewer2('Spectral')));
c = colorbar;
ylim([-inf -100]);
c.Label.String = 'log_{10}(\epsilon_{wwi}) [log W kg^{-1}]';
% title('IC3: Turbulence Production due to WWI');
exportgraphics(ax3a,'figures/main/param/' + mooring + '_epsWwi.png');

%% eps_wwi: zoomed in plots of time series

% We want to see if we can see spring-neap cycles or other tidal features
% when we zoom in by a certain amount.

depthToPlotFrom = 17;

t1 = 21000;
t2 = 26000;
t3 = 26500;
t4 = 32500;

ax6 = figure;
sgtitle('\epsilon_{wwi} (6-day mean): close-up');
subplot(1,2,1)
contourf(X2(depthToPlotFrom:end,t1:t2),Y2(depthToPlotFrom:end,t1:t2),epsWwiBinned(depthToPlotFrom:end,t1:t2),'LineColor','none');
datetick('x','dd/mm/yy','keeplimits');
c = colorbar;
colormap(flipud(cbrewer2('Spectral',NCbar)));
c.Label.String = '\epsilon_{wwi} [W kg^{-1}]';
c = colorbar('YTick',1:15,'TickLabels',{num2str(edges(1)), num2str(edges(2)), num2str(edges(3)), num2str(edges(4)), ... 
    num2str(edges(5)), num2str(edges(6)), num2str(edges(7)), num2str(edges(8)), num2str(edges(9)), ...
    num2str(edges(10)), num2str(edges(11)), num2str(edges(12)), num2str(edges(13)), num2str(edges(14)), num2str(edges(15))});
xlabel('time');
ylabel('depth [m]');

subplot(1,2,2)
contourf(X2(depthToPlotFrom:end,t3:t4),Y2(depthToPlotFrom:end,t3:t4),epsWwiBinned(depthToPlotFrom:end,t3:t4),'LineColor','none');
datetick('x','dd/mm/yy','keeplimits');
c = colorbar;
colormap(flipud(cbrewer2('Spectral',NCbar)));
c.Label.String = '\epsilon_{wwi} [W kg^{-1}]';
c = colorbar('YTick',1:15,'TickLabels',{num2str(edges(1)), num2str(edges(2)), num2str(edges(3)), num2str(edges(4)), ... 
    num2str(edges(5)), num2str(edges(6)), num2str(edges(7)), num2str(edges(8)), num2str(edges(9)), ...
    num2str(edges(10)), num2str(edges(11)), num2str(edges(12)), num2str(edges(13)), num2str(edges(14)), num2str(edges(15))});
xlabel('time');
ylabel('depth [m]');

clear depthToPlotFrom t1 t2 t3 t4;

% savefig('figures/main/param/' + mooring + '_epsWwi-subplots');
exportgraphics(ax6,'figures/main/param/' + mooring + '_epsWwi-subplots.png');

%% Save parameters for the next file

save Matfiles/epsWwi.mat eps_wwi sixDayEpsWwiM kurtEpsWwi X2 Y2;