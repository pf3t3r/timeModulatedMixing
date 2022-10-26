clc; clear; close all;
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [5 5 40 22]);
set(0,'defaultAxesFontSize',18);

% Add subfolders to path
addpath(genpath(pwd));

load('Data/Mixing_parameterization_fields.mat');
load('Matfiles/SA_CT_interpolated.mat');
load('Matfiles/N_is.mat');
load('Matfiles/p_i.mat');
load('Matfiles/IC3.mat');
load('Matfiles/M1.mat');
mooring = ["IC3","M1"];

load('Matfiles/E_hil');
load('Matfiles/Fzt.mat');

%% Compute epsHil(z,t) for IC3 and M1

rhoIC3 = gsw_rho(SAmid_IC3_i,CTmid_IC3_i,pmid_IC3);
rhoM1 = gsw_rho(SAmid_M1_i,CTmid_M1_i,pmid_M1);
% timeNum = meshgrid(datenum(time'),zmid_IC3(:,1));
clear CT SA;

for i=1:length(time)
    eps_hil_IC3(:,i) = Ehil6dm_IC3(i).*verticalStructureIC3Hil(:,i)./rhoIC3(:,i);
    eps_hil_M1(:,i) = Ehil6dm_M1(i).*verticalStructureM1Hil(:,i)./rhoM1(:,i);
end
epsHM_IC3 = mean(eps_hil_IC3,'omitnan');
epsHM_M1 = mean(eps_hil_M1,'omitnan');

sixDayEpsHM_IC3 = movmean(eps_hil_IC3,149,2);
sixDayEpsHM_M1 = movmean(eps_hil_M1,149,2);

%% Kurtosis for epsHil at IC3 and M1

for i=1:length(zmid_IC3(:,1))
    kurtEpsHil_IC3(i) = kurtosis(sixDayEpsHM_IC3(i,:));
    kurtEpsHillog10_IC3(i) = kurtosis(log10(sixDayEpsHM_IC3(i,:)));
%     kurtEpsHil_M1(i) = kurtosis(sixDayEpsHM_M1(i,:));
%     kurtEpsHillog10_M1(i) = kurtosis(log10(sixDayEpsHM_M1(i,:)));
end

for i = 1:length(zmid_M1(:,1))
    kurtEpsHil_M1(i) = kurtosis(sixDayEpsHM_M1(i,:));
    kurtEpsHillog10_M1(i) = kurtosis(log10(sixDayEpsHM_M1(i,:)));
end

ax1 = figure;
subplot(1,2,1)
plot(kurtEpsHil_IC3,zmid_IC3(:,1),'DisplayName','Kurt(\epsilon_{hil}(z,t))');
hold on
plot(kurtEpsHillog10_IC3,zmid_IC3(:,i),'DisplayName','Kurt(log(\epsilon_{hil}(z,t)))');
xline(3,':','DisplayName','Kurtosis = 3');
hold off
legend('Location','best');
xlabel('Kurtosis(\epsilon_{hil}(z,t))');
ylabel('Depth [m]');
title('IC3 (2014-2020)');

subplot(1,2,2)
plot(kurtEpsHil_M1,zmid_M1(:,1),'DisplayName','Kurt(\epsilon_{hil}(z,t))');
hold on
plot(kurtEpsHillog10_M1,zmid_M1(:,i),'DisplayName','Kurt(log(\epsilon_{hil}(z,t)))');
xline(3,':','DisplayName','Kurtosis = 3');
hold off
legend('Location','best');
xlabel('Kurtosis(\epsilon_{hil}(z,t))');
ylabel('Depth [m]');
title('M1 (2014-2020)');

savefig('figures/main/param/_kurtosisOfEpsHilAcrossDepth');
exportgraphics(ax1,'figures/main/param/_kurtosisOfEpsHilAcrossDepth.png');

%% IC3: epsHil distribution

ax2 = figure;
histogram(sixDayEpsHM_IC3);
title('\epsilon_{hil} (6dm) (before binning)');

savefig('figures/main/param/' + mooring(1) + '_epsHil_distribution');
exportgraphics(ax2,'figures/main/param/' + mooring(1) + '_epsHil_distribution.png');

%% M1: epsHil distribution

ax3 = figure;
histogram(sixDayEpsHM_M1);
title('M1: \epsilon_{hil} (6dm) (before binning)');

savefig('figures/main/param/' + mooring(2) + '_epsHil_distribution');
exportgraphics(ax3,'figures/main/param/' + mooring(2) + '_epsHil_distribution.png');

%% epsHil: initialise meshgrid for contour plots

[X2,Y2] = meshgrid(datenum(time'),zmid_IC3(:,1));
[X2a,Y2a] = meshgrid(datenum(time'),zmid_M1(:,1));

%% optimise epsHil binning for IC3 and M1

% Original optimisation - tuned according to histogram
% edges = [0 3.3e-10 5.5e-10 1e-9 1.4e-9 1.8e-9 2.1e-9 2.5e-9 3.5e-9 4.5e-9 5.5e-9 0.8e-8 1.7e-8 2.5e-8 14e-8];

% New optimisation: put all values < 1e-9 in one bin
edges = [0 1e-9 1.3e-9 1.6e-9 1.9e-9 2.2e-9 2.5e-9 3.0e-9 3.5e-9 4.5e-9 5.5e-9 0.8e-8 1.7e-8 2.5e-8 14e-8];

epsHilBinned_IC3 = discretize(sixDayEpsHM_IC3,edges);
epsHilBinned_M1 = discretize(sixDayEpsHM_M1,edges);

% Distribution after binning
ax4 = figure;
histogram(epsHilBinned_IC3);
title('IC3: \epsilon_{hil} (6dm) after binning');
savefig('figures/main/param/' + mooring(1) + '_epsHil_binnedDistribution');
exportgraphics(ax4,'figures/main/param/' + mooring(1) + '_epsHil_binnedDistribution.png');


ax5 = figure;
histogram(epsHilBinned_M1);
title('M1: \epsilon_{hil} (6dm) after binning');
savefig('figures/main/param/' + mooring(2) + '_epsHil_binnedDistribution');
exportgraphics(ax5,'figures/main/param/' + mooring(2) + '_epsHil_binnedDistribution.png');

%% IC3: epsHil (binned)
NCbar = length(edges)-1;

ax6 = figure;
contourf(X2,Y2,epsHilBinned_IC3,'LineColor','none');
datetick('x','yyyy mmm','keeplimits');
colormap(flipud(cbrewer2('Spectral',NCbar)));
c = colorbar;
c.Ticks = 1:(1-1/NCbar):length(edges);
c.TickLabels = {num2str(edges(1)), num2str(edges(2)), num2str(edges(3)), num2str(edges(4)), ... 
    num2str(edges(5)), num2str(edges(6)), num2str(edges(7)), num2str(edges(8)), num2str(edges(9)), ...
    num2str(edges(10)), num2str(edges(11)), num2str(edges(12)), num2str(edges(13)), num2str(edges(14)), num2str(edges(15))};
c.Label.String = '\epsilon_{hil} [W kg^{-1}]';
xlabel('time');
ylabel('depth [m]');
zlabel('\epsilon_{hil}');
% title('binned epsilon');

savefig('figures/main/param/' + mooring(1) + '_epsHil_Hovmoeller_Binned');
exportgraphics(ax6,'figures/main/param/' + mooring(1) + '_epsHil_Hovmoeller_Binned.png');

%% IC3: epsHil (logplot)

ax7 = figure;
contourf(X2,Y2,log10(eps_hil_IC3),-10:0.1:-7,'LineColor','none');
datetick('x','yyyy mmm','keeplimits');
% cmocean('tempo');
colormap(flipud(cbrewer2('Spectral')));
c = colorbar;
c.Label.String = 'log_{10}(\epsilon_{hil}) [log W kg^{-1}]';
ylabel('Depth [m]');
ylim([-inf -100]);
exportgraphics(ax7,'figures/main/param/' + mooring(1) + '_epsHil.png');

%% M1: epsHil (binned)

ax8 = figure;
contourf(X2a,Y2a,epsHilBinned_M1,'LineColor','none');
datetick('x','yyyy mmm','keeplimits');
colormap(flipud(cbrewer2('Spectral',NCbar)));
c = colorbar;
c.Ticks = 1:(1-1/NCbar):length(edges);
c.TickLabels = {num2str(edges(1)), num2str(edges(2)), num2str(edges(3)), num2str(edges(4)), ... 
    num2str(edges(5)), num2str(edges(6)), num2str(edges(7)), num2str(edges(8)), num2str(edges(9)), ...
    num2str(edges(10)), num2str(edges(11)), num2str(edges(12)), num2str(edges(13)), num2str(edges(14)), num2str(edges(15))};
c.Label.String = '\epsilon_{hil} [W kg^{-1}]';
xlabel('time');
ylabel('depth [m]');
zlabel('\epsilon_{hil}');
% title('M1: binned epsilon');

savefig('figures/main/param/' + mooring(2) + '_epsHil_Hovmoeller_Binned');
exportgraphics(ax8,'figures/main/param/' + mooring(2) + '_epsHil_Hovmoeller_Binned.png');

%% M1 epsHil (logplot)

ax9 = figure;
contourf(X2a,Y2a,log10(eps_hil_M1),-10:0.1:-7,'LineColor','none');
datetick('x','yyyy mmm','keeplimits');
% cmocean('tempo');
colormap(flipud(cbrewer2('Spectral')));
c = colorbar;
c.Label.String = 'log_{10}(\epsilon_{hil}) [log W kg^{-1}]';
ylabel('Depth [m]');
ylim([-inf -60]);
exportgraphics(ax9,'figures/main/param/' + mooring(2) + '_epsHil.png');

%% Same again: but this time with imagesc

% ax5 = figure;
% imagesc(epsHilBinned_IC3);
% xlabel('time [# of hrs]');
% ylabel('depth [levels]');
% title('\epsilon_{hil}: imagesc visualisation');
% 
% savefig('figures/main/param/' + mooring(1) + '_epsHil_Hovmoeller_imagesc');
% exportgraphics(ax5,'figures/main/param/' + mooring(1) + '_epsHil_Hovmoeller_imagesc.png');

%% IC3: Bottom over six years

ax10 = figure;
contourf(X2(36:end,:),Y2(36:end,:),epsHilBinned_IC3(36:end,:),'LineColor','none');
datetick('x','yyyy mmm','keeplimits');

colormap(flipud(cbrewer2('Spectral',NCbar)));
c = colorbar;
c.Ticks = 1:(1-1/NCbar):length(edges);
c.TickLabels = {num2str(edges(1)), num2str(edges(2)), num2str(edges(3)), num2str(edges(4)), ... 
    num2str(edges(5)), num2str(edges(6)), num2str(edges(7)), num2str(edges(8)), num2str(edges(9)), ...
    num2str(edges(10)), num2str(edges(11)), num2str(edges(12)), num2str(edges(13)), num2str(edges(14)), num2str(edges(15))};
c.Label.String = '\epsilon_{hil} [W kg^{-1}]';
xlabel('time');
ylabel('depth [m]');
zlabel('\epsilon_{hil}');
title('binned epsilon');

savefig('figures/main/param/' + mooring(1) + '_epsHil_Hovmoeller_Binned_zoom');
exportgraphics(ax10,'figures/main/param/' + mooring(1) + '_epsHil_Hovmoeller_Binned_zoom.png');


%% M1: Bottom over six years

ax11 = figure;
contourf(X2a(36:end,:),Y2a(36:end,:),epsHilBinned_M1(36:end,:),'LineColor','none');
datetick('x','yyyy mmm','keeplimits');
colormap(flipud(cbrewer2('Spectral',NCbar)));
c = colorbar;
c.Ticks = 1:(1-1/NCbar):length(edges);
c.TickLabels = {num2str(edges(1)), num2str(edges(2)), num2str(edges(3)), num2str(edges(4)), ... 
    num2str(edges(5)), num2str(edges(6)), num2str(edges(7)), num2str(edges(8)), num2str(edges(9)), ...
    num2str(edges(10)), num2str(edges(11)), num2str(edges(12)), num2str(edges(13)), num2str(edges(14)), num2str(edges(15))};
c.Label.String = '\epsilon_{hil} [W kg^{-1}]';
xlabel('time');
ylabel('depth [m]');
zlabel('\epsilon_{hil}');
title('M1: binned epsilon');

savefig('figures/main/param/' + mooring(2) + '_epsHil_Hovmoeller_Binned_zoom');
exportgraphics(ax11,'figures/main/param/' + mooring(2) + '_epsHil_Hovmoeller_Binned_zoom.png');


%% Parameters for zoom
depthToPlotFrom = 33;

t1 = 5000;
t2 = 10000;
t3 = 30000;
t4 = 36000;

%% IC3: Close-up of bottom in 2015 and 2018

% We want to see if we can see spring-neap cycles or other tidal features
% when we zoom in by a certain amount.

disp('IC3');
disp(NCbar);

ax12 = figure;
sgtitle('IC3: \epsilon_{hil} (6-day mean): close-up');
subplot(1,2,1)
contourf(X2(depthToPlotFrom:end,t1:t2),Y2(depthToPlotFrom:end,t1:t2),epsHilBinned_IC3(depthToPlotFrom:end,t1:t2),'LineColor','none');
datetick('x','dd/mm/yy','keeplimits');
colormap(flipud(cbrewer2('Spectral',NCbar)));
c = colorbar;
c.Ticks = 1:(1-1/NCbar):length(edges);
c.TickLabels = {num2str(edges(1)), num2str(edges(2)), num2str(edges(3)), num2str(edges(4)), ... 
    num2str(edges(5)), num2str(edges(6)), num2str(edges(7)), num2str(edges(8)), num2str(edges(9)), ...
    num2str(edges(10)), num2str(edges(11)), num2str(edges(12)), num2str(edges(13)), num2str(edges(14)), num2str(edges(15))};
c.Label.String = '\epsilon_{hil} [W kg^{-1}]';
xlabel('time');
ylabel('depth [m]');

subplot(1,2,2)
contourf(X2(depthToPlotFrom:end,t3:t4),Y2(depthToPlotFrom:end,t3:t4),epsHilBinned_IC3(depthToPlotFrom:end,t3:t4),'LineColor','none');
datetick('x','dd/mm/yy','keeplimits');
c = colorbar;
c.Ticks = 1:(1-1/NCbar):length(edges);
c.TickLabels = {num2str(edges(1)), num2str(edges(2)), num2str(edges(3)), num2str(edges(4)), ... 
    num2str(edges(5)), num2str(edges(6)), num2str(edges(7)), num2str(edges(8)), num2str(edges(9)), ...
    num2str(edges(10)), num2str(edges(11)), num2str(edges(12)), num2str(edges(13)), num2str(edges(14)), num2str(edges(15))};
c.Label.String = '\epsilon_{hil} [W kg^{-1}]';
xlabel('time');
ylabel('depth [m]');

savefig('figures/main/param/' + mooring(1) + '_epsHil-subplots');
exportgraphics(ax12,'figures/main/param/' + mooring(1) + '_epsHil-subplots.png');

%% M1: close-up of bottom in 2015 and 2018

disp('M1');
disp(NCbar);

ax13 = figure;
sgtitle('M1: \epsilon_{hil} (6-day mean): close-up');
subplot(1,2,1)
contourf(X2a(depthToPlotFrom:end,t1:t2),Y2a(depthToPlotFrom:end,t1:t2),epsHilBinned_M1(depthToPlotFrom:end,t1:t2),'LineColor','none');
datetick('x','dd/mm/yy','keeplimits');
colormap(flipud(cbrewer2('Spectral',NCbar)));
c = colorbar;
c.Ticks = 1:(1-1/NCbar):length(edges);
c.TickLabels = {num2str(edges(1)), num2str(edges(2)), num2str(edges(3)), num2str(edges(4)), ... 
    num2str(edges(5)), num2str(edges(6)), num2str(edges(7)), num2str(edges(8)), num2str(edges(9)), ...
    num2str(edges(10)), num2str(edges(11)), num2str(edges(12)), num2str(edges(13)), num2str(edges(14)), num2str(edges(15))};
c.Label.String = '\epsilon_{hil} [W kg^{-1}]';
xlabel('time');
ylabel('depth [m]');

subplot(1,2,2)
contourf(X2a(depthToPlotFrom:end,t3:t4),Y2a(depthToPlotFrom:end,t3:t4),epsHilBinned_M1(depthToPlotFrom:end,t3:t4),'LineColor','none');
datetick('x','dd/mm/yy','keeplimits');
colormap(flipud(cbrewer2('Spectral',NCbar)));
c = colorbar;
c.Ticks = 1:(1-1/NCbar):length(edges);
c.TickLabels = {num2str(edges(1)), num2str(edges(2)), num2str(edges(3)), num2str(edges(4)), ... 
    num2str(edges(5)), num2str(edges(6)), num2str(edges(7)), num2str(edges(8)), num2str(edges(9)), ...
    num2str(edges(10)), num2str(edges(11)), num2str(edges(12)), num2str(edges(13)), num2str(edges(14)), num2str(edges(15))};
c.Label.String = '\epsilon_{hil} [W kg^{-1}]';
xlabel('time');
ylabel('depth [m]');

savefig('figures/main/param/' + mooring(1) + '_epsHil-subplots');
exportgraphics(ax13,'figures/main/param/' + mooring(1) + '_epsHil-subplots.png');

clear depthToPlotFrom t1 t2 t3 t4;

%% Save parameters for the next file

save Matfiles/epsHil.mat eps_hil_IC3 sixDayEpsHM_IC3 kurtEpsHil_IC3 X2 Y2 ...
    eps_hil_M1 sixDayEpsHM_M1 kurtEpsHil_M1 X2a Y2a ...
    rhoIC3 rhoM1;
