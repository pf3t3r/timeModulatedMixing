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

% ax1 = figure;
% histogram(sixDayEpsWwiM);
% title('\epsilon_{wwi} (6dm): before binning');
% 
% % savefig('figures/main/param/' + mooring + '_epsWwi_distribution');
% exportgraphics(ax1,'figures/main/param/' + mooring + '_epsWwi_distribution.png');

%% Hovmoeller Diagrams of eps_wwi

depthsIC3 = zmid_IC3(:,1);
time2 = time';
timeNum = datenum(time2);

[X2,Y2] = meshgrid(timeNum,depthsIC3);

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

%% Save parameters for the next file

save Matfiles/epsWwi.mat eps_wwi sixDayEpsWwiM kurtEpsWwi X2 Y2;