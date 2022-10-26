clc; clear; close all;
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [5 2 40 22]);
set(0,'defaultAxesFontSize',18);

% Add subfolders to path
addpath(genpath(pwd));

load('Matfiles/E_hil');
load('Data/Mixing_parameterization_fields.mat');

load('Matfiles/N_is.mat');
load("Matfiles/M1.mat",'time');

depths_IC3 = load('Matfiles/IC3.mat').depths;
depths_M1 = load('Matfiles/M1.mat').depths_M1;
% load N2v2.mat;
% load IC3.mat;
mooring = ["IC3","M1"];

%% STUFF to do with buoyancy
%% Initialise dz

% dz_IC3 = cat(1,10*ones(5,1),15*ones(6,1),50*ones(30,1)); % THIS IS WRONG,
% SEE CORRECTED VERSION BELOW!!!
dz_IC3 = cat(1,10*ones(5,1),25*ones(6,1),50*ones(30,1));
dz_M1 = cat(1,10*ones(5,1),25*ones(6,1),50*ones(31,1));

% not used here ... used in N^2 -> put this in the right place!!
timeNum_IC3 = meshgrid(datenum(time'),zmid_IC3(:,1));
timeNum_M1 = meshgrid(datenum(time'),zmid_M1(:,1));
% depthLevels_IC3 = meshgrid(datenum(time'),depths_IC3);

%% Integrate N^2
for i=1:length(time)
    N2_int_overTime_IC3(i) = sum(N2_IC3(:,i).*dz_IC3,'omitnan');
    N2_int_overTime_M1(i) = sum(N2_M1(:,i).*dz_M1,'omitnan');
end

N2_int_3D_overTime_IC3 = repmat(N2_int_overTime_IC3,1,1,length(zmid_IC3(:,1)));
N2_int_3D_overTime_IC3 = squeeze(N2_int_3D_overTime_IC3(1,:,:))';

N2_int_3D_overTime_M1 = repmat(N2_int_overTime_M1,1,1,length(zmid_M1(:,1)));
N2_int_3D_overTime_M1 = squeeze(N2_int_3D_overTime_M1(1,:,:))';

%% IC3: N^2(t) [vertically-integrated N^2(t)]

ax1 = figure;
plot(time,N2_int_overTime_IC3);
xlabel('time');
ylabel('\int N^2 dz (t)');
title('\int N^2 dz (t)');

% savefig('figures/main/param/' + mooring(1) + '_N2_integrated_over_time');
exportgraphics(ax1,'figures/main/param/' + mooring(1) + '_N2_integrated_over_time.png');

%% M1: N^2(t) [vertically-integrated N^2(t)]

ax2 = figure;
plot(time,N2_int_overTime_M1);
xlabel('time');
ylabel('\int N^2 dz (t)');
title('\int N^2 dz (t)');

% savefig('figures/main/param/' + mooring(2) + '_N2_integrated_over_time');
exportgraphics(ax2,'figures/main/param/' + mooring(2) + '_N2_integrated_over_time.png');

%% IC3: N^2(z,t): histogram (before binning)

figure
% histogram(N2_IC3);
histogram(N_IC3);
title('IC3: N(z,t)');

%% M1: N^2(z,t): histogram (before binning)

figure
% histogram(N2_M1);
histogram(N_M1);
title('M1: N(z,t)');

%% IC3: N^2(z,t): binned contour plot

% edges = [0 1e-7 2e-7:2e-7:2e-6 3e-6:3e-6:3e-5];
% edges = linspace(0,3e-5,15);
edges = linspace(5e-5,6e-3,15);
binnedN_IC3 = discretize(N_IC3,edges);
NCbar = length(edges) - 1;

ax3 = figure;
contourf(timeNum_IC3,zmid_IC3,binnedN_IC3,'LineColor','none');
datetick('x','yyyy mmm','keeplimits');
cmap = colormap(flipud(cbrewer2('Spectral',NCbar)));
c = colorbar;
c.Ticks = 1:(1-1/NCbar):length(edges);
c.TickLabels = {num2str(edges(1)), num2str(edges(2)), num2str(edges(3)), num2str(edges(4)), ... 
    num2str(edges(5)), num2str(edges(6)), num2str(edges(7)), num2str(edges(8)), num2str(edges(9)), ...
    num2str(edges(10)), num2str(edges(11)), num2str(edges(12)), num2str(edges(13)), num2str(edges(14)), ...
    num2str(edges(15))};
% , num2str(edges(16)), num2str(edges(17)), num2str(edges(18)), num2str(edges(19)), ...
%     num2str(edges(20)), num2str(edges(21)), num2str(edges(22))};
% c = colorbar('YTick',1:length(edges), ...
%     'TickLabels',{num2str(edges(1)), num2str(edges(2)), num2str(edges(3)), num2str(edges(4)), ... 
%     num2str(edges(5)), num2str(edges(6)), num2str(edges(7)), num2str(edges(8)), num2str(edges(9)), ...
%     num2str(edges(10)), num2str(edges(11)), num2str(edges(12)), num2str(edges(13)), num2str(edges(14)), ...
%     num2str(edges(15)), num2str(edges(16)), num2str(edges(17)), num2str(edges(18)), num2str(edges(19)), ...
%     num2str(edges(20)), num2str(edges(21)), num2str(edges(22))});
c.Label.String = 'N^2(z,t) [s^{-2}]';
xlabel('time');
ylabel('depth [m]');
% title('IC3: N^2(z,t)');

% savefig('figures/main/param/' + mooring(1) + '_N2_zt_binned');
exportgraphics(ax3,'figures/main/param/' + mooring(1) + '_N2_zt_binned.png');

%% IC3: N(z,t) log plot

ax4 = figure;
contourf(timeNum_IC3,zmid_IC3,log10(N_IC3),-3.7:0.05:-2.4,'LineColor','none');
ylabel('Depth [m]');
datetick('x','yyyy mmm','keeplimits');
colormap(flipud(cbrewer2('Spectral')));
ylim([-inf -100]);
% cmocean('tempo');
c = colorbar;
c.Label.String = 'log_{10}(N(z,t)) [s^{-1}]';
exportgraphics(ax4,'figures/main/param/' + mooring(1) + '_N.png');

%% M1: N^2(z,t): binned contour plot

% edges = [0 1e-7 2e-7:2e-7:2e-6 3e-6:3e-6:3e-5];
edges = linspace(5e-5,1.4e-2,15);
% edges = [5e-5:0.8e-4:1e-3 5e-3 1e-2 1.5e-2];
binnedN_M1 = discretize(N_M1,edges);
NCbar = length(edges) - 1;

ax5 = figure;
contourf(timeNum_M1,zmid_M1,binnedN_M1,'LineColor','none');
datetick('x','yyyy mmm','keeplimits');
cmap = colormap(flipud(cbrewer2('Spectral',NCbar)));
c = colorbar;
c.Ticks = 1:(1-1/NCbar):length(edges);
c.TickLabels = {num2str(edges(1)), num2str(edges(2)), num2str(edges(3)), num2str(edges(4)), ... 
    num2str(edges(5)), num2str(edges(6)), num2str(edges(7)), num2str(edges(8)), num2str(edges(9)), ...
    num2str(edges(10)), num2str(edges(11)), num2str(edges(12)), num2str(edges(13)), num2str(edges(14)), ...
    num2str(edges(15))};
% , num2str(edges(16)), num2str(edges(17)), num2str(edges(18)), num2str(edges(19)), ...
%     num2str(edges(20)), num2str(edges(21)), num2str(edges(22))};
c.Label.String = 'N^2(z,t) [s^{-2}]';
xlabel('time');
ylabel('depth [m]');
% title('M1: N^2(z,t)');

% savefig('figures/main/param/' + mooring(2) + '_N2_zt_binned');
exportgraphics(ax5,'figures/main/param/' + mooring(2) + '_N2_zt_binned.png');

%% M1: N(z,t) log plot
ax6 = figure;
contourf(timeNum_M1,zmid_M1,log10(N_M1),-4.1:0.05:-2,'LineColor','none');
datetick('x','yyyy mmm','keeplimits');
% cmocean('tempo');
colormap(flipud(cbrewer2('Spectral')));
c = colorbar;
ylim([-inf -60]);
c.Label.String = 'log_{10}(N(z,t)) [s^{-1}]';
exportgraphics(ax6,'figures/main/param/' + mooring(2) + '_N.png');

%% IC3: N^2(z,t): histogram (after binning)

figure
histogram(binnedN_IC3);
title('IC3 after binning');

%% M1: N^2(z,t): histogram (after binning)

figure
histogram(binnedN_M1);
title('M1 after binning');

%% original Fzt file after this point

%% Modulate Hbot

% Hbot for IC3 and M1
Hbot_IC3 = decay_scale_bot(657,279);
Hbot_M1 = decay_scale_bot(660,279);

% Normalisation Factor 'nf'
nf = zeros(2,length(time));
for t = 1:length(time)
    nf(1,t) = 1/sqrt(Ehil6dm_IC3(t));
    nf(2,t) = 1/sqrt(Ehil6dm_M1(t));
end

% Tuning Factor 'tf'
tf = zeros(2,length(time));

for i = 1:2
    tf(i,:) = 1./mean((nf(i,:) - min(nf(i,:))) ./ (max(nf(i,:)) - min(nf(i,:))));

    % Normalise such that the mean of nf = E_hil(x,y)
    normalisation(i,:) = tf(i,:) .* (nf(i,:) - min(nf(i,:))) ./ (max(nf(i,:)) - min(nf(i,:)));
end

HbotIC3_norm = normalisation(1,:).*Hbot_IC3;
HbotM1_norm = normalisation(2,:).*Hbot_M1;
meanHBotIC3preNorm = mean(HbotIC3_norm);
meanHBotM1preNorm = mean(HbotM1_norm);

% Apply a floor to the modulated values for Hbot
for i=1:length(time)
    if HbotIC3_norm(i) < 100
        HbotIC3_norm(i) = 100;
    end
    if HbotM1_norm(i) < 100
        HbotM1_norm(i) = 100;
    end
end

ax6a = figure;
plot(time,HbotIC3_norm);
hold on
yline(meanHBotIC3preNorm,'LineWidth',2,'Color','blue');
yline(Hbot_IC3);
hold off
title('H_{bot}(t): IC3');
exportgraphics(ax6a,'figures/main/param/' + mooring(1) + '_HbotMod.png');

ax6b = figure;
plot(time,HbotM1_norm);
hold on
yline(meanHBotM1preNorm,'LineWidth',2,'Color','blue');
yline(Hbot_M1);
hold off
title('H_{bot}(t): M1');
exportgraphics(ax6b,'figures/main/param/' + mooring(2) + '_HbotMod.png');

%% IC3: Compute the vertical structure

% depthsIC3 = zq;
hab_IC3 = [];

% From Deployment Overview [De Jong & Fried]
IC3bathy = 1635;
for i = 1:length(zmid_IC3(:,1))
    hab_IC3(i) = IC3bathy + zmid_IC3(i,1);
end
hab_IC3(end) = 0;

if isfile('Matfiles/verticalStructureIC3.mat')
    disp('Time-dependent vertical structure already calculated');
    load('Matfiles/verticalStructureIC3.mat');
else
    disp('Calculating vertical structure...');
    r_bot = 0.86;
    verticalStructureIC3Hil = NaN(length(zmid_IC3(:,1)),length(time));
    for t = 1:length(time)
        for k=1:length(zmid_IC3(:,1))
%             T1(k,t) =  r_bot .* (1./(1 + (hab_IC3(k)./decay_scale_bot(657,279))).^2) ...
%                             .*(1./bathy(657,279) + 1./decay_scale_bot(657,279));
            T1(k,t) =  r_bot .* (1./(1 + (hab_IC3(k)./HbotIC3_norm(t))).^2) ...
                            .*(1./IC3bathy + 1./HbotIC3_norm(t));
            T2(k,t) = (1 - r_bot).*(N2_IC3(k,t)./N2_int_3D_overTime_IC3(k,t));
        end
    disp(t)
    end
    verticalStructureIC3Hil = T1 + T2;
    save Matfiles/verticalStructureIC3.mat verticalStructureIC3Hil T1 T2;
end


%% M1: Vertical Structure

% M1 coords = [660,279]

% M1 Bathymetry - from Deployment Overview [De Jong & Fried]
M1bathy = 1712;
for i = 1:length(zmid_M1(:,1))
    hab_M1(i) = M1bathy + zmid_M1(i,1);
end
hab_M1(end) = 0;

if isfile('Matfiles/verticalStructureM1.mat')
    disp('Time-dependent vertical structure for M1 already calculated');
    load('Matfiles/verticalStructureM1.mat');
else
    disp('Calculating vertical structure for M1...');
    r_bot = 0.86;
    %verticalStructureIC3Hil = NaN(length(depthsIC3),length(time));
    for t = 1:length(time)
        for k=1:length(zmid_M1(:,1))
%             T1(k,t) =  r_bot .* (1./(1 + (hab_M1(k)./decay_scale_bot(660,279))).^2) ...
%                             .*(1./bathy(660,279) + 1./decay_scale_bot(660,279));
            T1(k,t) =  r_bot .* (1./(1 + (hab_M1(k)./HbotM1_norm(t))).^2) ...
                            .*(1./M1bathy + 1./HbotM1_norm(t));
            T2(k,t) = (1 - r_bot).*(N2_M1(k,t)./N2_int_3D_overTime_M1(k,t));
        end
    disp(t)
    end
    verticalStructureM1Hil = T1 + T2;
    save Matfiles/verticalStructureM1.mat verticalStructureM1Hil T1 T2;
end
%% Set up 2D-grid for contourf

% time2 = time';
% timeNum = datenum(time2);
% [timeNum2,depths2] = meshgrid(timeNum,depths_IC3);

% load verticalStructure.mat;

%% IC3: Histogram: F(z,t)

% Use this as a justification for binning scheme.

ax7 = figure;
histogram(verticalStructureIC3Hil);
title('F(z,t): before binning');

savefig('figures/main/param/' + mooring(1) + '_hist_verticalStructure');
exportgraphics(ax7,'figures/main/param/' + mooring(1) + '_hist_verticalStructure.png');

%% M1: Histogram for F(z,t)

ax8 = figure;
histogram(verticalStructureIC3Hil);
title('F(z,t): before binning');

% savefig('figures/main/param/' + mooring(2) + '_hist_verticalStructure');
exportgraphics(ax8,'figures/main/param/' + mooring(2) + '_hist_verticalStructure.png');


%% IC3: F(z,t), binned

timeEnd = length(time);

edges = [0:0.9e-4:0.9e-3 1e-3:2e-3:7e-3];
NCbar = length(edges) - 1;
binnedEpsilon = discretize(verticalStructureIC3Hil,edges);

ax9 = figure;
contourf(timeNum_IC3(:,1:timeEnd)',zmid_IC3(:,1:timeEnd)',binnedEpsilon(:,1:timeEnd)','LineColor','none');
datetick('x','yy mmm','keeplimits','keepticks');
cmap = colormap(flipud(cbrewer2('Spectral',NCbar)));
c = colorbar;
c.Ticks = 1:(1-1/NCbar):length(edges);
c.TickLabels = {num2str(edges(1)), num2str(edges(2)), num2str(edges(3)), num2str(edges(4)), ... 
    num2str(edges(5)), num2str(edges(6)), num2str(edges(7)), num2str(edges(8)), num2str(edges(9)), ...
    num2str(edges(10)), num2str(edges(11)), num2str(edges(12)), num2str(edges(13)), num2str(edges(14)), num2str(edges(15))};
c.Label.String = 'F(z,t) [m^{-1}]';
xlabel('time');
ylabel('depth [m]');
% title('IC3: F(z,t), binned');
% savefig('figures/main/param/' + mooring(1) + '_verticalStructureFunctionNEW');
exportgraphics(ax9,'figures/main/param/' + mooring(1) + '_verticalStructureFunctionNEW.png');

%% IC3: F(z,t), not binned

ax9a = figure;
contourf(timeNum_IC3,zmid_IC3,log10(verticalStructureIC3Hil),-2.5:-0.1:-4.5,'LineColor','none');
datetick('x','yyyy mmm','keeplimits');
% cmocean('tempo');
colormap(flipud(cbrewer2('Spectral')));
c = colorbar;
c.Label.String = 'log_{10}(F(z,t)) [m^{-1}]';
ylim([-inf -100]);
ylabel('Depth [m]');
title('Vertical Structure');
exportgraphics(ax9a,'figures/main/param/' + mooring(1) + '_verticalStructure.png');

%% M1: F(z,t), binned

edges = [2e-4:0.8e-4:1.3e-3 3e-2];
NCbar = length(edges) - 1;
binnedEpsilonM1 = discretize(verticalStructureM1Hil,edges);

ax10 = figure;
contourf(timeNum_M1(:,1:timeEnd)',zmid_M1(:,1:timeEnd)',binnedEpsilonM1(:,1:timeEnd)','LineColor','none');
datetick('x','yyyy mmm','keeplimits');
cmap = colormap(flipud(cbrewer2('Spectral',NCbar)));
c = colorbar;
c.Ticks = 1:(1-1/NCbar):length(edges);
c.TickLabels = {num2str(edges(1)), num2str(edges(2)), num2str(edges(3)), num2str(edges(4)), ... 
    num2str(edges(5)), num2str(edges(6)), num2str(edges(7)), num2str(edges(8)), num2str(edges(9)), ...
    num2str(edges(10)), num2str(edges(11)), num2str(edges(12)), num2str(edges(13)), num2str(edges(14)), num2str(edges(15))};
% cmap = flipud(cbrewer2('div','Spectral',NCbar));
% colormap(cmap);
% c = colorbar('YTick',1:14,'TickLabels',{num2str(edges(1)), num2str(edges(2)), num2str(edges(3)), num2str(edges(4)), ... 
%     num2str(edges(5)), num2str(edges(6)), num2str(edges(7)), num2str(edges(8)), num2str(edges(9)), ...
%     num2str(edges(10)), num2str(edges(11)), num2str(edges(12)), num2str(edges(13)), num2str(edges(14))});
c.Label.String = 'F(z,t) [m^{-1}]';
xlabel('time');
ylabel('depth [m]');
% title('M1: F(z,t), binned');
% savefig('figures/main/param/' + mooring(2) + '_verticalStructure');
exportgraphics(ax10,'figures/main/param/' + mooring(2) + '_verticalStructure.png');

%% M1: F(z,t), not binned

ax10a = figure;
contourf(timeNum_M1,zmid_M1,log10(verticalStructureM1Hil),-2.5:-0.1:-4.5,'LineColor','none');
datetick('x','yyyy mmm','keeplimits');
% cmocean('tempo');
colormap(flipud(cbrewer2('Spectral')));
c = colorbar;
c.Label.String = 'log_{10}(F(z,t)) [m^{-1}]';
ylabel('Depth [m]');
ylim([-inf -60]);
title('Vertical Structure');
exportgraphics(ax10a,'figures/main/param/' + mooring(2) + '_verticalStructure.png');

%% Zoom parameters

% depthToPlotFrom = 33;
% t1 = 9000;
% t2 = 10000;
% t3 = 30000;
% t4 = 31000;

%% IC3: Zoomed-in F(z,t), binned

% ax11 = figure;
% % sgtitle('F(z,t): zoomed-in');
% subplot(1,2,1)
% contourf(timeNum_IC3(depthToPlotFrom:end,t1:t2)',zmid_IC3(depthToPlotFrom:end,t1:t2)',binnedEpsilon(depthToPlotFrom:end,t1:t2)','LineColor','none');
% datetick('x','yyyy mmm','keeplimits');
% cmap = flipud(cbrewer2('div','Spectral',NCbar));
% colormap(cmap);
% c = colorbar('YTick',1:14,'TickLabels',{num2str(edges(1)), num2str(edges(2)), num2str(edges(3)), num2str(edges(4)), ... 
%     num2str(edges(5)), num2str(edges(6)), num2str(edges(7)), num2str(edges(8)), num2str(edges(9)), ...
%     num2str(edges(10)), num2str(edges(11)), num2str(edges(12)), num2str(edges(13)), num2str(edges(14))});
% xlabel('time');
% ylabel('depth [m]');
% 
% subplot(1,2,2)
% contourf(timeNum_IC3(depthToPlotFrom:end,t3:t4)',zmid_IC3(depthToPlotFrom:end,t3:t4)',binnedEpsilon(depthToPlotFrom:end,t3:t4)','LineColor','none');
% datetick('x','yyyy mmm','keeplimits');
% cmap = flipud(cbrewer2('div','Spectral',NCbar));
% colormap(cmap);
% c = colorbar('YTick',1:14,'TickLabels',{num2str(edges(1)), num2str(edges(2)), num2str(edges(3)), num2str(edges(4)), ... 
%     num2str(edges(5)), num2str(edges(6)), num2str(edges(7)), num2str(edges(8)), num2str(edges(9)), ...
%     num2str(edges(10)), num2str(edges(11)), num2str(edges(12)), num2str(edges(13)), num2str(edges(14))});
% xlabel('time');
% ylabel('depth [m]');
% 
% % savefig('figures/main/param/' + mooring(1) + '_verticalStructureFunctionZoom');
% exportgraphics(ax11,'figures/main/param/' + mooring(1) + '_verticalStructureFunctionZoom.png');

%% M1: Zoomed-in F(z,t), binned

% ax12 = figure;
% sgtitle('F(z,t): zoomed-in');
% subplot(1,2,1)
% contourf(timeNum_M1(depthToPlotFrom:end,t1:t2)',zmid_M1(depthToPlotFrom:end,t1:t2)',binnedEpsilonM1(depthToPlotFrom:end,t1:t2)','LineColor','none');
% datetick('x','yyyy mmm','keeplimits');
% cmap = flipud(cbrewer2('div','Spectral',NCbar));
% colormap(cmap);
% c = colorbar('YTick',1:14,'TickLabels',{num2str(edges(1)), num2str(edges(2)), num2str(edges(3)), num2str(edges(4)), ... 
%     num2str(edges(5)), num2str(edges(6)), num2str(edges(7)), num2str(edges(8)), num2str(edges(9)), ...
%     num2str(edges(10)), num2str(edges(11)), num2str(edges(12)), num2str(edges(13)), num2str(edges(14))});
% xlabel('time');
% ylabel('depth [m]');
% 
% subplot(1,2,2)
% contourf(timeNum_M1(depthToPlotFrom:end,t3:t4)',zmid_M1(depthToPlotFrom:end,t3:t4)',binnedEpsilonM1(depthToPlotFrom:end,t3:t4)','LineColor','none');
% datetick('x','yyyy mmm','keeplimits');
% cmap = flipud(cbrewer2('div','Spectral',NCbar));
% colormap(cmap);
% c = colorbar('YTick',1:14,'TickLabels',{num2str(edges(1)), num2str(edges(2)), num2str(edges(3)), num2str(edges(4)), ... 
%     num2str(edges(5)), num2str(edges(6)), num2str(edges(7)), num2str(edges(8)), num2str(edges(9)), ...
%     num2str(edges(10)), num2str(edges(11)), num2str(edges(12)), num2str(edges(13)), num2str(edges(14))});
% xlabel('time');
% ylabel('depth [m]');
% 
% savefig('figures/main/param/' + mooring(2) + '_verticalStructureFunctionZoom');
% exportgraphics(ax12,'figures/main/param/' + mooring(2) + '_verticalStructureFunctionZoom.png');

%% clear
% clear depthToPlotFrom t1 t2 t3 t4;

%% IC3: Histogram after binning
ax13 = figure;
histogram(binnedEpsilon);
title('F(z,t): after binning');

% savefig('figures/main/param/' + mooring(1) + '_hist_verticalStructure_binned');
exportgraphics(ax13,'figures/main/param/' + mooring(1) + '_hist_verticalStructure_binned.png');
%% Kurtosis of F(z,t)

% kurtF = kurtosis(verticalStructureIC3Hil,0);
% 
% figure
% plot(time,kurtF);
% xlabel('time');
% ylabel('K_F (Kurtosis-value of F(z,t)');
% 
% 
% lessProneToOutliers = 0;
% moreProneToOutliers = 0;
% 
% for i=1:length(time)
%     if kurtF(i) < 3
%         lessProneToOutliers = lessProneToOutliers + 1;
%     end
%     if kurtF(i) > 3
%         moreProneToOutliers = moreProneToOutliers + 1;
%     end
% end
% 
% propProne = 100*moreProneToOutliers./length(time);

%% ADJUST LIMITS
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [5 2 25 22]);
set(0,'defaultAxesFontSize',18);

%% IC3: Plot of Vertical Structure components T1 and T2

T1 = load('Matfiles/verticalStructureIC3.mat').T1;
T2 = load('Matfiles/verticalStructureIC3.mat').T2;

ax14 = figure;

subplot(2,1,1)
contourf(timeNum_IC3,zmid_IC3,log10(T1),-2.5:-0.05:-4.5,'LineColor','none');
datetick('x','yyyy mmm','keeplimits');
c = colorbar;
c.Label.String = 'log_{10} (F(z,t)) [m^{-1}]';
colormap(flipud(cbrewer2('Spectral')));
% cmocean('tempo');
title('Topography Component');

subplot(2,1,2)
contourf(timeNum_IC3,zmid_IC3,log10(T2),-2.5:-0.05:-4.5,'LineColor','none');
datetick('x','yyyy mmm','keeplimits');
c = colorbar;
c.Label.String = 'log_{10} (F(z,t)) [m^{-1}]';
colormap(flipud(cbrewer2('Spectral')));
% cmocean('tempo');
title('Buoyancy Component');

% sgtitle('IC3: components contributing to F(z,t)');
% savefig('figures/main/param/' + mooring(1) + '_verticalStructureT1_T2');
exportgraphics(ax14,'figures/main/param/' + mooring(1) + '_verticalStructureT1_T2.png');

%% M1: Plot of Vertical Structure components T1 and T2

T1 = load('Matfiles/verticalStructureM1.mat').T1;
T2 = load('Matfiles/verticalStructureM1.mat').T2;

ax15 = figure;

subplot(2,1,1)
contourf(timeNum_M1,zmid_M1,log10(T1),-2.5:-0.1:-4.5,'LineColor','none');
datetick('x','yyyy mmm','keeplimits');
c = colorbar;
c.Label.String = 'log_{10} (F(z,t)) [m^{-1}]';
colormap(flipud(cbrewer2('Spectral')));
% cmocean('tempo');
title('Topography Component');

subplot(2,1,2)
contourf(timeNum_M1,zmid_M1,log10(T2),-2.5:-0.1:-4.5,'LineColor','none');
datetick('x','yyyy mmm','keeplimits');
c = colorbar;
c.Label.String = 'log_{10} (F(z,t)) [m^{-1}]';
colormap(flipud(cbrewer2('Spectral')));
% cmocean('tempo');
title('Buoyancy Component');

% sgtitle('M1: components contributing to F(z,t)');
% savefig('figures/main/param/' + mooring(2) + '_verticalStructureT1_T2');
exportgraphics(ax15,'figures/main/param/' + mooring(2) + '_verticalStructureT1_T2.png');

%% Save parameters for next file

save Matfiles/Fzt.mat verticalStructureIC3Hil verticalStructureM1Hil;