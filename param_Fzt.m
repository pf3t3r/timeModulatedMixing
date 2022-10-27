clc; clear; close all;
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [5 2 40 15]);
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

%% Buoyancy
%% Initialise dz

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

exportgraphics(ax1,'figures/main/param/' + mooring(1) + '_N2_integrated_over_time.png');

%% M1: N^2(t) [vertically-integrated N^2(t)]

ax2 = figure;
plot(time,N2_int_overTime_M1);
xlabel('time');
ylabel('\int N^2 dz (t)');
title('\int N^2 dz (t)');

exportgraphics(ax2,'figures/main/param/' + mooring(2) + '_N2_integrated_over_time.png');

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

%% M1: N(z,t) log plot
ax6 = figure;
contourf(timeNum_M1,zmid_M1,log10(N_M1),-4.1:0.05:-2,'LineColor','none');
datetick('x','yyyy mmm','keeplimits');
% cmocean('tempo');
colormap(flipud(cbrewer2('Spectral')));
c = colorbar;
ylim([-inf -100]);
c.Label.String = 'log_{10}(N(z,t)) [s^{-1}]';
exportgraphics(ax6,'figures/main/param/' + mooring(2) + '_N.png');

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

%% IC3: F(z,t), not binned

ax9a = figure;
contourf(timeNum_IC3,zmid_IC3,log10(verticalStructureIC3Hil),-2.5:-0.1:-4.5,'LineColor','none');
datetick('x','yyyy mmm','keeplimits');
colormap(flipud(cbrewer2('Spectral')));
c = colorbar;
c.Label.String = 'log_{10}(F(z,t)) [m^{-1}]';
ylim([-inf -100]);
ylabel('Depth [m]');
% title('Vertical Structure');
exportgraphics(ax9a,'figures/main/param/' + mooring(1) + '_verticalStructure.png');

%% M1: F(z,t), not binned

ax10a = figure;
contourf(timeNum_M1,zmid_M1,log10(verticalStructureM1Hil),-2.5:-0.1:-4.5,'LineColor','none');
datetick('x','yyyy mmm','keeplimits');
colormap(flipud(cbrewer2('Spectral')));
c = colorbar;
c.Label.String = 'log_{10}(F(z,t)) [m^{-1}]';
ylabel('Depth [m]');
ylim([-inf -100]);
% title('Vertical Structure');
exportgraphics(ax10a,'figures/main/param/' + mooring(2) + '_verticalStructure.png');

%% Kurtosis of F(z,t)
% % MAYBE I will come back to this...
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
title('Topography Component');

subplot(2,1,2)
contourf(timeNum_IC3,zmid_IC3,log10(T2),-2.5:-0.05:-4.5,'LineColor','none');
datetick('x','yyyy mmm','keeplimits');
c = colorbar;
c.Label.String = 'log_{10} (F(z,t)) [m^{-1}]';
colormap(flipud(cbrewer2('Spectral')));
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
title('Topography Component');

subplot(2,1,2)
contourf(timeNum_M1,zmid_M1,log10(T2),-2.5:-0.1:-4.5,'LineColor','none');
datetick('x','yyyy mmm','keeplimits');
c = colorbar;
c.Label.String = 'log_{10} (F(z,t)) [m^{-1}]';
colormap(flipud(cbrewer2('Spectral')));
title('Buoyancy Component');

% sgtitle('M1: components contributing to F(z,t)');
% savefig('figures/main/param/' + mooring(2) + '_verticalStructureT1_T2');
exportgraphics(ax15,'figures/main/param/' + mooring(2) + '_verticalStructureT1_T2.png');

%% Save parameters for next file

save Matfiles/Fzt.mat verticalStructureIC3Hil verticalStructureM1Hil;