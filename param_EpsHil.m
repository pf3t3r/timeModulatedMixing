clc; clear; close all;
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [5 5 40 15]);
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
% 
% ax1 = figure;
% subplot(1,2,1)
% plot(kurtEpsHil_IC3,zmid_IC3(:,1),'DisplayName','Kurt(\epsilon_{hil}(z,t))');
% hold on
% plot(kurtEpsHillog10_IC3,zmid_IC3(:,i),'DisplayName','Kurt(log(\epsilon_{hil}(z,t)))');
% xline(3,':','DisplayName','Kurtosis = 3');
% hold off
% legend('Location','best');
% xlabel('Kurtosis(\epsilon_{hil}(z,t))');
% ylabel('Depth [m]');
% title('IC3 (2014-2020)');
% 
% subplot(1,2,2)
% plot(kurtEpsHil_M1,zmid_M1(:,1),'DisplayName','Kurt(\epsilon_{hil}(z,t))');
% hold on
% plot(kurtEpsHillog10_M1,zmid_M1(:,i),'DisplayName','Kurt(log(\epsilon_{hil}(z,t)))');
% xline(3,':','DisplayName','Kurtosis = 3');
% hold off
% legend('Location','best');
% xlabel('Kurtosis(\epsilon_{hil}(z,t))');
% ylabel('Depth [m]');
% title('M1 (2014-2020)');
% 
% savefig('figures/main/param/_kurtosisOfEpsHilAcrossDepth');
% exportgraphics(ax1,'figures/main/param/_kurtosisOfEpsHilAcrossDepth.png');

%% IC3: epsHil distribution

% ax2 = figure;
% histogram(sixDayEpsHM_IC3);
% title('\epsilon_{hil} (6dm) (before binning)');
% 
% savefig('figures/main/param/' + mooring(1) + '_epsHil_distribution');
% exportgraphics(ax2,'figures/main/param/' + mooring(1) + '_epsHil_distribution.png');

%% M1: epsHil distribution

% ax3 = figure;
% histogram(sixDayEpsHM_M1);
% title('M1: \epsilon_{hil} (6dm) (before binning)');
% 
% savefig('figures/main/param/' + mooring(2) + '_epsHil_distribution');
% exportgraphics(ax3,'figures/main/param/' + mooring(2) + '_epsHil_distribution.png');

%% epsHil: initialise meshgrid for contour plots

[X2,Y2] = meshgrid(datenum(time'),zmid_IC3(:,1));
[X2a,Y2a] = meshgrid(datenum(time'),zmid_M1(:,1));

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

%% M1 epsHil (logplot)

ax9 = figure;
contourf(X2a,Y2a,log10(eps_hil_M1),-10:0.1:-7,'LineColor','none');
datetick('x','yyyy mmm','keeplimits');
% cmocean('tempo');
colormap(flipud(cbrewer2('Spectral')));
c = colorbar;
c.Label.String = 'log_{10}(\epsilon_{hil}) [log W kg^{-1}]';
ylabel('Depth [m]');
ylim([-inf -100]);
exportgraphics(ax9,'figures/main/param/' + mooring(2) + '_epsHil.png');

%% IC3 epsHil FFT

Ts = 3600;
fs = 1/Ts;
L = length(Ehil6dm_IC3);
f = (0:L/2-1)*fs/L; % per sec
fd = 86400*f; % per day

% M1
P2_M1 = abs(fft(mean(eps_hil_M1(39:42,:)))/L);
P1_M1 = P2_M1(1:floor(L/2));

% IC3
P2_IC3 = abs(fft(mean(eps_hil_IC3(38:41,:)))/L);
P1_IC3 = P2_IC3(1:floor(L/2));

ax10 = figure;
loglog(fd,P1_M1,'DisplayName','M1');
hold on
loglog(fd,P1_IC3,'DisplayName','IC3');
legend();
% xline(fqs(45)/(2*pi),'-','M_2','DisplayName','M_2','FontSize',16,'HandleVisibility','off');
% xline(1/7.5,'-','7.5 Days','DisplayName','7.5 Days','FontSize',16,'HandleVisibility','off');
xline(1/14,'-','14 Days','DisplayName','14 Days','FontSize',16,'HandleVisibility','off','LabelVerticalAlignment','middle');
xline(1/30,'-','30 Days','DisplayName','30 Days','FontSize',16,'HandleVisibility','off','LabelVerticalAlignment','middle');
xline(1/90,'-','90 Days','DisplayName','90 Days','FontSize',16,'HandleVisibility','off','LabelVerticalAlignment','middle');
xline(1/365,'-','1 Year','DisplayName','1 Year','FontSize',16,'HandleVisibility','off','LabelVerticalAlignment','middle');
hold off
xlabel('Frequency [cpd]');
ylabel('Power Density [W^2 m^{-4}]');
% title('FFT: \epsilon_{hil}(t) for IC3 and M1');

exportgraphics(ax10,'figures/main/param/_epsHil_FFT.png');

%% Save parameters for the next file

save Matfiles/epsHil.mat eps_hil_IC3 sixDayEpsHM_IC3 kurtEpsHil_IC3 X2 Y2 ...
    eps_hil_M1 sixDayEpsHM_M1 kurtEpsHil_M1 X2a Y2a ...
    rhoIC3 rhoM1;
