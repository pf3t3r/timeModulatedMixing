clc; clear; close all;
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [5 5 40 15]);
set(0,'defaultAxesFontSize',18);

% Add subfolders to path
addpath(genpath(pwd));

%% update each time with MOORING used - used for saving plots

mooring = "IC3";

%% Load data files

load('Mixing_parameterization_fields.mat');                     % Casimir's parametrisation
clear woce_lon woce_lat;
load('merged_hourly_unfiltered_data_20142020.mat', 'lat');      % Load latitude of mooring
load('IC3.mat','time');                                           % Load time
load('niwIcitEnergy.mat');

% Casimir's parametrisation contains the decay scales and depth-integrated
% dissipation rates needed for the time modulation.

%% Normalise the modulation

% Here we modulate the E_wwi signal according to the interpolated 6-month
% energy contained within the combined near-inertial and incoherent
% semidiurnal tidal bands.

% The normalisation factor = niwIcitInterp.
% This must be normalised such that the mean of the time-modulated E_wwi is
% equal to the single-valued E_wwi parametrisation.

tuningFactor = 1./mean((niwIcitInterp - min(niwIcitInterp)) ./ (max(niwIcitInterp) - min(niwIcitInterp)));
normalisation = tuningFactor*(niwIcitInterp - min(niwIcitInterp)) ./ (max(niwIcitInterp) - min(niwIcitInterp));

tf2 = 1./mean((niwIcitEnergy - min(niwIcitEnergy)) ./ (max(niwIcitEnergy) - min(niwIcitEnergy)));
norm2 = tf2*(niwIcitEnergy - min(niwIcitEnergy)) ./ (max(niwIcitEnergy) - min(niwIcitEnergy));

%% E_wwi: Depth-integrated IT dissipation rate due to WWI
% WWI = wave-wave interaction. This is the dissipation process that applies
% to near-inertial waves as well as incoherent tides.

% Parametrised energy input to mixing by near-inertial waves (Alford, 2020)
% W/m2
alfordIC3 = 3.3e-4;

% E_wwi_IC3_mod = normalisation*power_wwi(657,279);
E_wwi_IC3_mod = normalisation*alfordIC3;
% E_wwi_IC3_singlePoints = norm2*power_wwi(657,279);
E_wwi_IC3_singlePoints = norm2*alfordIC3;

E_wwi_IC3_mod_6dm = movmean(E_wwi_IC3_mod,149);

stdE = std(E_wwi_IC3_mod_6dm);
kurtE = kurtosis(E_wwi_IC3_mod_6dm);

% timePoints = [time(floor(0.5*52873/12)),time(floor(1.5*52873/12)),...
%     time(floor(2.5*52873/12)),time(floor(3.5*52873/12)),time(floor(4.5*52873/12)),...
%     time(floor(5.5*52873/12)),time(floor(6.5*52873/12)),time(floor(7.5*52873/12)),...
%     time(floor(8.5*52873/12)),time(floor(9.5*52873/12)),time(floor(10.5*52873/12)),...
%     time(floor(11.5*52873/12))];

% E_wwi(t) incl. mean and std
ax3 = figure;
plot(time,E_wwi_IC3_mod_6dm,'DisplayName','<E_{wwi}(t)> (6-day mean)');
hold on
% scatter(timePoints,E_wwi_IC3_singlePoints,'DisplayName','<E_{wwi}(t)> (non-interpolated)');
% yline(power_wwi(657,279),'-','DisplayName','E_{wwi}');
yline(alfordIC3,'-','DisplayName','E_{niw}');
% yline(mean(E_wwi_IC3_mod_6dm),'--','DisplayName','Mean <E_{hil(t)}>','LineWidth',3);
% yline(mean(E_wwi_IC3_mod_6dm)+stdE,'--','DisplayName','\mu + \sigma','Alpha',0.3);
% yline(mean(E_wwi_IC3_mod_6dm)-stdE,'--','DisplayName','\mu - \sigma','Alpha',0.3);
hold off
legend();
% xlabel('time');
ylabel('E_{wwi} [W m^{-2}]');
% title('depth-integrated dissipation due to wave-wave interaction');

% savefig('figures/main/param/' + mooring + '_EWwiModulation-6day-mean-sigma');
exportgraphics(ax3,'figures/main/param/' + mooring + '_EWwiModulation-6day-mean-sigma.png');

%% E_wwi(t): zoomed-in

t1 = 14300;
t2 = 15700;

ax4 = figure;
subplot(1,2,1)
plot(time,E_wwi_IC3_mod_6dm,'DisplayName','<E_{wwi}(t)> (6dm)','LineWidth',1.5);
legend();
xlabel('time');
ylabel('E_{wwi} [W m^{-2}]');
title('<E_{wwi}(t)> (6-day mean)');

subplot(1,2,2)
plot(time(t1:t2),E_wwi_IC3_mod_6dm(t1:t2),'DisplayName','<E_{wwi}(t)> (6dm)','LineWidth',1.5);
legend();
xlabel('time');
ylabel('E_{wwi} [W m^{-2}]');
title('<E_{wwi}(t)> (6-day mean): close-up');

clear time1 time2;

% savefig('figures/main/param/' + mooring + '_EWwi-subplots');
exportgraphics(ax4,'figures/main/param/' + mooring + '_EWwi-subplots.png');

%% V2. Histogram. Distribution of E.

ax5 = figure;
h = histogram(E_wwi_IC3_mod_6dm);
hold on
xline(mean(E_wwi_IC3_mod_6dm),':','\mu','LineWidth',1.5);
xline(mean(E_wwi_IC3_mod_6dm)+stdE,':','\mu + \sigma','LineWidth',1.5);
xline(mean(E_wwi_IC3_mod_6dm)-stdE,':','\mu - \sigma','LineWidth',1.5);
hold off
xlabel('E_{wwi}(t) [W m^{-2}]');
ylabel('No. of values');
title('Distribution of E_{wwi}(t)');

% savefig('figures/main/param/' + mooring + '_EWwi-histogram');
exportgraphics(ax5,'figures/main/param/' + mooring + '_EWwi-histogram.png');

%% Fourier Analysis of E_wwi.

Ts = 3600;
fs = 2*pi/Ts;
L = length(time);

fftEWwi = fft(E_wwi_IC3_mod_6dm);
P2 = abs(fftEWwi/L);
P1 = P2(1:floor(L/2));
f = (0:L/2-1)*fs/L;

[freqs,names] = frequencies_PF;
Omega = 7.2921e-5;

f_M2 = freqs(45)/(24*3600);
f_M4 = freqs(76)/(24*3600);
f_S1 = freqs(20)/(24*3600);
f_S3 = freqs(67)/(24*3600);
fC = 2*Omega*sin(deg2rad(lat(4)));

ax6 = figure;
semilogy(f,P1,'DisplayName','<E_{wwi}(t)>');
hold on
xline(f_M2,':',names(45),'DisplayName','M2');
xline(f_M4,':',names(76),'DisplayName','M4');
xline(fC,':','f','DisplayName','f');
xline(f_S1,':',names(20),'DisplayName','S1');
xline(f_S3,':',names(67),'DisplayName','S3');
hold off
legend();
xlabel('frequency [s^{-1}]');
ylabel('|<E_{hil}(t)>|^{2} [W^2 m^{-4}]');
title('DFT: U component of E');

% savefig('figures/main/param/' + mooring + '_EWwi-fft');
exportgraphics(ax6,'figures/main/param/' + mooring + '_EWwi-fft.png');

%% Evaluate how much of an effect smoothing has on the mean
windowSizes = 6:143:42906;
for k = 1 : length(windowSizes)
  smoothedSignal = movmean(E_wwi_IC3_mod, windowSizes(k));
  sad(k) = sum(abs(smoothedSignal - E_wwi_IC3_mod));
end

figure
plot(windowSizes, sad, 'b*-', 'LineWidth', 2);
xlabel('Window Size');
ylabel('SAD');

%% Save parameters for the next file

save Matfiles/E_wwi.mat E_wwi_IC3_mod_6dm E_wwi_IC3_mod;