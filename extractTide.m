% Add subfolders to path
addpath(genpath(pwd));

clear; close all; clc;
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [5 5 40 15]);

load('Data/merged_hourly_unfiltered_data_20142020.mat','lat');
load Matfiles/M1.mat;
load Matfiles/IC3.mat;

%% Label Data according to mooring name

mooring = ["IC3","M1"];
noOfMeters = [noOfMeters,noOfMetersM1];

%% Initialise time parameters
Ts = 3600;                                                  % Sampling time [s]  
% fs = 2*pi/Ts;                                               % Sampling frequency [/s]
L = length(time_int);

Omega = 7.2921e-5;
fCor = 2*Omega*sin(lat(4)*pi/180);

%% Import tidal parameters needed for LSHA

% Tidal Components
[freqs,names] = frequencies_PF;     % Tidal Frequencies [rad/day]

% Exclude longest periods (Sa+)
% Excluding even more ... now only Msf - HH
freqs = freqs(5:144);
names = names(5:144);

wn = freqs/(24*Ts);                 % Tidal Frequency [rad/s]
T_s = 2*pi./wn;                     % Tidal Period [s]
T_h = T_s/Ts;                       % Tidal Period [h]
T_d = T_s./86400;                   % Tidal Period [d]

%% Check to see which tides we can resolve (3 criteria)

% Lm = floor(L/12); % 6m
Lm = L; % Total time series
timeStepsLSHA = 1:Lm;

% 1. The longest period component that can be resolved is
% less than the extent of the time series. The longest period component in
% the above file is the Sa tide which has a period of one year. So this is
% technically OK but it is excluded (for now) because of other issues
% (contamination with the mesoscale signal).

for i=1:length(wn)
    if T_d(i) > Lm*(1/24)
        output = sprintf('%s can not be resolved and must be removed',string(names(i)));
        disp(output);
    end
end


% 2. Rayleigh Criterion: throw out frequencies that are too close to one
% another.
% t_obs = time_int; % is this correct?
breaksRayleighCriterion = zeros(length(wn),1);
for i=1:length(wn)-1
    if (wn(i+1) - wn(i))*(timeStepsLSHA*Ts) > 2*pi
        breaksRayleighCriterion(i+1) = 1;
    end
end

% If all values in breaksRayleighCriterion are zero, then the time series
% is long enough to resolve every component from its neighbour and none of
% them need to be removed. Any values equal to one break the Rayleigh
% Criterion and must be removed from the analysis.

% 3. Nyquist Criterion: throw out frequencies that are larger than the
% Nyquist Frequency fN. fN = (0.5)*fS.
T_fN = 2*Ts;
breaksNyquistCriterion = zeros(length(wn),1);

% Any values in breaksNyquistCriterion that are equal to one break the
% Nyquist criterion and must be removed from the analysis.

for i=1:length(wn)
    if 2*pi/wn(i) < T_fN
        breaksNyquistCriterion(i) = 1;
    end
end

% This code will show if any of the above three criteria have not been met.
if any(breaksRayleighCriterion,1)
    disp('Rayleigh Criterion not met.');
elseif any(breaksNyquistCriterion,1)
    disp('Nyquist Criterion not met.');
elseif any(T_s) > time_int(end)
    disp('Some components not resolved');
else
    disp('Rayleigh and Nyquist Criteria have been met for all components.')
    disp('The time series can resolve all components.')
end

%% Run the LSHA on full time series for u and v at IC3 and M1; extract residual for both

% Here any NaNs in the data will be converted into zeros by runLSHA.

% IC3
[A0,Amps,Phases,Xmean,Fit,R2,resi] = runLSHA(time_int,wn,u,v,noOfMeters(1));

% M1
[~,~,~,~,Fit_m1,~,~] = runLSHA(time_int,wn,u_M1',v_M1',30);

save Matfiles/IC3.mat Fit -append;
save Matfiles/M1.mat Fit_m1 -append;

%% Create barotropic tide

% Interpolation method is faulty => replace with method whereby I
% vertically-average across the depths of the original measurements.
% I may need to weight this averaging according to the depth.

bathy_M1 = 1712;
bathy_IC3 = 1635;

weighting_M1 = zeros(length(depths_M1),1);
weighting_M1(1) = 10;
for i = 1:length(depths_M1)-2
    weighting_M1(i+1) = (depths_M1(i+2) + depths_M1(i+1))/2 - (depths_M1(i+1) + depths_M1(i))/2;
end
% weighting_M1(end-1) = (depths_M1(end-1) + depths_M1(end))/2 - (depths_M1(end-2) + depths_M1(end-1))/2;
weighting_M1(end) = (bathy_M1 + depths_M1(end))/2 - (depths_M1(end) + depths_M1(end-1))/2;

weighting_IC3 = zeros(length(depths),1);
weighting_IC3(1) = 10;
for i = 1:length(depths)-2
    weighting_IC3(i+1) = (depths(i+2) + depths(i+1))/2 - (depths(i+1) + depths(i))/2;
end
weighting_IC3(end) = (bathy_IC3 + depths(end))/2 - (depths(end) + depths(end-1))/2; 

% Find the weighted average
for i = 1:length(weighting_M1)
    avgu(i,:) = weighting_M1(i)*(Fit_m1.u(i,:));
    avgv(i,:) = weighting_M1(i)*(Fit_m1.v(i,:));
end
barU_M1 = mean(avgu/(sum(weighting_M1/length(weighting_M1))));
barV_M1 = mean(avgv/(sum(weighting_M1/length(weighting_M1))));

for i = 1:length(weighting_IC3)
    avgu_i(i,:) = weighting_IC3(i)*(Fit.u(i,:));
    avgv_i(i,:) = weighting_IC3(i)*(Fit.v(i,:));
end
barU_IC3 = mean(avgu_i/(sum(weighting_IC3/length(weighting_IC3))));
barV_IC3 = mean(avgv_i/(sum(weighting_IC3/length(weighting_IC3))));



% barotropicU = mean(Fit.u);
barotropicU = barU_IC3;
% barotropicV = mean(Fit.v);
barotropicV = barV_IC3;
% barotropicU_M1 = mean(Fit_m1.u);
barotropicU_M1 = barU_M1;
% barotropicV_M1 = mean(Fit_m1.v);
barotropicV_M1 = barV_M1;

% The velocity vector U is created by the zonal and meridional vectors 
% bar_u and bar_v.
theta2 = atan2(barotropicV,barotropicU);
theta2_M1 = atan2(barotropicV_M1,barotropicU_M1);

% Find the velocity vector U that is the vector sum of bar_u and bar_v.
U = barotropicU./cos(theta2);
U_M1 = barotropicU_M1./cos(theta2_M1);

% Ensure consistency by cross-checking
U_check = barotropicV./sin(theta2);
U_check_M1 = barotropicV_M1./sin(theta2_M1);
for i = 1:8426
    assert((U(i) - U_check(i))^2 < 1e-10);
    assert((U_M1(i) - U_check_M1(i))^2 < 1e-10);
end

for i=1:length(time)
    if (U(i) - U_check(i))^2 > 1e-25
        disp(i)
    end
    if (U_M1(i) - U_check_M1(i))^2 > 1e-25
        disp(i)
    end
end

save Matfiles/IC3.mat U -append;
save Matfiles/M1.mat U_M1 -append;

%% M1

load Matfiles/M1.mat;
U_M1_smooth = movmean(U_M1,149);

set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [5 5 40 15]);
set(0,'defaultAxesFontSize',16);

ax1 = figure;
plot(time,100*U_M1,'LineWidth',1.2,'DisplayName','U_{M1}');
hold on
h1 = plot(time,100*U_M1_smooth,'LineWidth',1.2,'Color',[0.6 0.6 0.6],'DisplayName','U_{M1} (6DM)');
hold off
% h1.Color(4) = 0.6;
legend('Location','best');
ylabel('U [cm s^{-1}]');

exportgraphics(ax1,'figures/main/param/M1_barotropicTide.png');

save Matfiles/M1.mat U_M1_smooth -append;

%% IC3 
% gotta create it for IC3 at the same time!

% load('IC3.mat','U','time');

U_IC3 = load("Matfiles/IC3.mat").U;
U_IC3_smooth = movmean(U,149);

ax2 = figure;
plot(time,100*U,'LineWidth',1.2,'DisplayName','U_{IC3}');
hold on
h2 = plot(time,100*U_IC3_smooth,'LineWidth',1.2,'Color',[0.6 0.6 0.6],'DisplayName','U_{IC3} (6DM)');
hold off
% h2.Color(4) = 0.6;
ylabel('U [cm s^{-1}]');
legend('Location','best');
exportgraphics(ax2,'figures/main/param/IC3_barotropicTide.png');

save Matfiles/IC3.mat U_IC3_smooth -append;

%% FFT for U

bar2 = mean(Fit.u);

Ts = 3600;
fs = 1/Ts;
L = length(U);
f = (0:L/2-1)*fs/L; % per sec
fd = 86400*f; % per day

% M1
P2_M1 = abs(fft(U_M1)/L);
P1_M1 = P2_M1(1:floor(L/2));

% IC3
P2_IC3 = abs(fft(U)/L);
P1_IC3 = P2_IC3(1:floor(L/2));

ax2a = figure;
loglog(fd,P1_M1,'DisplayName','M1','LineWidth',1.5);
hold on
loglog(fd,P1_IC3,'DisplayName','IC3','LineWidth',1.5);
legend();
xline(freqs(72)/(2*pi),'--','M_4','DisplayName','M_4','FontSize',16,'HandleVisibility','off');
xline(freqs(41)/(2*pi),'--','M_2','DisplayName','M_2','FontSize',16,'HandleVisibility','off');
xline(freqs(17)/(2*pi),'--','K_1','DisplayName','K_1','FontSize',16,'HandleVisibility','off');
% xline(freqs(1)/(2*pi),'--','S_{sa}','DisplayName','S_{sa}','FontSize',16,'HandleVisibility','off');
% xline(freqs(5)/(2*pi),'--','M_{f}','DisplayName','M_{f}','FontSize',16,'HandleVisibility','off');
xline(1/7,'--','7 Days','DisplayName','7 Days','FontSize',16,'HandleVisibility','off');
xline(1/14,'--','14 Days','DisplayName','14 Days','FontSize',16,'HandleVisibility','off');
xline(1/30,'--','1 Month','DisplayName','1 Month','FontSize',16,'HandleVisibility','off');
xline(1/365.25,'--','1 Year','DisplayName','1 Year','FontSize',16,'HandleVisibility','off');
hold off
xlabel('Frequency [cpd]');
ylabel('Power Density [W^2 m^{-4}]');
% title('FFT: U(t) for IC3 and M1');

exportgraphics(ax2a,'figures/main/param/_fft_barotropicTide.png');

%% Kurtosis

kurt_U_IC3 = kurtosis(U);
kurt_U_M1 = kurtosis(U_M1);