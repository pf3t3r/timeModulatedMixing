% Attempting to combine tidalAnalysis A - C.

clear; close all; clc;
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [5 5 40 20]);

% Add subfolders to path
addpath(genpath(pwd));

load('merged_hourly_unfiltered_data_20142020.mat','lat');
% load IC3-v2.mat;
load IC3.mat;

%% Label Data according to mooring name
mooring = "IC3";

%% Initialise time parameters
Ts = 3600;                                                  % Sampling time [s]  
fs = 2*pi/Ts;                                               % Sampling frequency [/s]
L = length(time_int);

Omega = 7.2921e-5;
fCor = 2*Omega*sin(lat(4)*pi/180);

%% Import tidal parameters needed for LSHA

% Tidal Components
[freqs,names] = frequencies_PF;     % Tidal Frequencies [rad/day]

% Exclude longest periods (Sa+)
freqs = freqs(2:144);
names = names(2:144);

wn = freqs/(24*Ts);                 % Tidal Frequency [rad/s]
T_s = 2*pi./wn;                     % Tidal Period [s]
T_h = T_s/Ts;                       % Tidal Period [h]
T_d = T_s./86400;                   % Tidal Period [d]

%% Check to see which tides we can resolve (3 criteria)

Lm = floor(L/12); % 6m
Lm = floor(L/72); % 1m
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

%% All time periods (every month)

for i=1:72
    [A0(i), Amps(i), Phases(i), Xmean(i), Fit(i), R2(i), resi(i)] = runLSHA(time_int((i-1)*Lm+1:i*Lm),wn,u(:,(i-1)*Lm+1:i*Lm),v(:,(i-1)*Lm+1:i*Lm),noOfMeters);
end

convertingFitStructToArray = cell2mat(struct2cell(Fit));
Fitu = convertingFitStructToArray(1:41,:,:);
Fitv = convertingFitStructToArray(42:end,:,:);

convertingResiStructToArray = cell2mat(struct2cell(resi));
Resiu = convertingResiStructToArray(1:41,:,:);
Resiv = convertingResiStructToArray(42:end,:,:);

clear A0 Amps Phases Xmean Fit R2 resi;

%% FFT on 6-month windows

P2 = [];
P1 = [];
f = [];

for i=1:72
    [P2(i,:,:), P1(i,:,:), f(i,:,:)] = runDFT(noOfMeters,Lm,fs,u(:,(i-1)*Lm+1:i*Lm),v(:,(i-1)*Lm+1:i*Lm),Fitu(:,:,i),Fitv(:,:,i),Resiu(:,:,i),Resiv(:,:,i));
end

f = squeeze(f(1,:,:));

%% Extract residual from FFT file

n = noOfMeters;

% resi
P1_resi_u = P1(:,4*n+1:5*n,:);
P1_resi_v = P1(:,5*n+1:end,:);

%% find depth-averaged residual

P1_resi_u_dptAvg = squeeze(mean(P1_resi_u,2));
P1_resi_v_dptAvg = squeeze(mean(P1_resi_v,2));

%% Make fits

% This code creates a fit for the 12 residuals in u and the 12 residuals in
% v.

pu = [];
pv = [];

% Warning code for following polyfit
% Trying to reduce the degree (from 7 to 5)
for i=1:72
    pu(i,:) = polyfit(f,P1_resi_u_dptAvg(i,:),7);
    pv(i,:) = polyfit(f,P1_resi_v_dptAvg(i,:),7);
    fit_u(i,:) = polyval(pu(i,:),f);
    fit_v(i,:) = polyval(pv(i,:),f);
end

%% test expo func
y = exp(-900000*f)+0.0004;

combinedFitExpU = fit_u;
combinedFitExpV = fit_v;

for i=1:72
    for j=1:floor(Lm/2)
        if fit_u(i,j) < y(j)
            combinedFitExpU(i,j) = y(j);
            combinedFitExpV(i,j) = y(j);
        end
    end
end

% figure
% plot(f,y);

%% Plot DFT according to time period specified by Lm

periodsToPlot = 1:72;

ax1 = figure;
% plot(f,P1_resi_u_dptAvg,'Color',[0.7 0.7 0.7]);
% hold on
plot(f,mean(P1_resi_u_dptAvg));
hold on
% plot(f,y,'Color','green');
% plot(f,combinedFitExpU(periodsToPlot,:),'Color','cyan');
plot(f,mean(combinedFitExpU(periodsToPlot,:)),'Color','cyan');
% plot(f,fit_u,'Color','magenta');
xline(fCor,':','Coriolis Frequency (IC3)');
xline(0.8*fCor,':','near-inertial band');
xline(1.2*fCor,':','near-inertial band');
xline(1.4052e-4,':',{'M_2'}); 
xlim([0.7*fCor 1.3*fCor]);
hold off
title('u, depth-averaged, 6-year mean');

% savefig('figures/main/tideEx/' + mooring + '_u_resi_energy_m1');
exportgraphics(ax1,'figures/main/Extract/' + mooring + '_u_resi_energy_m1.png');

ax2 = figure;
% plot(f,P1_resi_v_dptAvg(periodsToPlot,:),'Color',[0.7 0.7 0.7]);
plot(f,mean(P1_resi_v_dptAvg),'Color',[0.7 0.7 0.7]);
hold on
% plot(f,fit_v,'Color','cyan');
% plot(f,combinedFitExpV(periodsToPlot,:),'Color','cyan');
plot(f,mean(combinedFitExpV),'Color','cyan');
xline(fCor,':','Coriolis Frequency (IC3)');
xline(0.8*fCor,':','near-inertial band');
xline(1.2*fCor,':','near-inertial band');
xline(1.4052e-4,':',{'M_2'}); 
xlim([.7*fCor 1.3*fCor]);
hold off
title('v, depth-averaged, 6-year mean');

% savefig('figures/main/tideEx/' + mooring + '_v_resi_energy_m1');
exportgraphics(ax2,'figures/main/Extract/' + mooring + '_v_resi_energy_m1.png');

NIW_icIT_band = P1_resi_u_dptAvg - fit_u;

for i = 1:72
    for j = 1:367
        if NIW_icIT_band(i,j) < 0
            NIW_icIT_band(i,j) = 0;
        end
    end
end

% NIW_icIT_band = P1_resi_u_dptAvg - combinedFitExpU;
% NIW_icIT_band(NIW_icIT_band<0) = 0;
M2 = 1.4052e-4;

ax3 = figure;
% plot(f,NIW_icIT_band);
% hold on
plot(f,mean(NIW_icIT_band));
hold on
xline(fCor,':','Coriolis Frequency (IC3)');
xline(0.8*fCor,':','near-inertial band');
xline(1.2*fCor,':','near-inertial band');
xline(M2,':',{'M_2'}); 
% xline(0.8*M2,':',{'incoherent semidiurnal band'}); 
% xline(1.2*M2,':',{'incoherent semidiurnal band'}); 
hold off
% xlim([0.7*fCor 1.3*M2]);
xlim([.7*fCor 1.3*fCor]);
title('NIW Energy');

% savefig('figures/main/Extract/' + mooring + '_energyInCombinedIcitNiwBand_m1');
exportgraphics(ax3,'figures/main/Extract/' + mooring + '_energyInCombinedIcitNiwBand_m1.png');

%% Find the area under the curve in NIW_icIT_band

% This defines the NIW band. Note that due to insufficient resolution the
% entire NIW band is not covered. We assume that the NIW energy is evenly
% distributed and that therefore the half-peak we use here is
% representative of NIW energy modulation.

limit1 = find(f>0.798*fCor & f<0.802*fCor);
% limit2 = find(f>1.199*fCor & f<1.202*fCor);
limit2 = find(f>0.997*fCor & f<1.01*fCor);

for i=1:72
    niwIcitEnergy(i) = trapz(f(limit1:limit2),NIW_icIT_band(i,limit1:limit2));
end

% xt = linspace(1,52873,72);
xt = floor(734/2):734:52873;

%% Interpolate to find full time series of values

xq = 1:1:52873;
niwIcitInterp = pchip(xt,niwIcitEnergy,xq);

%% Plot the modulation of NIW/ICIT energy over time

ax4 = figure;
scatter(time(floor(Lm/2):Lm:end),niwIcitEnergy,'MarkerEdgeColor','none','MarkerFaceColor',[0.1 0.5 0.4]);
hold on
plot(time,niwIcitInterp);
hold off
% datetick
% month = {'apr19', 'may19', 'june19', 'july19', 'aug19', 'sept19', 'oct19', 'nov19', 'dec19', 'jan20', 'feb20', 'mar20'};
% set(gca,'xticklabels',month);
% set(gca,'XTick',xt,'XTickLabel',{'14 JASOND','15 JFMAMJ','15 JASOND',...
%     '16 JFMAMJ','16 JASOND','17 JFMAMJ','17 JASOND','18 JFMAMJ','18 JASOND'...
%     '19 JFMAMJ','19 JASOND','20 JFMAMJ','20 JASOND'});
xlabel('Time');
ylabel('Combined near-inertial and incoherent semidiurnal energy [m^2 s^{-2}]');
title('Modulation of WWI');

exportgraphics(ax4,'figures/main/Extract/' + mooring + '_WwiModulation6mth.png');

save Matfiles/niwIcitEnergy.mat niwIcitEnergy niwIcitInterp;