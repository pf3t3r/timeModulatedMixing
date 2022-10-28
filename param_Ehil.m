clc; clear; close all;
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [5 2 40 15]);
set(0,'defaultAxesFontSize',18);

% Add subfolders to path
addpath(genpath(pwd));

%% LOAD parametrisation (de Lavergne et al, 2020) & mooring latitude

load('Data/Mixing_parameterization_fields.mat');                 % Parametrisation
clear woce_lon woce_lat;
load('Data/merged_hourly_unfiltered_data_20142020.mat', 'lat');  % Mooring latitude
% load('WOCE_climatology.mat','zax');                       % Load zax

% DEPTH-INTEGRATED POWER FROM PARAMETRISATION [De Lavergne et al, 2020]
powerhilM1 = power_hil(660,279);
powerhilIC3 = power_hil(657,279);

% ALTERNATIVE CALCULATION
hrms_IC3 = 146;
hrms_M1 = 180;
Ehil_alt_IC3 = mean(N2_bot_IC3)*mean(U)*hrms_IC3.^2;
Ehil_alt_M1 = mean(N2_bot_M1)*mean(U_M1)*hrms_M1.^2;

%% Some notes
% Casimir's parametrisation contains the decay scales and depth-integrated
% dissipation rates needed for the time modulation.

% NOTE that for the time modulation based on the barotropic tide, we use
% eps_hil (which is approximately equal to eps_bot). For the other
% processes, i.e. near-inertial waves and incoherent tide, we will use
% eps_wwi.

% Note also that the bottom buoyancy frequency estimate is based on the
% vertical average of the three bottom-most estimates of the buoyancy
% frequency. This figure was chosen in order to align with the average
% value for H_bot (150m). An alternative method is to use the midpoint
% value of the two lowest sensors. Sensitivity analysis may need to be
% carried out on N_bot.

%% LOAD barotropic tide U, buoyancy frequency N, depths, and pressures

% MOORING NAME: refer to this for saving figures
mooring = ["IC3","M1"];

% BAROTROPIC TIDE for IC3 and M1
M1 = load('Matfiles/M1.mat');
U_M1 = M1.U_M1;
IC3 = load('Matfiles/IC3.mat');
U_IC3 = IC3.U;

% TIME
time = M1.time;

% BUOYANCY, SALINITY, TEMPERATURE, PRESSURE for IC3 and M1
load('Matfiles/N_is.mat');               % Load N^2 and N^2_bot (IC3 and M1)
load('Matfiles/SA_CT_interpolated.mat'); % Load interpolated and stabilised SA and CT (IC3 and M1)
load('Matfiles/p_i.mat');                % Load interpolated pressure (IC3 and M1)

% Depths
depths_IC3 = IC3.depths;
depths_M1 = M1.depths_M1;

%% Check stratification

anySuperLowN2inIC3 = ~find(any(N2_IC3<1e-11));
anySuperLowN2inM1 = ~find(any(N2_M1<1e-11));

if isempty(anySuperLowN2inM1)
    disp('N2 in M1 OK...');
end

if isempty(anySuperLowN2inIC3)
    disp('N2 in IC3 OK...');
end

clear anySuperLowN2inIC3 anySuperLowN2inM1;

%% Show tide

figure
plot(time,U_IC3*100);
ylabel('U [cms^{-1}]');

figure
plot(time,U_M1*100);
ylabel('U [cms^{-1}]');

%% M1: find E_hil
[Ehil_M1, Ehil6dm_M1, E_hil_M1_Nbot2,E_hil_M1_Nbot2_6dm,E_hil_M1_U,E_hil_M1_U_6dm,stdE_M1,kurtE_M1] = EhilCalculator(U_M1',N2_M1,N2_bot_M1,powerhilM1);

%% IC3: find E_hil
[Ehil_IC3,Ehil6dm_IC3,E_hil_IC3_Nbot2,E_hil_IC3_Nbot2_6dm,E_hil_IC3_U,E_hil_IC3_U_6dm,stdE_IC3,kurtE_IC3] = EhilCalculator(U_IC3',N2_IC3,N2_bot_IC3,powerhilIC3);

%% Sensitivity Analysis: Approach (1)

% In this analysis we check what happens to vertically-integrated
% dissipation rate if Nbot is double or half its initial value. We expect
% that since the signal is normalised that this will have a negligible
% effect.

% Double Nbot for IC3 and M1
[Eh_IC3_Dbl,~,~,~,~,~,~,~] = EhilCalculator(U_IC3',N2_IC3,2*N2_bot_IC3,powerhilIC3);
[Eh_M1_Dbl, ~, ~,~,~,~,~,~] = EhilCalculator(U_M1',N2_M1,2*N2_bot_M1,powerhilM1);

figure
% plot(Ehil6dm_IC3-Eh6dm_IC3_Dbl);
plot(Ehil_IC3-Eh_IC3_Dbl);
title('IC3 Sensitivity Analysis I.a: 2x Nbot');

figure
% plot(Ehil6dm_M1-Eh6dm_M1_Dbl);
plot(Ehil_M1 - Eh_M1_Dbl);
title('M1 Sensitivity Analysis I.a: 2x Nbot');

clear Eh_IC3_Dbl Eh6dm_IC3_Dbl Eh_IC3_Nbot2_Dbl Eh_IC3_Nbot2_6dm_Dbl Eh_IC3_U_Dbl Eh_IC3_U_6dm_Dbl stdE_IC3_Dbl kurtE_IC3_Dbl;
clear Eh_M1_Dbl Eh6dm_M1_Dbl Eh_M1_Nbot2_Dbl Eh_M1_Nbot2_6dm_Dbl Eh_M1_U_Dbl Eh_M1_U_6dm_Dbl stdE_M1_Dbl kurtE_M1_Dbl;

% Halve Nbot for IC3 and M1
[Eh_IC3_Hlf,~,~,~,~,~,~,~] = EhilCalculator(U_IC3',N2_IC3,0.5*N2_bot_IC3,powerhilIC3);
[Eh_M1_Hlf, ~, ~,~,~,~,~,~] = EhilCalculator(U_M1',N2_M1,0.5*N2_bot_M1,powerhilM1);

figure
% plot(Ehil6dm_IC3-Eh6dm_IC3_Hlf);
plot(Ehil_IC3 - Eh_IC3_Hlf);
title('IC3 Sensitivity Analysis I.b: 0.5x Nbot');

figure
% plot(Ehil6dm_M1-Eh6dm_M1_Hlf);
plot(Ehil_M1 - Eh_M1_Hlf);
title('M1 Sensitivity Analysis I.b: 0.5x Nbot');

clear Eh_IC3_Hlf Eh6dm_IC3_Hlf Eh_IC3_Nbot2_Hlf Eh_IC3_Nbot2_6dm_Hlf Eh_IC3_U_Hlf Eh_IC3_U_6dm_Hlf stdE_IC3_Hlf kurtE_IC3_Hlf;
clear Eh_M1_Hlf Eh6dm_M1_Hlf Eh_M1_Nbot2_Hlf Eh_M1_Nbot2_6dm_Hlf Eh_M1_U_Hlf Eh_M1_U_6dm_Hlf stdE_M1_Hlf kurtE_M1_Hlf;

%% Sensitivity Analysis: Approach (2)

% I may need to rewrite sections in extractBuoyancy for this.

% First value is the DEFAULT
% Other values are monotonically increasing depth ranges from minimum
% (40:41) up to a maximum of (x:41)
% We choose x by visually inspecting Hovmoller Diagrams of N(z,t) and
% determining a point where there is a significant shift in N.
% We find that this occurs for x = 26 in IC3 and x = 31 in M1.
possibleN2botsIC3 = [mean(N2_IC3(40:41,:)); mean(N2_IC3(39:41,:));...
    mean(N2_IC3(38:41,:));mean(N2_IC3(37:41,:));mean(N2_IC3(36:41,:));...
    mean(N2_IC3(35:41,:));mean(N2_IC3(34:41,:));mean(N2_IC3(33:41,:));...
    mean(N2_IC3(32:41,:));mean(N2_IC3(31:41,:));mean(N2_IC3(30:41,:));...
    mean(N2_IC3(29:41,:));mean(N2_IC3(28:41,:));mean(N2_IC3(27:41,:));...
    mean(N2_IC3(26:41,:))];

possibleN2botsM1 = [mean(N2_M1(40:41,:));mean(N2_M1(39:41,:));...
    mean(N2_M1(38:41,:));mean(N2_M1(37:41,:));mean(N2_M1(36:41,:));...
    mean(N2_M1(35:41,:));mean(N2_M1(34:41,:));mean(N2_M1(33:41,:));...
    mean(N2_M1(32:41,:));mean(N2_M1(31:41,:))];

%%
% 's' for sensitivity analysis
for i=1:length(possibleN2botsIC3(:,1))
    [Ehil_IC3s(i,:), Ehil6dm_IC3s(i,:), ~,~,~,~,stdE_IC3s(i),kurtE_IC3s(i)] = EhilCalculator(U_IC3',N2_IC3,possibleN2botsIC3(i,:),powerhilIC3);
end

for i = 1:length(possibleN2botsM1(:,1))
    [Ehil_M1s(i,:),Ehil6dm_M1s(i,:),~,~,~,~,stdE_M1s(i),kurtE_M1s(i)] = EhilCalculator(U_M1',N2_M1,possibleN2botsM1(i,:),powerhilM1);
end

%% SAD

Ehil_IC3s(1,1) = Ehil_IC3s(1,2);
Ehil_IC3s(2,1) = Ehil_IC3s(2,2);

for i=1:15
    testCheck(i,:) = ((Ehil_IC3s(1,:) - Ehil_IC3s(i,:))./Ehil_IC3s(1,:))*100;
    testCheckMean(i) = mean(testCheck(i,:),'omitnan');
end

%%
figure
plot(time,Ehil_IC3s);

figure;
yyaxis left;
plot(kurtE_IC3s,'DisplayName','Kurt(E)');hold on;
yyaxis right;
plot(stdE_IC3s,'DisplayName','SD(E)');
title('IC3: STD and Kurtosis of E_{hil} when N_{bot} is defined over different ranges');

figure;
yyaxis left;
plot(kurtE_M1s,'DisplayName','Kurt(E)');hold on;
yyaxis right;
plot(stdE_M1s,'DisplayName','SD(E)');
title('M1: STD and Kurtosis of E_{hil} when N_{bot} is defined over different ranges');

%%
monthTicks = [457];
for i = 1:70
    a = 457 + i*744;
    monthTicks = [monthTicks a];
end

%% M1: E_hil(t) with mean and STD

ax1 = figure;
plot(time,Ehil6dm_M1,'DisplayName','<E_{hil}(t)> (6DM)','LineWidth',1.2);
hold on
% yline(powerhilM1,'-','DisplayName','E_{hil}');
yline(mean(Ehil6dm_M1),'--','DisplayName','mean <E_{hil}(t)>','LineWidth',3);
% yline(mean(Ehil6dm_M1)+stdE_M1,'--','DisplayName','\mu + \sigma','Alpha',0.5);
% yline(mean(Ehil6dm_M1)-stdE_M1,'--','DisplayName','\mu - \sigma','Alpha',0.5);
hold off
% xticks(time(monthTicks));
% xtickformat('yy mm');
% ax = gca;
% ax.XAxis.MinorTickValues = time(monthTicks);
% ax.XAxis.MinorTick = 'on';
legend();
% xlabel('time');
ylabel('Energy Available for Mixing E_{hil} [W m^{-2}]');
% title('M1: depth-integrated internal tide dissipation rate');

savefig('figures/main/param/' + mooring(2) + '_EHil-6dm-sigma');
exportgraphics(ax1,'figures/main/param/' + mooring(2) + '_EHil-6dm-sigma.png');


%% IC3: E_hil(t) with mean and STD

ax2 = figure;
plot(time,Ehil6dm_IC3,'DisplayName','<E_{hil}(t)> (6DM)','LineWidth',1.2);
hold on
% yline(powerhilIC3,'-','DisplayName','E_{hil}');
yline(mean(Ehil6dm_IC3),'--','DisplayName','mean <E_{hil}(t)>','LineWidth',3);
% yline(mean(Ehil6dm_IC3)+stdE_IC3,'--','DisplayName','\mu + \sigma','Alpha',0.5);
% yline(mean(Ehil6dm_IC3)-stdE_IC3,'--','DisplayName','\mu - \sigma','Alpha',0.5);
hold off
legend();
% xlabel('time');
ylabel('Energy Available for Mixing E_{hil} [W m^{-2}]');
% title('IC3: depth-integrated internal tide dissipation rate');

savefig('figures/main/param/' + mooring(1) + '_EHil-6dm-sigma');
exportgraphics(ax2,'figures/main/param/' + mooring(1) + '_EHil-6dm-sigma.png');

%% FFT for M1 and IC3

[fqs,names] = frequencies_PF;

Ts = 3600;
fs = 1/Ts;
L = length(Ehil6dm_M1);
f = (0:L/2-1)*fs/L; % per sec
fd = 86400*f; % per day

% M1
P2_M1 = abs(fft(Ehil6dm_M1)/L);
P1_M1 = P2_M1(1:floor(L/2));

% IC3
P2_IC3 = abs(fft(Ehil6dm_IC3)/L);
P1_IC3 = P2_IC3(1:floor(L/2));

ax2a = figure;
loglog(fd,P1_M1,'DisplayName','M1');
hold on
loglog(fd,P1_IC3,'DisplayName','IC3');
legend();
xline(fqs(45)/(2*pi),'-','M_2','DisplayName','M_2','FontSize',16,'HandleVisibility','off');
xline(1/7.5,'-','7.5 Days','DisplayName','7.5 Days','FontSize',16,'HandleVisibility','off');
xline(1/14,'-','14 Days','DisplayName','14 Days','FontSize',16,'HandleVisibility','off');
xline(1/50,'-','50 Days','DisplayName','50 Days','FontSize',16,'HandleVisibility','off');
hold off
xlabel('Frequency [cpd]');
ylabel('Power Density [W^2 m^{-4}]');
title('FFT: E_{hil}(t) for IC3 and M1');

exportgraphics(ax2a,'figures/main/param/_Ehil_FFT.png');


%% define time segment for comparison of N_bot^2 and U

t1 = 14300;
t2 = 15700;

%% Update fontsize for comparative influence plots
set(0,'defaultAxesFontSize',16);

%% IC3: E_hil(t), comparative influence of N_bot^2 and U

ax3 = figure;
subplot(3,2,1)
plot(time,Ehil6dm_IC3,'DisplayName','<E_{hil}(t)> (6dm)','LineWidth',1.5);
legend('Location','northeast','FontSize',10);
ylabel('E_{hil} [W m^{-2}]');
% title('<E_{hil}(t)> (6-day mean)');

subplot(3,2,3)
plot(time,E_hil_IC3_Nbot2_6dm,'DisplayName','<N_{bot}^2(t)> (6dm)','Color','black','LineWidth',1.5);
legend('Location','northeast','FontSize',10);
ylabel('E_{hil} [W m^{-2}]');
% title('<N_{bot}^2(t)> (6-day mean)');

subplot(3,2,5)
plot(time,E_hil_IC3_U_6dm,'DisplayName','<U(t)> (6dm)','Color','green','LineWidth',1.5);
legend('Location','northeast','FontSize',10);
ylabel('E_{hil} [W m^{-2}]');
% title('<U(t)> (6-day mean)');

subplot(3,2,2)
plot(time(t1:t2),Ehil6dm_IC3(t1:t2),'DisplayName','<E_{hil}(t)> (6dm)','LineWidth',1.5);
legend('Location','northeast','FontSize',10);

% title('<E_{hil}(t)> (6-day mean): close-up');

subplot(3,2,4)
plot(time(t1:t2),E_hil_IC3_Nbot2_6dm(t1:t2),'DisplayName','<N_{bot}^2(t)> (6dm)','Color','black','LineWidth',1.5);
legend('Location','northeast','FontSize',10);
% title('<N_{bot}^2(t)> (6-day mean): close-up');

subplot(3,2,6)
plot(time(t1:t2),E_hil_IC3_U_6dm(t1:t2),'DisplayName','<U(t)> (6dm)','Color','green','LineWidth',1.5);
legend('Location','northeast','FontSize',10);
% title('<U(t)> (6-day mean): close-up');

% sgtitle('IC3: Influence of N^2 and U on E_{hil}');
savefig('figures/main/param/' + mooring(1) + '_EHil-subplots');
exportgraphics(ax3,'figures/main/param/' + mooring(1) + '_EHil-subplots.png');

%% same as above but without zoom
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [5 2 40 22]);

ax3a = figure;
subplot(3,1,1)
plot(time,Ehil6dm_IC3,'DisplayName','<E_{hil}(t)> (6DM)','LineWidth',1.5);
legend();
ylabel({'Energy Available for';' Mixing E_{hil} [Wm^{-2}]'});
% title('<E_{hil}(t)>');

subplot(3,1,2)
plot(time,E_hil_IC3_Nbot2_6dm,'DisplayName','<N_{bot}^2(t)> (6DM)','Color','black','LineWidth',1.5);
legend();
ylabel({'Energy Available for';' Mixing E_{hil} [Wm^{-2}]'});
% title('<N_{bot}^2(t)>');

subplot(3,1,3)
plot(time,E_hil_IC3_U_6dm,'DisplayName','<U(t)> (6DM)','Color','green','LineWidth',1.5);
legend();
ylabel({'Energy Available for';' Mixing E_{hil} [Wm^{-2}]'});
% title('<U(t)>');

% sgtitle('IC3: Influence of N^2 and U on E_{hil}');
savefig('figures/main/param/' + mooring(1) + '_EHil-subplotsNoZoom');
exportgraphics(ax3a,'figures/main/param/' + mooring(1) + '_EHil-subplotsNoZoom.png');

%% AS ABOVE but in ONE PLOT
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [5 2 40 22]);

ax3aa = figure;
plot(time,Ehil6dm_IC3,'DisplayName','<E_{hil}(t)> (6DM)','LineWidth',1.5);
hold on
plot(time,E_hil_IC3_Nbot2_6dm,'DisplayName','<N_{bot}^2(t)> (6DM)','Color','black','LineWidth',1);
plot(time,E_hil_IC3_U_6dm,'DisplayName','<U(t)> (6DM)','Color','green','LineWidth',1);
hold off
legend();
ylabel({'Energy Available for';' Mixing E_{hil} [Wm^{-2}]'});
exportgraphics(ax3aa,'figures/main/param/' + mooring(1) + '_EHil_drivers.png');

%% IC3 comp influence of Nbot2 and U: zoom

set(0,'defaultAxesFontSize',15);

ax3b = figure;
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [5 2 40 15]);
ax3b = figure;
subplot(2,1,1)
yyaxis left
plot(time(t1:t2),Ehil6dm_IC3(t1:t2),'DisplayName','<E_{hil}(t)> (6DM)','LineWidth',1.5);
ylabel('<E_{hil}(t)> [W m^{-2}]');
hold on
yyaxis right
plot(time(t1:t2),movmean(N2_bot_IC3(t1:t2),149),'DisplayName','N_{bot}^2(t) (6DM)','LineWidth',1.5);
% plot(time(t1:t2),E_hil_IC3_Nbot2_6dm(t1:t2),'DisplayName','<N_{bot}^2(t)> (6DM)','Color','black','LineWidth',1.5);
ylabel('N_{bot}^2 [s^{-2}]');
hold off
legend();
subplot(2,1,2)
plot(time(t1:t2),movmean(U_IC3(t1:t2),149),'DisplayName','U(t) (6DM)','LineWidth',1.5);
% plot(time(t1:t2),E_hil_IC3_U_6dm(t1:t2),'DisplayName','<U(t)> (6DM)','LineWidth',1.5);
ylabel('U [m s^{-1}]');
legend();

exportgraphics(ax3b,'figures/main/param/' + mooring(1) + '_EHil_Zoom.png');

%% M1: E_hil(t), comparative influence of N_bot^2 and U
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [5 2 40 22]);
set(0,'defaultAxesFontSize',18);

ax4 = figure;
subplot(3,2,1)
plot(time,Ehil6dm_M1,'DisplayName','<E_{hil}(t)> (6dm)','LineWidth',1.5);
legend();
xlabel('time');
ylabel('E_{hil} [W m^{-2}]');
title('<E_{hil}(t)> (6-day mean)');

subplot(3,2,3)
plot(time,E_hil_M1_Nbot2_6dm,'DisplayName','<N_{bot}^2(t)> (6dm)','Color','black','LineWidth',1.5);
legend();
xlabel('time');
ylabel('E_{hil} [W m^{-2}]');
title('<N_{bot}^2(t)> (6-day mean)');

subplot(3,2,5)
plot(time,E_hil_M1_U_6dm,'DisplayName','<U(t)> (6dm)','Color','green','LineWidth',1.5);
legend();
xlabel('time');
ylabel('E_{hil} [W m^{-2}]');
title('<U(t)> (6-day mean)');

subplot(3,2,2)
plot(time(t1:t2),Ehil6dm_M1(t1:t2),'DisplayName','<E_{hil}(t)> (6dm)','LineWidth',1.5);
legend();
xlabel('time');
ylabel('E_{hil} [W m^{-2}]');
title('<E_{hil}(t)> (6-day mean): close-up');

subplot(3,2,4)
plot(time(t1:t2),E_hil_M1_Nbot2_6dm(t1:t2),'DisplayName','<N_{bot}^2(t)> (6dm)','Color','black','LineWidth',1.5);
legend();
xlabel('time');
ylabel('E_{hil} [W m^{-2}]');
title('<N_{bot}^2(t)> (6-day mean): close-up');

subplot(3,2,6)
plot(time(t1:t2),E_hil_M1_U_6dm(t1:t2),'DisplayName','<U(t)> (6dm)','Color','green','LineWidth',1.5);
legend();
xlabel('time');
ylabel('E_{hil} [W m^{-2}]');
title('<U(t)> (6-day mean): close-up');

sgtitle('M1: Influence of N^2 and U on E_{hil}');
savefig('figures/main/param/' + mooring(2) + '_EHil-subplots');
exportgraphics(ax4,'figures/main/param/' + mooring(2) + '_EHil-subplots.png');


%% same as above but without zoom

%% M1: E_hil(t), comparative influence of N_bot^2 and U
set(0,'defaultAxesFontSize',14);

ax4a = figure;
subplot(3,1,1)
plot(time,Ehil6dm_M1,'DisplayName','<E_{hil}(t)> (6DM)','LineWidth',1.5);
legend();
ylabel({'Energy Available for';' Mixing E_{hil} [Wm^{-2}]'});
% title('<E_{hil}(t)>');

subplot(3,1,2)
plot(time,E_hil_M1_Nbot2_6dm,'DisplayName','<N_{bot}^2(t)> (6DM)','Color','black','LineWidth',1.5);
legend();
ylabel({'Energy Available for';' Mixing E_{hil} [Wm^{-2}]'});
% title('<N_{bot}^2(t)>');

subplot(3,1,3)
plot(time,E_hil_M1_U_6dm,'DisplayName','<U(t)> (6DM)','Color','green','LineWidth',1.5);
legend();
ylabel({'Energy Available for';' Mixing E_{hil} [Wm^{-2}]'});
% title('<U(t)>');


% sgtitle('M1: Influence of N^2 and U on E_{hil}');
savefig('figures/main/param/' + mooring(2) + '_EHil-subplotsNoZoom');
exportgraphics(ax4a,'figures/main/param/' + mooring(2) + '_EHil-subplotsNoZoom.png');

%% M1 Ehil: only ZOOM

set(0,'defaultAxesFontSize',15);

set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [5 2 40 15]);
ax4b = figure;
subplot(2,1,1)
yyaxis left
plot(time(t1:t2),Ehil6dm_M1(t1:t2),'DisplayName','<E_{hil}(t)> (6DM)','LineWidth',1.5);
ylabel('<E_{hil}(t)> [W m^{-2}]');
hold on
yyaxis right
plot(time(t1:t2),movmean(N2_bot_M1(t1:t2),149),'DisplayName','N_{bot}^2(t) (6DM)','LineWidth',1.5);
% plot(time(t1:t2),E_hil_M1_Nbot2_6dm(t1:t2),'DisplayName','<N_{bot}^2(t)> (6DM)','Color','black','LineWidth',1.5);
ylabel('N_{bot}^2 [s^{-2}]');
hold off
legend();
subplot(2,1,2)
plot(time(t1:t2),movmean(U_M1(t1:t2),149),'DisplayName','U(t) (6DM)','LineWidth',1.5);
% plot(time(t1:t2),E_hil_M1_U_6dm(t1:t2),'DisplayName','<U(t)> (6DM)','LineWidth',1.5);
ylabel('U [m s^{-1}]');
legend();

exportgraphics(ax4b,'figures/main/param/' + mooring(2) + '_EHil_Zoom.png');

%% Clear time segment
clear time1 time2;

%% Update fontsize again: back to big letters
set(0,'defaultAxesFontSize',20);
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [5 2 40 15]);

%% IC3: Histogram. Distribution of E_hil.

ax5 = figure;
histogram(Ehil6dm_IC3);
hold on
xline(mean(Ehil6dm_IC3),':','\mu','LineWidth',1.5,'FontSize',24,'LabelHorizontalAlignment','left');
xline(mean(Ehil6dm_IC3)+stdE_IC3,':','\mu + \sigma','LineWidth',1.5,'FontSize',24,'LabelHorizontalAlignment','left');
xline(mean(Ehil6dm_IC3)-stdE_IC3,':','\mu - \sigma','LineWidth',1.5,'FontSize',24,'LabelHorizontalAlignment','left');
hold off
xlabel('E_{hil}(t) [W m^{-2}]');
ylabel('No. of values');
% title('IC3: Distribution of E_{hil}(t)');

savefig('figures/main/param/' + mooring(1) + '_EHil-histogram');
exportgraphics(ax5,'figures/main/param/' + mooring(1) + '_EHil-histogram.png');

%% IC3: Histogram. Distribution of log10(E_hil(t))

ax5a = figure;
histogram(log10(Ehil6dm_IC3));
xlabel('log_{10}(E_{hil}(t)) [log W m^{-2}] ');
ylabel('No. of values');

exportgraphics(ax5a,'figures/main/param/' + mooring(1) + '_log10EHil-histogram.png');

%% M1: Histogram. Distribution of E_hil(t).

ax6 = figure;
histogram(Ehil6dm_M1);
hold on
xline(mean(Ehil6dm_M1),':','\mu','LineWidth',1.5,'FontSize',24,'LabelHorizontalAlignment','left');
xline(mean(Ehil6dm_M1)+stdE_M1,':','\mu + \sigma','LineWidth',1.5,'FontSize',24,'LabelHorizontalAlignment','left');
xline(mean(Ehil6dm_M1)-stdE_M1,':','\mu - \sigma','LineWidth',1.5,'FontSize',24,'LabelHorizontalAlignment','left');
hold off
xlabel('E_{hil}(t) [W m^{-2}]');
ylabel('No. of values');
% title('M1: Distribution of E_{hil}(t)');

savefig('figures/main/param/' + mooring(2) + '_EHil-histogram');
exportgraphics(ax6,'figures/main/param/' + mooring(2) + '_EHil-histogram.png');

%% M1: Histogram. Distribution of log10(E_hil(t))

ax6a = figure;
histogram(log10(Ehil6dm_M1));
xlabel('log_{10}(E_{hil}(t)) [log W m^{-2}] ');
ylabel('No. of values');

exportgraphics(ax6a,'figures/main/param/' + mooring(2) + '_log10EHil-histogram.png');

%% IC3: Correlation of normalised U and N2 with E

% 6DM correlations
[corr_E_N_IC3,pvalEN_IC3] = corr(E_hil_IC3_Nbot2_6dm,Ehil6dm_IC3);          % OK
[corr_E_U_IC3,pvalEU_IC3] = corr(E_hil_IC3_U_6dm,Ehil6dm_IC3);              % OK

[corr_E_N_IC3_KEN,pvalEN_IC3_KEN] = corr(E_hil_IC3_Nbot2_6dm,Ehil6dm_IC3,'type','Kendall');          % OK
[corr_E_U_IC3_KEN,pvalEU_IC3_KEN] = corr(E_hil_IC3_U_6dm,Ehil6dm_IC3,'type','Kendall');              % OK

% %% Correlations for hrly, mthly, multi-month means

% Hourly (original data)
[corr_E_N_IC3_hrly,pvalEN_IC3_hrly] = corr(E_hil_IC3_Nbot2,Ehil_IC3);
[corr_E_U_IC3_hrly,pvalEU_IC3_hrly] = corr(E_hil_IC3_U,Ehil_IC3);

%% M1: Correlation of <E_hil(t)> with <U(t)> and <N^2_bot(t)>

[corr_E_N_M1,pvalEN_M1] = corr(E_hil_M1_Nbot2_6dm,Ehil6dm_M1);
[corr_E_U_M1,pvalEU_M1] = corr(E_hil_M1_U_6dm,Ehil6dm_M1);

% Hourly
[corr_E_N_M1_hrly,pvalEN_M1_hrly] = corr(E_hil_M1_Nbot2,Ehil_M1);
[corr_E_U_M1_hrly,pvalEU_M1_hrly] = corr(E_hil_M1_U,Ehil_M1);

%% Save parameters for the next file

save Matfiles/E_hil.mat Ehil_IC3 Ehil_M1 Ehil6dm_IC3 Ehil6dm_M1;