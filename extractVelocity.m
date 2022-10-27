clear; close all; clc;
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [5 5 40 20]);
% Add subfolders to path
addpath(genpath(pwd));

%% Load Irminger Current Mooring Data
input = load('Data/merged_hourly_unfiltered_data_20142020.mat');

%% Initialise time parameters
time = datetime(input.time_grid,'ConvertFrom','datenum');   % time [h], datetime format
Ts = 3600;                                                  % Sampling time [s]  
fs = 2*pi/Ts;                                               % Sampling frequency [/s]
time_int = 0:Ts:length(input.time_grid)*Ts-1;               % Time [s], integer format
L = length(time_int);

Omega = 7.2921e-5;
fCor = 2*Omega*sin(input.lat(4)*pi/180);

%% FROM THIS POINT: only run the code for whichever mooring we want processed


%% IC3: extract levels in which there are observations

% Extract velocity from the CMs
u_CMs = cat(1,squeeze(input.U(4,:,72)),squeeze(input.U(4,:,73)),squeeze(input.U(4,:,74)),squeeze(input.U(4,:,99)),squeeze(input.U(4,:,100)),squeeze(input.U(4,:,102)),squeeze(input.U(4,:,145)),squeeze(input.U(4,:,155)),squeeze(input.U(4,:,156)));
v_CMs = cat(1,squeeze(input.V(4,:,72)),squeeze(input.V(4,:,73)),squeeze(input.V(4,:,74)),squeeze(input.V(4,:,99)),squeeze(input.V(4,:,100)),squeeze(input.V(4,:,102)),squeeze(input.V(4,:,145)),squeeze(input.V(4,:,155)),squeeze(input.V(4,:,156)));

u_CMs(isnan(u_CMs)) = 0;
v_CMs(isnan(v_CMs)) = 0;

% Average the CM velocities
% I need to ACTUALLY average these CM velocities
u_CMs = cat(1,u_CMs(1,:) + u_CMs(2,:) + u_CMs(3,:), u_CMs(4,:) + u_CMs(5,:) + u_CMs(6,:), u_CMs(7,:) + u_CMs(8,:) + u_CMs(9,:));
v_CMs = cat(1,v_CMs(1,:) + v_CMs(2,:) + v_CMs(3,:), v_CMs(4,:) + v_CMs(5,:) + v_CMs(6,:), v_CMs(7,:) + v_CMs(8,:) + v_CMs(9,:));

% Extract velocity from the ADCPs
u_adcp = squeeze(input.U(4,:,8:45))';
v_adcp = squeeze(input.V(4,:,8:45))';

% Create one array for u and v from the combination of CM and ADCP data
u = cat(1, u_adcp, u_CMs);
v = cat(1, v_adcp, v_CMs);

% save depths of the columns extracted above
dCM1 = (input.DGRID(72,4) + input.DGRID(73,4) + input.DGRID(74,4))/3;
dCM2 = (input.DGRID(99,4) + input.DGRID(100,4) + input.DGRID(102,4))/3;
dCM3 = (input.DGRID(145,4) + input.DGRID(155,4) + input.DGRID(156,4))/3;
depths = cat(1,input.DGRID(8:45,4),dCM1,dCM2,dCM3);

clear dCM1 dCM2 dCM3;
noOfMeters = length(depths);

%% Fill in missing values


% Fill missing values according to [double] the minimum window size needed
% to bridge the gap. This factor should smooth the value a little...
factor = 588;

u(u==0) = NaN;
v(v==0) = NaN;

for i = 1:length(depths)
    u_filled(i,:) = fillmissing(u(i,:),"movmean",2*factor);
    v_filled(i,:) = fillmissing(v(i,:),"movmean",2*factor);
end

% test to see whether all NaN values have been filled
u_filled(isnan(u_filled)) = 0;
v_filled(isnan(v_filled)) = 0;

figure
plot(time,v(30,:),'DisplayName','u','LineStyle','-');
hold on
plot(time,v_filled(30,:),'DisplayName','u_filled','LineStyle',':');
hold off
legend();

% 
u = u_filled;
v = v_filled;

%% Save the data
save Matfiles/IC3.mat u v depths time time_int noOfMeters;

%% IC2: input

IC2u = squeeze(input.U(3,:,:));
IC2v = squeeze(input.V(3,:,:));

%% IC2: Which levels actually have values

for i=1:length(IC2u(1,:))
    if ~isnan(IC2u(1,i))
        disp(i)
    end
end

%% IC2:ADCPs from 6-46

u_IC2_adcp = IC2u(:,6:46);
v_IC2_adcp = IC2v(:,6:46);

% %%
% figure
% plot(time,IC2u(:,47));

%% IC2: CM1 (62,71,72)

figure
% plot(time,IC2u(:,60:78));
% hold on
plot(time,IC2u(:,62));
hold on
plot(time,IC2u(:,71));
plot(time,IC2u(:,72));
hold off

u_IC2_CM1 = cat(1,IC2u(:,62) + IC2u(:,71) + IC2u(:,72));
v_IC2_CM1 = cat(1,IC2v(:,62) + IC2v(:,71) + IC2v(:,72));
dIC2_CM1 = (input.DGRID(62,3) + input.DGRID(71,3) + input.DGRID(72,3))/3;

%% IC2: CM2 (87,88,97,98)

figure
plot(time,IC2u(:,87));
hold on
plot(time,IC2u(:,88));
plot(time,IC2u(:,97));
plot(time,IC2u(:,98));
hold off

u_IC2_CM2 = cat(1,IC2u(:,87) + IC2u(:,88) + IC2u(:,97) + IC2u(:,98));
v_IC2_CM2 = cat(1,IC2v(:,87) + IC2v(:,88) + IC2v(:,97) + IC2v(:,98));
dIC2_CM2 = (input.DGRID(87,3) + input.DGRID(88,3) + input.DGRID(97,3) + input.DGRID(98,3))/4;

%% IC2: CM3 (134,135,138,139)

figure
plot(time,IC2u(:,134));
hold on
plot(time,IC2u(:,135));
plot(time,IC2u(:,138));
plot(time,IC2u(:,139));
hold off

u_IC2_CM3 = cat(1,IC2u(:,134) + IC2u(:,135) + IC2u(:,138) + IC2u(:,139));
v_IC2_CM3 = cat(1,IC2v(:,134) + IC2v(:,135) + IC2v(:,138) + IC2v(:,139));
dIC2_CM3 = (input.DGRID(134,3) + input.DGRID(135,3) + input.DGRID(138,3) + input.DGRID(139,3))/4;

%% IC2: CM4 (144,147,149,150)
% 14-15 missing

figure
plot(time,IC2u(:,144));
hold on
plot(time,IC2u(:,147));
plot(time,IC2u(:,149));
plot(time,IC2u(:,150));
hold off

u_IC2_CM4 = cat(1,IC2u(:,144) + IC2u(:,147) + IC2u(:,149) + IC2u(:,150));
v_IC2_CM4 = cat(1,IC2v(:,144) + IC2v(:,147) + IC2v(:,149) + IC2v(:,150));
dIC2_CM4 = (input.DGRID(144,3) + input.DGRID(147,3) + input.DGRID(149,3) + input.DGRID(150,3))/4;

%% IC2: CM5 (154,156,160)
figure
plot(time,IC2u(:,154));
hold on
plot(time,IC2u(:,156));
plot(time,IC2u(:,160));
hold off

u_IC2_CM5 = cat(1,IC2u(:,154) + IC2u(:,156) + IC2u(:,160));
v_IC2_CM5 = cat(1,IC2v(:,154) + IC2v(:,156) + IC2v(:,160));
dIC2_CM5 = (input.DGRID(154,3) + input.DGRID(156,3) + input.DGRID(160,3))/3;

%%

u_IC2 = cat(2, u_IC2_adcp, u_IC2_CM1, u_IC2_CM2, u_IC2_CM3, u_IC2_CM4, u_IC2_CM5);
v_IC2 = cat(2, v_IC2_adcp, v_IC2_CM1, v_IC2_CM2, v_IC2_CM3, v_IC2_CM4, v_IC2_CM5);
depths_IC2 = cat(1,input.DGRID(6:46,3),dIC2_CM1,dIC2_CM2,dIC2_CM3,dIC2_CM4,dIC2_CM5);
noOfMetersIC2 = length(depths_IC2);

%% Fill missing values: IC2

% Fill missing values according to [double] the minimum window size needed
% to bridge the gap. This factor should smooth the value a little...
factor = 588;

u_IC2(u_IC2==0) = NaN;
v_IC2(v_IC2==0) = NaN;

u_filled_IC2 = [];
v_filled_IC2 = [];

for i = 1:length(depths_IC2)
    u_filled_IC2(:,i) = fillmissing(u_IC2(:,i),"movmean",2*factor);
    v_filled_IC2(:,i) = fillmissing(v_IC2(:,i),"movmean",2*factor);
end

% test to see whether all NaN values have been filled
u_filled_IC2(isnan(u_filled_IC2)) = 0;
v_filled_IC2(isnan(v_filled_IC2)) = 0;

figure
plot(time,u_IC2(:,35),'DisplayName','u','LineStyle','-');
hold on
plot(time,u_filled_IC2(:,35),'DisplayName','u_filled','LineStyle',':');
hold off
legend();

% 
u_IC2 = u_filled_IC2;
v_IC2 = v_filled_IC2;

%% Save the data
save Matfiles/IC2.mat u_IC2 v_IC2 depths_IC2 time time_int noOfMetersIC2;


%% M1: input

M1u = squeeze(input.U(6,:,:));
M1v = squeeze(input.V(6,:,:));

%% M1: Which levels actually have values

valuesInLevel = [];

for i=1:length(M1u(1,:))
    if ~isnan(M1u(1,i))
        disp(i)
        X = i;
        valuesInLevel = [valuesInLevel X];
    end
end

%% M1: ADCP

figure
subplot(2,1,1)
plot(time,M1u(:,5:30));
subplot(2,1,2)
plot(time,M1v(:,5:30));

%% M1: Nortek Aquadopps at 700, 1200, 1430, 1645

% NA1: 700m

figure
subplot(2,2,1)
plot(time,M1u(:,72));
title('Aquadopp 700m: u');

% NA2: 1200m
subplot(2,2,2)
plot(time,M1u(:,123));
hold on
plot(time,M1u(:,122));
hold off
title('Aquadopp 1200m: u');

% NA3: 1430m
subplot(2,2,3)
plot(time,M1u(:,141));
title('Aquadopp 1430m: u');

% NA4: 1645m
subplot(2,2,4)
plot(time,M1u(:,156));
hold on
plot(time,M1u(:,157));
hold off
title('Aquadopp 1645m: u');


%% M1: Combine
M1u(isnan(M1u)) = 0;
M1v(isnan(M1v)) = 0;

u_m1 = cat(2,M1u(:,5:30),M1u(:,72),(M1u(:,123)+M1u(:,122)),...
    M1u(:,141),(M1u(:,156)+M1u(:,157)));

v_m1 = cat(2,M1v(:,5:30),M1v(:,72),(M1v(:,123)+M1v(:,122)),...
    M1v(:,141),(M1v(:,156)+M1v(:,157)));

depths_M1= cat(1,input.DGRID(5:30,6),...
    input.DGRID(72,6),...
    (4*input.DGRID(123,6) + 2*input.DGRID(122,6))/6,...
    input.DGRID(141,6),...
    (4*input.DGRID(156,6) + 2*input.DGRID(157,6))/6);
noOfMetersM1 = length(depths_M1);

%% M1: fill missing values

u_filled_M1 = fillMissingVelocity(noOfMetersM1,u_m1,time,360);
v_filled_M1 = fillMissingVelocity(noOfMetersM1,v_m1,time,360);

u_M1 = u_filled_M1;
v_M1 = v_filled_M1;

save Matfiles/M1.mat u_M1 v_M1 depths_M1 noOfMetersM1 time time_int;

% Fill missing values according to [double] the minimum window size needed
% to bridge the gap. This factor should smooth the value a little...
% factor = 588;
% 
% u_M1(u_M1==0) = NaN;
% v_M1(v_M1==0) = NaN;
% 
% u_filled_M1 = [];
% v_filled_M1 = [];
% 
% for i = 1:length(depths_IC2)
%     u_filled_M1(:,i) = fillmissing(u_IC2(:,i),"movmean",2*factor);
%     v_filled_IC2(:,i) = fillmissing(v_IC2(:,i),"movmean",2*factor);
% end

% % test to see whether all NaN values have been filled
% u_filled_IC2(isnan(u_filled_IC2)) = 0;
% v_filled_IC2(isnan(v_filled_IC2)) = 0;
% 
% figure
% plot(time,u_IC2(:,35),'DisplayName','u','LineStyle','-');
% hold on
% plot(time,u_filled_IC2(:,35),'DisplayName','u_filled','LineStyle',':');
% hold off
% legend();
% 
% % 
% u_IC2 = u_filled_IC2;
% v_IC2 = v_filled_IC2;
