clc; close all; clear;
load('Matfiles/IC3.mat','time');
load('Matfiles/SA_IC3.mat');
load('Matfiles/CT_IC3.mat');
depth_SA = load('Matfiles/depthsS.mat','depths_IC3').depths_IC3;
depth_CT = load('Matfiles/depthsT.mat','depths_T_IC3').depths_T_IC3([1,3,5,6]);
t1 = 33000;
t2 = 37500;

figure;
subplot(2,2,1);plot(time(t1:t2),absoluteSalinity(1,t1:t2));title(depth_SA(1));
subplot(2,2,2);plot(time(t1:t2),absoluteSalinity(2,t1:t2));title(depth_SA(2));
subplot(2,2,3);plot(time(t1:t2),absoluteSalinity(3,t1:t2));title(depth_SA(3));
subplot(2,2,4);plot(time(t1:t2),absoluteSalinity(4,t1:t2));title(depth_SA(4));

% SUBPLOT 1
% before July 2019 = level 5; after July 2018 = level 19

% SUBPLOT 2
% before July 2018 = level 48; after July 2018 = level 49;

% SUBPLOT 3
% before July 2018 = level 100, after July 2018 = level 99

% SUBPLOT 4
% before July 2018 = level 147, after July 2018 = level 156

figure;
subplot(2,2,1);plot(time(t1:t2),conservativeTemperature(1,t1:t2));title(depth_CT(1));
subplot(2,2,2);plot(time(t1:t2),conservativeTemperature(2,t1:t2));title(depth_CT(2));
subplot(2,2,3);plot(time(t1:t2),conservativeTemperature(3,t1:t2));title(depth_CT(3));
subplot(2,2,4);plot(time(t1:t2),conservativeTemperature(4,t1:t2));title(depth_CT(4));

% IC3: four sensors at the following depth levels
% 1: 5,9 --->> no change July 2018
% 3: 47,48,49 ---> no change July 2018
% 5: 99,100 ---> 100/99
% 6: 147,155,156 ---> 147/156

T = load('Data/merged_hourly_unfiltered_data_20142020.mat').T;
%%

figure;
plot(time,squeeze(T(4,:,47)));
hold on
plot(time,squeeze(T(4,:,48)));
plot(time,squeeze(T(4,:,49)));
hold off