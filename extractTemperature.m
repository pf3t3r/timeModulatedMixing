% Here, the following is done:
% - observational data from Irminger Sea is loaded
% - temperature (T) is extracted for EACH MOORING
% - conservative temperature (Theta) is calculated and output into a file.
% If these files already exist this method can be skipped.
% - a plot is made of T vs Theta. This is just for illustrative purposes.

% Add subfolders to path
addpath(genpath(pwd));

clear; clc; close all;
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [5 5 40 20]);

%% Load data
input = load('Data/merged_hourly_unfiltered_data_20142020.mat');
time = datetime(input.time_grid,'ConvertFrom','datenum');   % time [h], datetime format

latitude = input.lat;
lon = input.lon';
dgrid = input.DGRID;

T = input.T;

% For each of the five moorings IC1, IC2, IC3, IC4, and M1, we extract the
% levels at which data is recorded, merge data that was recorded by a
% single sensor to create a full six-year time series, and find the
% (weighted) average depth that corresponds to each sensor time series.

%% Extract IC1 Temperature

% IC1: 7 sensors at the following depth levels
% 1: 3,7,9,15
% 2: 15,20,32
% 3: 43,45,46,48
% 4: 73,74,76
% 5: 98,100,101
% 6: 134,136
% 7: 161

% NOTE that sensors 4, 5, and 7 are each missing two years of T data
% Sensor 4 is missing 18-20
% Sensor 5 is missing 16-18
% Sensor 6 is missing 18-20
% Sensor 7 is missing 16-18

T_IC1 = squeeze(T(2,:,:));

IC1_15a = nan(length(T_IC1(:,15)),1);
IC1_15b = nan(length(T_IC1(:,15)),1);
IC1_15a(1:8500) = T_IC1(1:8500,15);
IC1_15b(8501:end) = T_IC1(8501:end,15);

for i=1:length(time)
    T_IC1_1(i) = sum([T_IC1(i,3),T_IC1(i,7),T_IC1(i,9),IC1_15a(i)],'omitnan');
    T_IC1_2(i) = sum([IC1_15b(i),T_IC1(i,20),T_IC1(i,32)],'omitnan');
    T_IC1_3(i) = sum([T_IC1(i,43),T_IC1(i,45),T_IC1(i,46),T_IC1(i,48)],'omitnan');
    T_IC1_4(i) = sum([T_IC1(i,73),T_IC1(i,74),T_IC1(i,76)],'omitnan');
    T_IC1_5(i) = sum([T_IC1(i,98),T_IC1(i,100),T_IC1(i,101)],'omitnan');
    T_IC1_6(i) = sum([T_IC1(i,134),T_IC1(i,136)],'omitnan');
    T_IC1_7(i) = sum([T_IC1(i,161)],'omitnan');
end

T_IC1_avg = cat(1,T_IC1_1,T_IC1_2,T_IC1_3,T_IC1_4,T_IC1_5,T_IC1_6,T_IC1_7);
T_IC1_avg(T_IC1_avg==0) = NaN;

% Weighted average depth
depths_T_IC1 = [(dgrid(3,2) + 2*dgrid(7,2) + 2*dgrid(9,2) + dgrid(15,2))/6 ...
    (dgrid(15,2) + 4*dgrid(20,2) + dgrid(32,2))/6 ...
    (2*dgrid(43,2) + 2*dgrid(45,2) + dgrid(46,2) + dgrid(48,2))/4 ...
    (dgrid(73,2) + 2*dgrid(74,2) + dgrid(76,2))/4 ...
    (dgrid(98,2) + dgrid(100,2) + 2*dgrid(101,2))/4 ...
    (2*dgrid(134,2) + 2*dgrid(136,2))/4 ...
    dgrid(161,2)];

clear T_IC1_1 T_IC1_2 T_IC1_3 T_IC1_4 T_IC1_5 T_IC1_6 T_IC1_7;

%% Extract IC2 Temperature

% IC2: 7 sensors at the following depths
% 1: 3,11
% 2: 17,18,28,36
% 3: 39,46,47,50
% 4: 62,63,71,73
% 5: 88,98
% 4: 134,135,139
% 5: 156,160

% NOTE that sensor 1 is missing data from 14-15 and 16-18. Additionally,
% it's possible that I have some levels mixed up between sensor 2 and 3.

T_IC2 = squeeze(T(3,:,:));

IC1_62a = nan(length(T_IC2(:,62)),1);
IC1_62b = nan(length(T_IC2(:,62)),1);
IC1_62a(1:35200) = T_IC2(1:35200,62);
IC1_62b(35201:end) = T_IC2(35201:end,62);

for i=1:length(time)
    T_IC2_1(i) = sum([T_IC1(i,3),T_IC1(i,11)],'omitnan');
    T_IC2_2(i) = sum([T_IC1(i,17),T_IC1(i,18),T_IC1(i,28),T_IC1(i,36)],'omitnan');
    T_IC2_3(i) = sum([T_IC1(i,39),T_IC1(i,46),T_IC1(i,47),T_IC1(i,50)],'omitnan');
    T_IC2_4(i) = sum([IC1_62b(i),T_IC1(i,63),T_IC1(i,71),T_IC1(i,73)],'omitnan');
    T_IC2_5(i) = sum([T_IC1(i,88),T_IC1(i,98)],'omitnan');
    T_IC2_6(i) = sum([T_IC1(i,134),T_IC1(i,135),T_IC1(i,139)],'omitnan');
    T_IC2_7(i) = sum([T_IC1(i,156),T_IC1(i,160)],'omitnan');
end

T_IC2_avg = cat(1,T_IC2_1,T_IC2_2,T_IC2_3,T_IC2_4,T_IC2_5,T_IC2_6,T_IC2_7);
T_IC2_avg(T_IC2_avg==0) = NaN;

% Weighted average depth
depths_T_IC2 = [(dgrid(3,3) + 2*dgrid(11,3))/3 ...
    (dgrid(17,3) + dgrid(18,3) + 2*dgrid(28,3) + 2*dgrid(36,3))/6 ...
    (2*dgrid(39,3) + dgrid(46,3) + dgrid(47,3) + 2*dgrid(50,3))/6 ...
    (2*dgrid(62,3) + 2*dgrid(63,3) + dgrid(71,3) + dgrid(73,3))/6 ...
    (4*dgrid(88,3) + 2*dgrid(98,3))/6 ...
    (2*dgrid(134,3) + 2*dgrid(135,3) + 2*dgrid(139,3))/6 ...
    (4*dgrid(156,3) + 2*dgrid(160,3))/6];

clear T_IC2_1 T_IC2_2 T_IC2_3 T_IC2_4 T_IC2_5 T_IC2_6 T_IC2_7;

%% Extract IC3 Temperature

% IC3: four sensors at the following depth levels
% 1: 5,9
% 2: 18,19,23
% 3: 47,48,49
% 4: 72,73,74
% 5: 99,100
% 6: 147,155,156

% NOTE that sensor 1 is missing data from 14-15

T_IC3 = squeeze(T(4,:,:));

for i=1:length(time)
    T_IC3_1(i) = sum([T_IC3(i,5),T_IC3(i,9)],'omitnan');
    T_IC3_2(i) = sum([T_IC3(i,18),T_IC3(i,19),T_IC3(i,23)],'omitnan');
    T_IC3_3(i) = sum([T_IC3(i,47),T_IC3(i,48),T_IC3(i,49)],'omitnan');
    T_IC3_4(i) = sum([T_IC3(i,72),T_IC3(i,73),T_IC3(i,74)],'omitnan');
    T_IC3_5(i) = sum([T_IC3(i,99),T_IC3(i,100)],'omitnan');
    T_IC3_6(i) = sum([T_IC3(i,147),T_IC3(i,155),T_IC3(i,156)],'omitnan');
end

T_IC3_avg = cat(1,T_IC3_1,T_IC3_2,T_IC3_3,T_IC3_4,T_IC3_5,T_IC3_6);
T_IC3_avg(T_IC3_avg==0) = NaN;

% Weighted average depth
depths_T_IC3 = [(dgrid(5,4) + 4*dgrid(9,4))/5 ...
    (dgrid(18,4) + dgrid(19,4) + 4*dgrid(23,4))/6 ...
    (dgrid(47,4) + dgrid(48,4) + 4*dgrid(49,4))/6 ...
    (dgrid(72,4) + dgrid(73,4) + 4*dgrid(74,4))/6 ...
    (2*dgrid(99,4) + 4*dgrid(100,4))/6 ...
    (2*dgrid(147,4) + dgrid(155,4) + 3*dgrid(156,4))/6];

clear T_IC3_1 T_IC3_2 T_IC3_3 T_IC3_4 T_IC3_5 T_IC3_6;

%% Extract IC4 Temperature

% IC4: six sensors at the following depth levels
% 1: 5,8,11
% 2: 18,20,21,31
% 3: 48,49
% 4: 74
% 5: 100
% 6: 136,137,154

% NOTE that sensor 1 is missing data from 14-15

T_IC4 = squeeze(T(5,:,:));

for i=1:length(time)
    T_IC4_1(i) = sum([T_IC4(i,5),T_IC4(i,8),T_IC4(i,11)],'omitnan');
    T_IC4_2(i) = sum([T_IC4(i,18),T_IC4(i,20),T_IC4(i,21),T_IC4(i,31)],'omitnan');
    T_IC4_3(i) = sum([T_IC4(i,48),T_IC4(i,49)],'omitnan');
    T_IC4_4(i) = sum([T_IC4(i,74)],'omitnan');
    T_IC4_5(i) = sum([T_IC4(i,100)],'omitnan');
    T_IC4_6(i) = sum([T_IC4(i,136),T_IC4(i,137),T_IC4(i,154)],'omitnan');
end

T_IC4_avg = cat(1,T_IC4_1,T_IC4_2,T_IC4_3,T_IC4_4,T_IC4_5,T_IC4_6);
T_IC4_avg(T_IC4_avg==0) = NaN;

% Weighted average depth
depths_T_IC4 = [(dgrid(5,5) + 2*dgrid(8,5) + 2*dgrid(11,5))/5 ...
    (dgrid(18,5) + dgrid(20,5) + 2*dgrid(21,5) + 2*dgrid(31,5))/6 ...
    (2*dgrid(48,5) + 4*dgrid(49,5))/6 ...
    dgrid(74,5) ...
    dgrid(100,5) ...
    (dgrid(136,5) + 3*dgrid(137,5) + 2*dgrid(154,5))/6];

clear T_IC4_1 T_IC4_2 T_IC4_3 T_IC4_4 T_IC4_5 T_IC4_6;

%% Extract M1 Temperature

% M1: ten sensors at the following depth levels
% 1: 6
% 2: 11,12
% 3: 21,22
% 4: 36,37
% 5: 51,52
% 6: 72
% 7: 92
% 8: 122,123
% 9: 141
% 10: 156,157

T_M1 = squeeze(T(6,:,:));

for i=1:length(time)
    T_M1_1(i) = sum([T_M1(i,6)],'omitnan');
    T_M1_2(i) = sum([T_M1(i,11),T_M1(i,12)],'omitnan');
    T_M1_3(i) = sum([T_M1(i,21),T_M1(i,22)],'omitnan');
    T_M1_4(i) = sum([T_M1(i,36),T_M1(i,37)],'omitnan');
    T_M1_5(i) = sum([T_M1(i,51),T_M1(i,52)],'omitnan');
    T_M1_6(i) = sum([T_M1(i,72)],'omitnan');
    T_M1_7(i) = sum([T_M1(i,92)],'omitnan');
    T_M1_8(i) = sum([T_M1(i,122),T_M1(i,123)],'omitnan');
    T_M1_9(i) = sum([T_M1(i,141)],'omitnan');
    T_M1_10(i) = sum([T_M1(i,156),T_M1(i,157)],'omitnan');
end

T_M1_avg = cat(1,T_M1_1,T_M1_2,T_M1_3,T_M1_4,T_M1_5,T_M1_6,T_M1_7,T_M1_8,T_M1_9,T_M1_10);
T_M1_avg(T_M1_avg==0) = NaN;

% Weighted average depth
depths_T_M1 = [dgrid(6,6) ...
    (3*dgrid(11,6) + 3*dgrid(12,6))/6 ...
    (2*dgrid(21,6) + 4*dgrid(22,6))/6 ...
    (2*dgrid(36,6) + 4*dgrid(37,6))/6 ...
    (2*dgrid(51,6) + 4*dgrid(52,6))/6 ...
    dgrid(72,6) ...
    dgrid(92,6) ...
    (2*dgrid(122,6) + 4*dgrid(123,6))/6 ...
    dgrid(141,6) ...
    (2*dgrid(156,6) + 4*dgrid(157,6))/6];

clear T_M1_1 T_M1_2 T_M1_3 T_M1_4 T_M1_5 T_M1_6 T_M1_7 T_M1_8 T_M1_9 T_M1_10;

%% find the levels

% figure
% plot(time,T_M1(:,155:161),'LineWidth',2);
% legend();

%% Save Temperature from the five moorings

T_IC1 = T_IC1_avg;
T_IC2 = T_IC2_avg;
T_IC3 = T_IC3_avg;
T_IC4 = T_IC4_avg;
T_M1 = T_M1_avg;

clear T_IC1_avg T_IC2_avg T_IC3_avg T_IC4_avg T_M1_avg;

save Matfiles/T.mat T_IC1 T_IC2 T_IC3 T_IC4 T_M1;

%% Temperature: fill in missing values

% There is a 2-year gap in 3/7 of levels in this time series
T_IC1_f = fillMissingTemperature(7,T_IC1',time,18000);

% There is a 3-year gap in 1/7 of levels in this time series
% actually this is BROKEN
% T_IC2_f = fillMissingTemperature(7,T_IC2',time,35000);

% There is a 1-year gap for 1/6 of the levels in this time series
T_IC3_f = fillMissingTemperature(6,T_IC3',time,9000);

% There is a 1-year gap for 1/6 of the levels in this time series
T_IC4_f = fillMissingTemperature(6,T_IC4',time,9000);

% There are no major gaps in this time series
T_M1_f = fillMissingTemperature(10,T_M1',time);

save Matfiles/T_filled.mat T_IC1_f T_IC3_f T_IC4_f T_M1_f;

%% FIG: IC1-T

ax1 = figure;
hold on
for i=1:length(depths_T_IC1)
    plot(time,T_IC1_f(i,:));
end
lgd = legend({num2str(depths_T_IC1(1)),num2str(depths_T_IC1(2)),...
    num2str(depths_T_IC1(3)),num2str(depths_T_IC1(4)),...
    num2str(depths_T_IC1(5)),num2str(depths_T_IC1(6)),...
    num2str(depths_T_IC1(7))});
title(lgd,'Sensor Depth (m)');
title('IC1: T(z,t)');
xlabel('Time');
ylabel('Temperature (C)');

mooring = "IC1";
savefig('figures/main/extract/' + mooring +'_temperature');
exportgraphics(ax1,'figures/main/extract/' + mooring + '_temperature.png');

%% FIG: IC3 T

ax2 = figure;
hold on
for i=1:length(depths_T_IC3)
    plot(time,T_IC3_f(i,:));
end
lgd = legend({num2str(depths_T_IC3(1)),num2str(depths_T_IC3(2)),...
    num2str(depths_T_IC3(3)),num2str(depths_T_IC3(4)),...
    num2str(depths_T_IC3(5)),num2str(depths_T_IC3(6)),...
    });
title(lgd,'Sensor Depth (m)');
title('IC3: T(z,t)');
xlabel('Time');
ylabel('Temperature (C)');

mooring = "IC3";
savefig('figures/main/extract/' + mooring +'_temperature');
exportgraphics(ax2,'figures/main/extract/' + mooring + '_temperature.png');

%% FIG: IC4: T

ax3 = figure;
hold on
for i=1:length(depths_T_IC4)
    plot(time,T_IC4_f(i,:));
end
lgd = legend({num2str(depths_T_IC4(1)),num2str(depths_T_IC4(2)),...
    num2str(depths_T_IC4(3)),num2str(depths_T_IC4(4)),...
    num2str(depths_T_IC4(5)),num2str(depths_T_IC4(6)),...
    });
title(lgd,'Sensor Depth (m)');
title('IC4: T(z,t)');
xlabel('Time');
ylabel('Temperature (C)');

mooring = "IC4";
savefig('figures/main/extract/' + mooring +'_temperature');
exportgraphics(ax3,'figures/main/extract/' + mooring + '_temperature.png');

%% FIG: M1 T

ax4 = figure;
hold on
for i=1:length(depths_T_M1)
    plot(time,T_M1_f(i,:));
end
lgd = legend({num2str(depths_T_M1(1)),num2str(depths_T_M1(2)),...
    num2str(depths_T_M1(3)),num2str(depths_T_M1(4)),...
    num2str(depths_T_M1(5)),num2str(depths_T_M1(6)),...
    num2str(depths_T_M1(7)),num2str(depths_T_M1(8)),...
    num2str(depths_T_M1(9)),num2str(depths_T_M1(10)),...
    });
title(lgd,'Sensor Depth (m)');
title('M1: T(z,t)');
xlabel('Time');
ylabel('Temperature (C)');

mooring = "M1";
savefig('figures/main/extract/' + mooring +'_temperature');
exportgraphics(ax4,'figures/main/extract/' + mooring + '_temperature.png');

%% Depths: save for later

save Matfiles/depthsT depths_T_IC1 depths_T_IC2 depths_T_IC3 depths_T_IC4 depths_T_M1;

%% Pressures: needed for T -> CT conversion

p_IC1 = sw_pres(depths_T_IC1,latitude(2))';
p_IC2 = sw_pres(depths_T_IC2,latitude(3))';
p_IC3 = sw_pres(depths_T_IC3,latitude(4))';
p_IC4 = sw_pres(depths_T_IC4,latitude(5))';
p_M1 = sw_pres(depths_T_M1,latitude(6))';

save Matfiles/p_T.mat p_IC1 p_IC2 p_IC3 p_IC4 p_M1;
%% Convert Temperature (T) into Conservative Temperature (CT/Theta)


%% T1->CT1
if isfile('CT_IC1.mat')
    disp('CT already calculated for IC1');
    load CT_IC1.mat;
else
    disp('Calculating CT for IC1');
    load('SA_IC1.mat','absoluteSalinity');
    [CT_IC1] = convertIntoConservativeTemperature(absoluteSalinity,T_IC1_f([1,3,5,6,7],:),p_IC1([1,3,5,6,7]),time,"Matfiles/CT_IC1");
end

%% T2->CT2
% SEE PROBLEM WITH T_IC2_f above: won't work...
% 
% if isfile('CT_IC2.mat')
%     disp('CT already calculated for IC2');
%     load CT_IC2.mat;
% else
%     disp('Calculating CT for IC2');
%     load SA_IC2.mat;
%     [CT_IC2] = convertIntoConservativeTemperature(absoluteSalinity,T_IC2_f([1,3,5,6,7],:),p_IC2([1,3,5,6,7]),time,"CT_IC2");
% end

%% T3->CT3
if isfile('CT_IC3.mat')
    disp('CT already calculated for IC3');
    load CT_IC3.mat;
else
    disp('Calculating CT for IC3');
    load('SA_IC3.mat','absoluteSalinity');
    [CT_IC3] = convertIntoConservativeTemperature(absoluteSalinity,T_IC3_f([1,3,5,6],:),p_IC3([1,3,5,6]),time,"Matfiles/CT_IC3");
end

%% T4->CT4
if isfile('CT_IC4.mat')
    disp('CT already calculated for IC4');
    load CT_IC4.mat;
else
    disp('Calculating CT for IC4');
    load('SA_IC4.mat','absoluteSalinity');
    [CT_IC4] = convertIntoConservativeTemperature(absoluteSalinity,T_IC4_f([1,3,5,6],:),p_IC4([1,3,5,6]),time,"Matfiles/CT_IC4");
end

%% T-M1->CT-M1
if isfile('CT_M1.mat')
    disp('CT already calculated for M1');
    load CT_M1.mat;
else
    disp('Calculating CT for M1');
    load('SA_M1.mat','absoluteSalinity');
    [CT_M1] = convertIntoConservativeTemperature(absoluteSalinity,T_M1_f,p_M1,time,"Matfiles/CT_M1");
end