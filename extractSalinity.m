% Here, the following is done:
% - observational data from Irminger Sea is loaded
% - practical salinity (PS) is extracted for EACH MOORING
% - absolute salinity (SA) is calculated and output into a file. If these
% files already exist this method can be skipped.
% - a plot is made of SP vs SA. This is just for illustrative purposes.

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

SP = input.S;

% For each of the five moorings IC1, IC2, IC3, IC4, and M1, we extract the
% levels at which data is recorded, merge data that was recorded by a
% single sensor to create a full six-year time series, and find the
% (weighted) average depth that corresponds to each sensor time series.

%% Extract IC1 Salinity

% IC1: 5 sensors at the following depth levels
% 1: 3,7,9,15
% 2: 43,45,46,48
% 3: 98,100,101
% 4: 134,136
% 5: 161

% NOTE that sensors 3, 4, and 5 are each missing two years of Sp data
% Sensor 3 is missing 16-18
% Sensor 4 is missing 18-20
% Sensor 5 is missing 16-18

SP_IC1 = squeeze(SP(2,:,:));
for i=1:length(time)
    SP_IC1_1(i) = sum([SP_IC1(i,3),SP_IC1(i,7),SP_IC1(i,9),SP_IC1(i,15)],'omitnan');
    SP_IC1_2(i) = sum([SP_IC1(i,43),SP_IC1(i,45),SP_IC1(i,46),SP_IC1(i,48)],'omitnan');
    SP_IC1_3(i) = sum([SP_IC1(i,98),SP_IC1(i,100),SP_IC1(i,101)],'omitnan');
    SP_IC1_4(i) = sum([SP_IC1(i,136)],'omitnan');
    SP_IC1_5(i) = sum([SP_IC1(i,161)],'omitnan');
end

SP_IC1_avg = cat(1,SP_IC1_1,SP_IC1_2,SP_IC1_3,SP_IC1_4,SP_IC1_5);
SP_IC1_avg(SP_IC1_avg==0) = NaN;

% Weighted average depth
depths_IC1 = [(dgrid(3,2) + 2*dgrid(7,2) + 2*dgrid(9,2) + dgrid(15,2))/6 ...
    (2*dgrid(43,2) + 2*dgrid(45,2) + dgrid(46,2) + dgrid(48,2))/6 ...
    (dgrid(98,2) + dgrid(100,2) + 2*dgrid(101,2))/4 ...
    (dgrid(136,2)) ...
    (dgrid(161,2))];

clear SP_IC1_1 SP_IC1_2 SP_IC1_3 SP_IC1_4 SP_IC1_5;

%% Extract IC2 Salinity

% IC2: 5 sensors at the following depths
% 1: 3,11,18
% 2: 36,39,45,46
% 3: 88,98
% 4: 134,135,139
% 5: 156,160

% NOTE that sensor 1 is missing data from 16-18. It is possible that this
% data is at level 62...

SP_IC2 = squeeze(SP(3,:,:));
for i=1:length(time)
    SP_IC2_1(i) = sum([SP_IC2(i,3),SP_IC2(i,11),SP_IC2(i,18)],'omitnan');
    SP_IC2_2(i) = sum([SP_IC2(i,36),SP_IC2(i,39),SP_IC2(i,45),SP_IC2(i,46)],'omitnan');
    SP_IC2_3(i) = sum([SP_IC2(i,88),SP_IC2(i,98)],'omitnan');
    SP_IC2_4(i) = sum([SP_IC2(i,134),SP_IC2(i,135),SP_IC2(i,139)],'omitnan');
    SP_IC2_5(i) = sum([SP_IC2(i,156),SP_IC2(i,160)],'omitnan');
end

SP_IC2_avg = cat(1,SP_IC2_1,SP_IC2_2,SP_IC2_3,SP_IC2_4,SP_IC2_5);
SP_IC2_avg(SP_IC2_avg==0) = NaN;

% Weighted average depth
depths_IC2 = [(dgrid(3,3) + 2*dgrid(11,3) + dgrid(18,3))/4 ...
    (2*dgrid(36,3) + 2*dgrid(39,3) + dgrid(45,3) + dgrid(46,3))/6 ...
    (4*dgrid(88,3) + 2*dgrid(98,3))/6 ...
    (2*dgrid(134,3) + 2*dgrid(135,3) + 2*dgrid(139,3))/6 ...
    (4*dgrid(156,3) + 2*dgrid(160,3))/6];

clear SP_IC2_1 SP_IC2_2 SP_IC2_3 SP_IC2_4 SP_IC2_5;

%% Extract IC3 Salinity

% IC3: four sensors at the following depth levels
% 1: 5,9,19
% 2: 47,48,49
% 3: 99,100
% 4: 147,155,156

SP_IC3 = squeeze(SP(4,:,:));
for i=1:length(time)
    SP_IC3_1(i) = sum([SP_IC3(i,5),SP_IC3(i,9),SP_IC3(i,19)],'omitnan');
    SP_IC3_2(i) = sum([SP_IC3(i,47),SP_IC3(i,48),SP_IC3(i,49)],'omitnan');
    SP_IC3_3(i) = sum([SP_IC3(i,99),SP_IC3(i,100)],'omitnan');
    SP_IC3_4(i) = sum([SP_IC3(i,147),SP_IC3(i,155),SP_IC3(i,156)],'omitnan');
end

SP_IC3_avg = cat(1,SP_IC3_1,SP_IC3_2,SP_IC3_3,SP_IC3_4);
SP_IC3_avg(SP_IC3_avg==0) = NaN;

% Weighted average depth
depths_IC3 = [(dgrid(5,4) + 4*dgrid(9,4) + dgrid(19,4))/6 ...
    (dgrid(47,4) + dgrid(48,4) + 4*dgrid(49,4))/6 ...
    (2*dgrid(99,4) + 4*dgrid(100,4))/6 ...
    (2*dgrid(147,4) + dgrid(155,4) + 3*dgrid(156,4))/6];

clear SP_IC3_1 SP_IC3_2 SP_IC3_3 SP_IC3_4;

%% Extract IC4 Salinity

% IC4: four sensors at the following depth levels
% 1: 5,8,11,20
% 2: 48,49
% 3: 100
% 4: 136,137,154

SP_IC4 = squeeze(SP(5,:,:));
for i=1:length(time)
    SP_IC4_1(i) = sum([SP_IC4(i,5),SP_IC4(i,8),SP_IC4(i,11),SP_IC4(i,20)],'omitnan');
    SP_IC4_2(i) = sum([SP_IC4(i,48),SP_IC4(i,49)],'omitnan');
    SP_IC4_3(i) = sum([SP_IC4(i,100)],'omitnan');
    SP_IC4_4(i) = sum([SP_IC4(i,136),SP_IC4(i,137),SP_IC4(i,154)],'omitnan');
end

SP_IC4_avg = cat(1,SP_IC4_1,SP_IC4_2,SP_IC4_3,SP_IC4_4);
SP_IC4_avg(SP_IC4_avg==0) = NaN;

% Weighted average depth
depths_IC4 = [(dgrid(5,5) + 2*dgrid(8,5) + 2*dgrid(11,2) + dgrid(20,2))/6 ...
    (2*dgrid(48,5) + 4*dgrid(49,5))/6 ...
    (6*dgrid(100,5))/6 ...V
    (dgrid(136,5) + 3*dgrid(137,5) + 2*dgrid(154,5))/6];

clear SP_IC4_1 SP_IC4_2 SP_IC4_3 SP_IC4_4;

%% Extract M1 Salinity

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

SP_M1 = squeeze(SP(6,:,:));
for i=1:length(time)
    SP_M1_1(i) = sum([SP_M1(i,6)],'omitnan');
    SP_M1_2(i) = sum([SP_M1(i,11),SP_M1(i,12)],'omitnan');
    SP_M1_3(i) = sum([SP_M1(i,21),SP_M1(i,22)],'omitnan');
    SP_M1_4(i) = sum([SP_M1(i,36),SP_M1(i,37)],'omitnan');
    SP_M1_5(i) = sum([SP_M1(i,51),SP_M1(i,52)],'omitnan');
    SP_M1_6(i) = sum([SP_M1(i,72)],'omitnan');
    SP_M1_7(i) = sum([SP_M1(i,92)],'omitnan');
    SP_M1_8(i) = sum([SP_M1(i,122),SP_M1(i,123)],'omitnan');
    SP_M1_9(i) = sum([SP_M1(i,141)],'omitnan');
    SP_M1_10(i) = sum([SP_M1(i,156),SP_M1(i,157)],'omitnan');
end

SP_M1_avg = cat(1,SP_M1_1,SP_M1_2,SP_M1_3,SP_M1_4,SP_M1_5,SP_M1_6,SP_M1_7,SP_M1_8,SP_M1_9,SP_M1_10);
SP_M1_avg(SP_M1_avg==0) = NaN;

% Weighted average depth
depths_M1 = [dgrid(6,6) ...
    (3*dgrid(11,6) + 3*dgrid(12,6))/6 ...
    (2*dgrid(21,6) + 4*dgrid(22,6))/6 ...
    (2*dgrid(36,6) + 4*dgrid(37,6))/6 ...
    (2*dgrid(51,6) + 4*dgrid(52,6))/6 ...
    dgrid(72,6) ...
    dgrid(92,6) ...
    (2*dgrid(122,6) + 4*dgrid(123,6))/6 ...
    dgrid(141,6) ...
    (2*dgrid(156,6) + 4*dgrid(157,6))/6];

clear SP_M1_1 SP_M1_2 SP_M1_3 SP_M1_4 SP_M1_5 SP_M1_6 SP_M1_7 SP_M1_8 SP_M1_9 SP_M1_10;

%% Save Practical Salinity from the five moorings

SP_IC1 = SP_IC1_avg;
SP_IC2 = SP_IC2_avg;
SP_IC3 = SP_IC3_avg;
SP_IC4 = SP_IC4_avg;
SP_M1 = SP_M1_avg;

clear SP_IC1_avg SP_IC2_avg SP_IC3_avg SP_IC4_avg SP_M1_avg;

save Matfiles/Sp.mat SP_IC1 SP_IC2 SP_IC3 SP_IC4 SP_M1;

%% Practical Salinity: fill in missing values

% There is a 2-year gap for 3/5 of levels in this time series
SP_IC1_f = fillMissingSalinity(5,SP_IC1',time,18000);

% There is a 2-year gap for 1/5 of the levels in this time series
SP_IC2_f = fillMissingSalinity(5,SP_IC2',time,18000);

% There are no major gaps in the following time series
SP_IC3_f = fillMissingSalinity(4,SP_IC3',time);
SP_IC4_f = fillMissingSalinity(4,SP_IC4',time);
SP_M1_f = fillMissingSalinity(10,SP_M1',time);

save Matfiles/Sp_filled.mat SP_IC1_f SP_IC2_f SP_IC3_f SP_IC4_f SP_M1_f;

%% Depths: save for later

save Matfiles/depthsS depths_IC1 depths_IC2 depths_IC3 depths_IC4 depths_M1;

%% Pressures: needed for Sp->SA conversion

p_IC1 = sw_pres(depths_IC1,latitude(2))';
p_IC2 = sw_pres(depths_IC2,latitude(3))';
p_IC3 = sw_pres(depths_IC3,latitude(4))';
p_IC4 = sw_pres(depths_IC4,latitude(5))';
p_M1 = sw_pres(depths_M1,latitude(6))';

save Matfiles/p_S.mat p_IC1 p_IC2 p_IC3 p_IC4 p_M1;

%% Convert Practical Salinity (SP) into Absolute Salinity (SA)
% Note that the convertIntoAbsoluteSalinity function will save a mat file
% according to the filename entered.

% SA must be found before CT because CT is dependent on SA.

%% SP1->SA1
if isfile('SA_IC1.mat')
    disp('SA already calculated for IC1');
    load SA_IC1.mat;
else
    disp('Calculating SA for IC1');
    [SA_IC1] = convertIntoAbsoluteSalinity(SP_IC1_f,time,p_IC1,lon(2),latitude(2),"Matfiles/SA_IC1");
end

%% SP2->SA2
if isfile('SA_IC2.mat')
    disp('SA already calculated for IC2');
    load SA_IC2.mat;
else
    disp('Calculating SA for IC2');
    SA_IC2 = convertIntoAbsoluteSalinity(SP_IC2_f,time,p_IC2,lon(3),latitude(3),"Matfiles/SA_IC2");
end

%% SP3->SA3
if isfile('SA_IC3.mat')
    disp('SA already calculated for IC3');
    load SA_IC3.mat;
else
    disp('Calculating SA for IC3');
    SA_IC3 = convertIntoAbsoluteSalinity(SP_IC3_f,time,p_IC3,lon(4),latitude(4),"Matfiles/SA_IC3");
end

%% SP4->SA4
if isfile('SA_IC4.mat')
    disp('SA already calculated for IC4');
    load SA_IC4.mat;
else
    disp('Calculating SA for IC4');
    SA_IC4 = convertIntoAbsoluteSalinity(SP_IC4_f,time,p_IC4,lon(5),latitude(5),"Matfiles/SA_IC4");
end

%% SP-M1 -> SA-M1
if isfile('SA_M1.mat')
    disp('SA already calculated for M1');
    load SA_M1.mat;
else
    disp('Calculating SA for M1');
    SA_M1 = convertIntoAbsoluteSalinity(SP_M1_f,time,p_M1,lon(6),latitude(6),"Matfiles/SA_M1");
end
