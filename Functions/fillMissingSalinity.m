function [salinityFilled] = fillMissingSalinity(noOfSensors,salinityWithGaps,time,window)
% fillMissingSalinity
% This function takes as input an array of salinity values that
% are missing some values and also the number of sensors across the water
% column. It then fills the values using MATLAB's fillmissing function.

if (~exist('window', 'var'))
        window = 588;
end

salinityProcessing = salinityWithGaps;

% If the array starts with a NaN reassign it as a mean of SP at that depth
% over a two-month period.
for i=1:noOfSensors
    if isnan(salinityWithGaps(1,i))
        salinityProcessing(1,i) = mean(salinityProcessing(1:24*60,i),'omitnan');
    end
end

% Replace unphysical values of SP (defined as < 25g/kg) with a mean value
% for SA at that depth.
for i=1:length(time)
    for j=1:noOfSensors
        if salinityWithGaps(i,j) < 25
            salinityProcessing(i,j) = mean(salinityProcessing(:,j),'omitnan');
        end
    end
end

% Fill missing values according to [double] the minimum window size needed
% to bridge the gap. This factor should smooth the value a little...
salinityFilled = fillmissing(salinityProcessing,"movmean",2*window);

% test to see whether all NaN values have been filled
salinityFilled(isnan(salinityFilled)) = 0;

figure
plot(time,salinityWithGaps(:,3),'DisplayName','raw');
hold on
plot(time,salinityFilled(:,3),':','DisplayName','filled');
hold off
legend();
xlabel('time');
ylabel('S_p [g/kg]');
title('filling in missing S_p values');

% Ensure that there are no missing values in the new Sp array
zerosSP = find(~salinityFilled);
assert(isempty(zerosSP));

salinityFilled = salinityFilled';
end