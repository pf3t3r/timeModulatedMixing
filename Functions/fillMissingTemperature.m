function [temperatureFilled] = fillMissingTemperature(noOfSensors,temperatureWithGaps,time,window)
%filledSalinityTimeSeries takes as input an array of salinity values that
%are missing some values and also the number of sensors across the water
%column. It then fills the values using MATLAB's fillmissing function.

if (~exist('window', 'var'))
        window = 588;
end

temperatureProcessing = temperatureWithGaps;

% If the array starts with a NaN reassign it as a mean of T at that depth
% over a two-month period.
for i=1:noOfSensors
    if isnan(temperatureWithGaps(1,i))
        temperatureProcessing(1,i) = mean(temperatureProcessing(1:24*60,i),'omitnan');
    end
end

% Replace unphysical values of T (defined as < 1K) with a mean value
% for T at that depth.
for i=1:length(time)
    for j=1:noOfSensors
        if temperatureWithGaps(i,j) < 1
            temperatureProcessing(i,j) = mean(temperatureProcessing(:,j),'omitnan');
        end
    end
end

% Fill missing values according to [double] the minimum window size needed
% to bridge the gap. This factor should smooth the value a little...
temperatureFilled = fillmissing(temperatureProcessing,"movmean",2*window);

% test to see whether all NaN values have been filled
temperatureFilled(isnan(temperatureFilled)) = 0;

figure
plot(time,temperatureWithGaps(:,2),'DisplayName','raw');
hold on
plot(time,temperatureFilled(:,2),':','DisplayName','filled');
hold off
legend();
xlabel('time');
ylabel('T [g/kg]');
title('filling in missing T values');

% Ensure that there are no missing values in the new Sp array
zerosT = find(~temperatureFilled);
assert(isempty(zerosT));

temperatureFilled = temperatureFilled';
end