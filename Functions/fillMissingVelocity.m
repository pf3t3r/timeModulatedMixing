function [velocityFilled] = fillMissingVelocity(noOfSensors,velocityWithGaps,time,window)
% fillMissingVelocity
% This function takes a velocity time series with missing values as an
% input and outputs a velocity time series where the gaps have been
% filled by a moving mean with a certain window width

% Default window width = 588; This gives a good fit to the two week period
% of missing data typical of redeployment
if (~exist('window', 'var'))
        window = 588;
end

velocityWithGaps(velocityWithGaps==0) = NaN;

velocityFilled = [];

for i = 1:noOfSensors
    velocityFilled(:,i) = fillmissing(velocityWithGaps(:,i),"movmean",2*window);
end

% test to see whether all NaN values have been filled
velocityFilled(isnan(velocityFilled)) = 0;
zerosV = find(~velocityFilled);
assert(isempty(zerosV));

figure
plot(time,velocityWithGaps(:,27),'DisplayName','input velocity','LineStyle','-');
hold on
plot(time,velocityFilled(:,27),'DisplayName','filled velocity','LineStyle',':');
hold off
legend();
end